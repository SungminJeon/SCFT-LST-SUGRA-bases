// sugra_generator.h
// SUGRA base construction from LST catalog
//
// Pipeline:
//   1. Load LST catalog → reconstruct IF (via theory_sugra.h)
//   2. Attach external curves (hardcoded rules)
//   3. Optionally glue dummy LSTs (A-type: 1222...221)
//   4. NHC check: all non-(-1) clusters must be recognized NHCs
//   5. H, V from NHC table → gravitational anomaly: H - V + 29T - 273 = 0

#pragma once

#include "theory_sugra.h"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>  // std::gcd
#include <unordered_set>
#include <cstdint>
#include <string_view>

// ############################################################################
// PART 0: 128-bit dedup key (binary fingerprint + double hash)
// ############################################################################
// Two independent 64-bit hashes (FNV-1a forward and std::hash<string_view>)
// give effective ~128-bit collision resistance: P(collision) ≈ N²/2^128, which
// is ~10⁻²² for N=10⁹ — safely below numerical noise at any catalog size we
// will run.

struct DedupKey128 {
    uint64_t lo, hi;
    bool operator==(const DedupKey128& o) const noexcept {
        return lo == o.lo && hi == o.hi;
    }
};

struct DedupKey128Hash {
    size_t operator()(const DedupKey128& k) const noexcept {
        // mix lo and hi; both halves already well-distributed
        return (size_t)(k.lo ^ (k.hi * 0x9E3779B97F4A7C15ULL));
    }
};

inline uint64_t fnv1a_64(const uint8_t* data, size_t len) {
    const uint64_t basis = 0xcbf29ce484222325ULL;
    const uint64_t prime = 0x100000001b3ULL;
    uint64_t h = basis;
    for (size_t i = 0; i < len; i++) {
        h ^= data[i];
        h *= prime;
    }
    return h;
}

// Two-hash combination: FNV-1a forward + std::hash<string_view> (typically
// MurmurHash2/3 family). Different hash families → independent enough for
// practical 128-bit collision resistance.
// Second 64-bit half of the 128-bit dedup key. Uses FNV-1a with a different
// offset basis from the first half so the two hashes are statistically
// independent. std::hash<string_view> on some platforms (e.g. libc++'s
// __murmur2_or_cityhash with process-randomized seed) is not guaranteed to be
// run-to-run deterministic; using a pure FNV-1a here guarantees that
// dedup_key_v2_binary returns the same key for the same buffer across
// processes and across OpenMP runs.
inline uint64_t fnv1a_64_alt(const uint8_t* data, size_t len) {
    // Alternate basis (~ FNV-1a basis XOR'd with a 64-bit constant) to make the
    // second hash effectively independent of the first.
    const uint64_t basis = 0xcbf29ce484222325ULL ^ 0xa5a5a5a5a5a5a5a5ULL;
    const uint64_t prime = 0x100000001b3ULL;
    uint64_t h = basis;
    for (size_t i = 0; i < len; i++) {
        h ^= data[i];
        h *= prime;
    }
    return h;
}

inline DedupKey128 hash128(const std::vector<uint8_t>& buf) {
    DedupKey128 k;
    k.lo = fnv1a_64(buf.data(), buf.size());
    k.hi = fnv1a_64_alt(buf.data(), buf.size());
    return k;
}

inline void buf_append_i32(std::vector<uint8_t>& buf, int32_t v) {
    for (int i = 0; i < 4; i++) buf.push_back((uint8_t)((v >> (i*8)) & 0xff));
}
inline void buf_append_i64(std::vector<uint8_t>& buf, int64_t v) {
    for (int i = 0; i < 8; i++) buf.push_back((uint8_t)((v >> (i*8)) & 0xff));
}
inline void buf_append_bytes(std::vector<uint8_t>& buf, const void* data, size_t len) {
    const uint8_t* p = (const uint8_t*)data;
    buf.insert(buf.end(), p, p + len);
}

// ############################################################################
// PART 1: GAUGE DATA + CENTRAL CHARGE
// ############################################################################

struct GaugeInfo {
    std::string name;
    int dim;
    int h_dual;
    int rank;
};

inline const GaugeInfo GAUGE_NONE = {"none", 0, 0, 0};
inline const GaugeInfo GAUGE_SU2 = {"su2",  3, 2, 1};
inline const GaugeInfo GAUGE_SP1 = {"sp1",  3, 2, 1};
inline const GaugeInfo GAUGE_SU3 = {"su3",  8, 3, 2};
inline const GaugeInfo GAUGE_G2  = {"g2",  14, 4, 2};
inline const GaugeInfo GAUGE_SO7 = {"so7", 21, 5, 3};
inline const GaugeInfo GAUGE_SO8 = {"so8", 28, 6, 4};
inline const GaugeInfo GAUGE_F4  = {"f4",  52, 9, 4};
inline const GaugeInfo GAUGE_E6  = {"e6",  78, 12, 6};
inline const GaugeInfo GAUGE_E7  = {"e7", 133, 18, 7};
inline const GaugeInfo GAUGE_E8  = {"e8", 248, 30, 8};
inline const GaugeInfo GAUGE_SO16 = {"so16", 120, 14, 8};
inline const GaugeInfo GAUGE_SU8 = {"su8", 63, 8, 7};
inline const GaugeInfo GAUGE_SU16 = {"su16", 255, 16, 15};

inline double central_charge(const GaugeInfo& g, int k) {
    if (g.dim == 0) return 0.0;
    return (double)(k * g.dim) / (double)(k + g.h_dual);
}

inline GaugeInfo gauge_from_si(int si) {
    switch (si) {
        case -3:  return GAUGE_SU3;
        case -4:  return GAUGE_SO8;
        case -5:  return GAUGE_F4;
        case -6:  return GAUGE_E6;
        case -7: case -8: return GAUGE_E7;
        case -12: return GAUGE_E8;
        default:  return GAUGE_NONE;
    }
}

// ############################################################################
// PART 2: NHC TABLE
// ############################################################################
//
// Each NHC has:
//   - self_ints: the chain of self-intersections
//   - gauges: gauge algebra per curve
//   - H, V: total hypermultiplet / vector multiplet contribution
//   - spectrum: eigenvalues for matching (auto-computed)
//
// After removing -1 curves, every connected component of non-(-1)
// curves must match one of these NHCs.
// Total anomaly: sum(H_i) - sum(V_i) + 29T - 273 = 0

struct NHCDef {
    std::string name;
    std::vector<int> self_ints;
    std::vector<GaugeInfo> gauges;
    int H;                              // hypermultiplet contribution
    int V;                              // vector multiplet contribution
    std::vector<double> spectrum;       // auto-computed
};

inline std::vector<double> eigenvalue_spectrum(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return {};
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    std::vector<double> s(n);
    for (int i = 0; i < n; i++) s[i] = es.eigenvalues()(i);
    std::sort(s.begin(), s.end());
    return s;
}

inline bool spectra_equal(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-9) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++)
        if (std::abs(a[i] - b[i]) > tol) return false;
    return true;
}

inline const std::vector<NHCDef>& get_nhc_table() {
    static std::vector<NHCDef> table = []() {
        std::vector<NHCDef> t;
        
        auto add = [&](const std::string& name,
                       std::vector<int> si,
                       std::vector<GaugeInfo> g,
                       int H, int V,
                       std::vector<int> ints = {}) {
            NHCDef def;
            def.name = name;
            def.self_ints = si;
            def.gauges = g;
            def.H = H;
            def.V = V;
            int n = si.size();
            Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
            for (int i = 0; i < n; i++) {
                IF(i, i) = si[i];
                if (i + 1 < n) {
                    int k = (i < (int)ints.size()) ? ints[i] : 1;
                    IF(i, i+1) = k; IF(i+1, i) = k;
                }
            }
            def.spectrum = eigenvalue_spectrum(IF);
            t.push_back(def);
        };
        
        // ┌─────────────────────────────────────────────────────────┐
        // │ NHC Table: name, {self_ints}, {gauges}, H, V           │
        // │                                                         │
        // │ USER: fill in correct H, V for each NHC                │
        // │ H=0, V=0 are placeholders                              │
        // └─────────────────────────────────────────────────────────┘
        
        // Single curves
        add("-3",     {-3},         {GAUGE_SU3},                    0, 8);
        add("-4",     {-4},         {GAUGE_SO8},                    0, 28);
        add("-5",     {-5},         {GAUGE_F4},                     0, 52);
        add("-6",     {-6},         {GAUGE_E6},                     0, 78);
        add("-7",     {-7},         {GAUGE_E7},                     28, 133);
        add("-8",     {-8},         {GAUGE_E7},                     0, 133);
        add("-12",    {-12},        {GAUGE_E8},                     0, 248);
        
        // Multi-curve chains
        add("-2-3",   {-2, -3},     {GAUGE_SU2, GAUGE_G2},         8, 17);
        add("-2-3-2", {-2, -3, -2}, {GAUGE_SU2, GAUGE_SO7, GAUGE_SU2}, 16, 27);
        add("-2-2-3", {-2, -2, -3}, {GAUGE_NONE, GAUGE_SP1, GAUGE_G2}, 8, 17);
        
        // USER: add more NHC chains as needed
        // add("name", {self_ints...}, {gauges...}, H, V);
	add("-2-4", {-2, -4}, {GAUGE_SU8, GAUGE_SO16},             128,183, {2});

        
        return t;
    }();
    return table;
}

// ############################################################################
// PART 3: NHC CHECK → CLUSTER ASSIGNMENT
// ############################################################################

struct ClusterInfo {
    int cluster_id;
    std::vector<int> curve_indices;
    const NHCDef* nhc;             // nullptr for pure -2 clusters
    int H, V;
};

struct NHCResult {
    bool passes;
    std::string fail_reason;
    std::vector<ClusterInfo> clusters;
    
    // Per-curve info
    std::vector<GaugeInfo> curve_gauges;    // indexed by curve in IF
    std::vector<int> curve_cluster;         // cluster_id per curve (-1 for -1 curves)
    
    int total_H() const {
        int h = 0;
        for (auto& c : clusters) h += c.H;
        return h;
    }
    int total_V() const {
        int v = 0;
        for (auto& c : clusters) v += c.V;
        return v;
    }
};

inline std::vector<int> walk_chain(const Eigen::MatrixXi& IF, const std::vector<int>& cluster) {
    int m = cluster.size();
    if (m <= 1) return cluster;
    std::map<int, std::vector<int>> adj;
    for (int a : cluster)
        for (int b : cluster)
            if (a != b && IF(a, b) != 0) adj[a].push_back(b);
    int start = cluster[0];
    for (int c : cluster) if ((int)adj[c].size() <= 1) { start = c; break; }
    std::vector<int> ordered;
    std::set<int> vis;
    int cur = start;
    while (true) {
        ordered.push_back(cur);
        vis.insert(cur);
        bool found = false;
        for (int nb : adj[cur])
            if (vis.find(nb) == vis.end()) { cur = nb; found = true; break; }
        if (!found) break;
    }
    return ordered;
}

inline NHCResult check_nhc(const Eigen::MatrixXi& IF, double cc_budget = 8.0) {
    NHCResult result;
    result.passes = true;
    int n = IF.rows();
    if (n == 0) return result;
    
    result.curve_gauges.assign(n, GAUGE_NONE);
    result.curve_cluster.assign(n, -1);
    
    // Separate -1 and non-(-1) curves
    std::vector<int> minus1, non_minus1;
    for (int i = 0; i < n; i++) {
        if (IF(i, i) == -1) minus1.push_back(i);
        else non_minus1.push_back(i);
    }
    if (non_minus1.empty()) return result;
    
    // BFS: connected components of non-(-1) curves
    std::vector<std::vector<int>> comps;
    std::vector<bool> visited(n, false);
    for (int start : non_minus1) {
        if (visited[start]) continue;
        std::vector<int> comp;
        std::vector<int> queue = {start};
        visited[start] = true;
        while (!queue.empty()) {
            int cur = queue.back(); queue.pop_back();
            comp.push_back(cur);
            for (int j : non_minus1)
                if (!visited[j] && IF(cur, j) != 0) { visited[j] = true; queue.push_back(j); }
        }
        std::sort(comp.begin(), comp.end());
        comps.push_back(comp);
    }
    
    // Match each component to NHC table
    for (int ci = 0; ci < (int)comps.size(); ci++) {
        auto& comp = comps[ci];
        int m = comp.size();
        
        // Extract submatrix
        Eigen::MatrixXi sub(m, m);
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                sub(a, b) = IF(comp[a], comp[b]);
        
        auto spec = eigenvalue_spectrum(sub);
        
        // Check pure -2
        bool pure2 = true;
        for (int c : comp) if (IF(c, c) != -2) { pure2 = false; break; }
        
        if (pure2) {
            ClusterInfo cl;
            cl.cluster_id = ci;
            cl.curve_indices = comp;
            cl.nhc = nullptr;
            cl.H = 0;
            cl.V = 0;
            result.clusters.push_back(cl);
            for (int c : comp) result.curve_cluster[c] = ci;
            continue;
        }
        
        // NHC lookup
        const NHCDef* nhc = nullptr;
        for (const auto& def : get_nhc_table())
            if (spectra_equal(spec, def.spectrum)) { nhc = &def; break; }
        
        if (!nhc) {
            result.passes = false;
            std::ostringstream ss;
            ss << "L1: unknown cluster {";
            for (int a = 0; a < m; a++) { if (a) ss << ","; ss << IF(comp[a], comp[a]); }
            ss << "}";
            result.fail_reason = ss.str();
            return result;
        }
        
        // Assign gauges per curve
        auto ordered = walk_chain(IF, comp);
        std::vector<int> chain_si;
        for (int c : ordered) chain_si.push_back(IF(c, c));
        
        bool forward = (chain_si == nhc->self_ints);
        std::vector<int> rev(chain_si.rbegin(), chain_si.rend());
        bool reverse = (!forward && rev == nhc->self_ints);
        
        if (forward) {
            for (size_t i = 0; i < ordered.size(); i++) {
                result.curve_gauges[ordered[i]] = nhc->gauges[i];
                result.curve_cluster[ordered[i]] = ci;
            }
        } else if (reverse) {
            int sz = ordered.size();
            for (size_t i = 0; i < ordered.size(); i++) {
                result.curve_gauges[ordered[i]] = nhc->gauges[sz - 1 - i];
                result.curve_cluster[ordered[i]] = ci;
            }
        } else {
            // Fallback: assign by self-int
            for (int c : comp) {
                result.curve_gauges[c] = gauge_from_si(IF(c, c));
                result.curve_cluster[c] = ci;
            }
        }
        
        ClusterInfo cl;
        cl.cluster_id = ci;
        cl.curve_indices = comp;
        cl.nhc = nhc;
        cl.H = nhc->H;
        cl.V = nhc->V;
        result.clusters.push_back(cl);
    }
    
    // Central charge check on each -1 curve
    // Include (-1) neighbors that already have gauge (e.g. hat1 with su(8))
    for (int m1 : minus1) {
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == m1 || IF(m1, j) == 0) continue;
            if (IF(j, j) == -1 && result.curve_gauges[j].dim == 0) continue;
            int k = std::abs(IF(m1, j));
            total_c += central_charge(result.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            result.passes = false;
            std::ostringstream ss;
            ss << "L2: c.c.=" << std::fixed << std::setprecision(2)
               << total_c << " > " << cc_budget << " at -1 curve " << m1;
            result.fail_reason = ss.str();
            return result;
        }
    }
    
    return result;
}

// ── External -2 → -1 gauge enhancement ──
// When an external -2 curve is attached to a -1 curve, it carries su(2)
// gauge algebra: V=3, H=4 (4 half-hypers in fundamental).
// This must be applied AFTER check_nhc and BEFORE compute_anomaly.
inline void enhance_external_m2_gauge(NHCResult& nhc,
                                       const Eigen::MatrixXi& IF,
                                       int ext_idx, int target_si,
                                       double cc_budget = 8.0) {
    if (IF(ext_idx, ext_idx) != -2 || target_si != -1) return;
    // Skip if this -2 already carries gauge from check_nhc (part of an NHC chain
    // like -2-3, -2-3-2, -2-2-3, -2-4). Enhancement only applies to a standalone
    // -2 external attached purely to (-1) target.
    if (ext_idx >= 0 && ext_idx < (int)nhc.curve_gauges.size()
        && nhc.curve_gauges[ext_idx].dim != 0) return;

    // Set gauge on the external -2
    nhc.curve_gauges[ext_idx] = GAUGE_SU2;
    
    // Update cluster H, V
    for (auto& cl : nhc.clusters) {
        for (int c : cl.curve_indices) {
            if (c == ext_idx) {
                cl.V += 3;
                cl.H += 4;
                goto done_cluster;
            }
        }
    }
    done_cluster:
    
    // Re-check central charge on neighboring -1 curves
    // Include (-1) neighbors that have gauge (e.g. hat1 with su(8))
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (IF(i, i) != -1 || IF(i, ext_idx) == 0) continue;
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i || IF(i, j) == 0) continue;
            if (IF(j, j) == -1 && nhc.curve_gauges[j].dim == 0) continue;
            int k = std::abs(IF(i, j));
            total_c += central_charge(nhc.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            nhc.passes = false;
            nhc.fail_reason = "L2(ext): c.c. exceeded after su(2) enhancement";
            return;
        }
    }
}

// ── hat{1} gauge enhancement ──
// hat{1} (C²=-1, b₀·C=-1) carries gauge algebra depending on target:
//   → -1 target: su(8),  V=63,  H=36  (1 symmetric hyper)
//   → -2 target: su(16), V=255, H=264 (1 symmetric + 8 fundamental hypers)
// This must be applied AFTER check_nhc and BEFORE compute_anomaly.
inline void enhance_hat1_gauge(NHCResult& nhc,
                                const Eigen::MatrixXi& IF,
                                int hat1_idx, int target_si,
                                double cc_budget = 8.0) {
    if (IF(hat1_idx, hat1_idx) != -1) return;
    
    GaugeInfo gauge;
    int dV, dH;
    if (target_si == -1) {
        gauge = GAUGE_SU8; dV = 63; dH = 36;
    } else if (target_si == -2) {
        gauge = GAUGE_SU16; dV = 255; dH = 264;
    } else {
        return;
    }
    
    nhc.curve_gauges[hat1_idx] = gauge;
    
    // Update cluster H, V
    for (auto& cl : nhc.clusters) {
        for (int c : cl.curve_indices) {
            if (c == hat1_idx) {
                cl.V += dV;
                cl.H += dH;
                goto done_hat1;
            }
        }
    }
    // hat1 is a -1 curve, not in any cluster — add a new one
    {
        ClusterInfo cl;
        cl.cluster_id = (int)nhc.clusters.size();
        cl.curve_indices = {hat1_idx};
        cl.nhc = nullptr;
        cl.H = dH;
        cl.V = dV;
        nhc.clusters.push_back(cl);
        nhc.curve_cluster[hat1_idx] = cl.cluster_id;
    }
    done_hat1:
    
    // Re-check central charge on neighboring -1 curves
    // Include all (-1) neighbors that have gauge (hat1 with su(8), etc.)
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (i == hat1_idx || IF(i, i) != -1 || IF(i, hat1_idx) == 0) continue;
        double total_c = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i || IF(i, j) == 0) continue;
            if (IF(j, j) == -1 && nhc.curve_gauges[j].dim == 0) continue;
            int k = std::abs(IF(i, j));
            total_c += central_charge(nhc.curve_gauges[j], k);
        }
        if (total_c > cc_budget + 1e-9) {
            nhc.passes = false;
            nhc.fail_reason = "L2(hat1): c.c. exceeded after gauge enhancement";
            return;
        }
    }
}

// ############################################################################
// PART 4: SIGNATURE + DETERMINANT
// ############################################################################

struct SigInfo {
    int sig_pos, sig_neg, sig_zero;
    long long det;
};

// Fast signature check using LDLT decomposition (no eigenvalues computed).
// Much faster than full eigensolve when only sig_pos/sig_neg/sig_zero are needed.
// Does NOT compute det or eigenvalues — use compute_sig for those.
// Fast signature pre-filter using LDLT decomposition.
// Uses a conservative (small) threshold so that borderline eigenvalues
// are counted as positive/negative rather than zero. This ensures we
// never falsely reject (sig_pos may be overestimated, never underestimated).
// Returns true if sig_pos > 1 is definite (safe to reject).
inline bool sig_pos_exceeds_one_fast(const Eigen::MatrixXi& IF, int T_max) {
    int n = IF.rows();
    if (n == 0) return false;
    Eigen::MatrixXd Fd = IF.cast<double>();
    // Hybrid: fast LDLT for the common case, eigensolve fallback when LDLT's pivots
    // are near zero (catastrophic cancellation → unreliable D signs). Threshold 1e-6
    // is well above the catastrophic regime (~1e-16) but well below typical valid |D|.
    Eigen::LDLT<Eigen::MatrixXd> ldlt(Fd);
    auto D = ldlt.vectorD();
    double min_abs = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        double a = std::abs(D(i));
        if (a < min_abs) min_abs = a;
    }
    if (min_abs < 1e-6) {
        // Unreliable LDLT → fall back to permutation-invariant eigensolve
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.compute(Fd, Eigen::EigenvaluesOnly);
        int pos = 0, neg = 0;
        for (int i = 0; i < n; i++) {
            double e = es.eigenvalues()(i);
            if (e > 1e-8) pos++;
            else if (e < -1e-8) neg++;
        }
        return (pos > 1 || neg > T_max);
    }
    int pos = 0, neg = 0;
    for (int i = 0; i < n; i++) {
        if (D(i) > 1e-12) pos++;
        else if (D(i) < -1e-12) neg++;
    }
    return (pos > 1 || neg > T_max);
}

// Check extended IF (with b0) signature.
// For valid SUGRA base: extended IF must have sig_pos == 1.
// h1_set: indices of hat1 curves (b0·C = -1 instead of C²+2)
inline bool check_extended_sig(const Eigen::MatrixXi& IF, int T,
                                const std::set<int>& h1_set = {}) {
    int n = IF.rows();
    Eigen::MatrixXi ext = Eigen::MatrixXi::Zero(n+1, n+1);
    ext.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        int b = h1_set.count(i) ? -1 : (IF(i, i) + 2);
        ext(i, n) = b;
        ext(n, i) = b;
    }
    ext(n, n) = 9 - T;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ext.cast<double>());
    int pos = 0;
    for (int i = 0; i < n + 1; i++)
        if (es.eigenvalues()(i) > 1e-8) pos++;
    return pos == 1;
}

// Sentinel returned by exact_det when the true determinant does not fit in the
// supported range. It is >1 in magnitude AND prime (never a perfect square), so
// the |det|<=1 and is_perfect_square filters both conservatively reject it.
static constexpr long long DET_TOOLARGE = 9223372036854775783LL; // largest prime < 2^63

// Exact integer determinant via fraction-free (Bareiss) Gaussian elimination.
// Intersection forms are exact integer matrices, so this avoids the precision
// loss of Fd.determinant() (double, inexact past 2^53) and the int truncation
// of casting it back. Intermediate minors are bounded by Hadamard's inequality
// and grow with n; if any minor would exceed the __int128 budget we bail to
// DET_TOOLARGE (such a |det| is astronomically large — neither unimodular nor a
// small perfect square, so the conservative reject is physically correct).
inline long long exact_det_impl(const Eigen::MatrixXi& IF) {
    int n = IF.rows();
    if (n == 0) return 1;  // det of empty matrix = 1 (unused, but well-defined)
    std::vector<__int128> a((size_t)n * n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            a[(size_t)i * n + j] = (__int128)IF(i, j);
    auto at = [&](int i, int j) -> __int128& { return a[(size_t)i * n + j]; };
    // Keep every minor within ±GUARD so that products of two minors stay inside
    // __int128 (2^62 * 2^62 = 2^124 < 2^127). GUARD < 2^63 ⇒ final det fits int64.
    const __int128 GUARD = ((__int128)1) << 62;
    __int128 prev = 1;
    int sign = 1;
    for (int k = 0; k < n - 1; k++) {
        if (at(k, k) == 0) {
            int p = -1;
            for (int i = k + 1; i < n; i++) if (at(i, k) != 0) { p = i; break; }
            if (p < 0) return 0;  // zero column at/below pivot ⇒ singular ⇒ det 0
            for (int j = 0; j < n; j++) std::swap(at(k, j), at(p, j));
            sign = -sign;
        }
        __int128 piv = at(k, k);
        for (int i = k + 1; i < n; i++) {
            __int128 aik = at(i, k);
            for (int j = k + 1; j < n; j++) {
                __int128 v = (at(i, j) * piv - aik * at(k, j)) / prev;  // Bareiss: exact
                if (v > GUARD || v < -GUARD) return DET_TOOLARGE;
                at(i, j) = v;
            }
            at(i, k) = 0;
        }
        prev = piv;
    }
    return (long long)((__int128)sign * at(n - 1, n - 1));
}

// Public exact determinant (currently unused in the hot path — compute_sig uses
// det_cap(∏spectrum); kept for a future lazy-exact --det-sq path).
inline long long exact_det(const Eigen::MatrixXi& IF) { return exact_det_impl(IF); }

inline bool is_perfect_square(long long n) {
    if (n < 0) return false;
    if (n == 0) return true;  // det=0: valid (fiber class exists)
    long long s = (long long)std::llround(std::sqrt((double)n));
    // Check neighbours to absorb any double rounding near the boundary.
    for (long long t = (s > 0 ? s - 1 : 0); t <= s + 1; t++)
        if ((__int128)t * t == (__int128)n) return true;
    return false;
}

// Capped determinant from a floating product of the spectrum (eigenvalues, or
// LDLT diagonal D since det(A)=∏D_i). O(n) — replaces the O(n^3) __int128 Bareiss
// in the hot path. Same double precision as the legacy Fd.determinant(); values
// beyond double's exact-integer range (~2^52) are not needed by any filter
// (|det|<=1 fails; perfect-square only matters for small dets) → DET_TOOLARGE.
inline long long det_cap(double pdet) {
    return (std::abs(pdet) > 4.5e15) ? DET_TOOLARGE : (long long)std::llround(pdet);
}

inline SigInfo compute_sig(const Eigen::MatrixXi& IF) {
    SigInfo s = {0, 0, 0, 0};
    int n = IF.rows();
    if (n == 0) return s;
    Eigen::MatrixXd Fd = IF.cast<double>();
    // EigenvaluesOnly: ~2× faster than full diagonalization (we don't need eigenvectors)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(Fd, Eigen::EigenvaluesOnly);
    double pdet = 1.0;
    for (int i = 0; i < n; i++) {
        double ev = es.eigenvalues()(i);
        pdet *= ev;
        if (std::abs(ev) < 1e-8) s.sig_zero++;
        else if (ev > 0) s.sig_pos++;
        else s.sig_neg++;
    }
    s.det = det_cap(pdet);
    return s;
}

// Hybrid signature: LDLT fast path + eigensolve fallback.
// Sylvester's law: D vector from LDLT gives sig_pos/sig_neg/sig_zero counts
// directly (matches sig_pos_exceeds_one_fast's reliability gate). For n=193
// this avoids the per-leaf eigensolve in the common (reliable LDLT) case.
inline SigInfo compute_sig_fast(const Eigen::MatrixXi& IF) {
    SigInfo s = {0, 0, 0, 0};
    int n = IF.rows();
    if (n == 0) return s;
    Eigen::MatrixXd Fd = IF.cast<double>();

    Eigen::LDLT<Eigen::MatrixXd> ldlt(Fd);
    auto D = ldlt.vectorD();
    double min_abs = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++) {
        double a = std::abs(D(i));
        if (a < min_abs) min_abs = a;
    }

    if (min_abs >= 1e-6) {
        // LDLT reliable — extract signature + det.
        double pdet = 1.0;
        for (int i = 0; i < n; i++) {
            double d = D(i);
            pdet *= d;  // det(A) = ∏D_i (permutation determinant squares to 1)
            if (std::abs(d) < 1e-8) s.sig_zero++;
            else if (d > 0) s.sig_pos++;
            else s.sig_neg++;
        }
        s.det = det_cap(pdet);
        return s;
    }

    // Fallback: eigensolve when LDLT pivots are near zero (catastrophic
    // cancellation in D signs is possible). Threshold 1e-6 matches the
    // existing sig_pos_exceeds_one_fast gate.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(Fd, Eigen::EigenvaluesOnly);
    double pdet = 1.0;
    for (int i = 0; i < n; i++) {
        double ev = es.eigenvalues()(i);
        pdet *= ev;
        if (std::abs(ev) < 1e-8) s.sig_zero++;
        else if (ev > 0) s.sig_pos++;
        else s.sig_neg++;
    }
    s.det = det_cap(pdet);
    return s;
}

// Overload: also output eigenvalues for reuse (e.g. in dedup_key_v2)
inline SigInfo compute_sig(const Eigen::MatrixXi& IF, std::vector<double>& eigenvalues_out) {
    SigInfo s = {0, 0, 0, 0};
    int n = IF.rows();
    eigenvalues_out.clear();
    if (n == 0) return s;
    Eigen::MatrixXd Fd = IF.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(Fd, Eigen::EigenvaluesOnly);
    eigenvalues_out.resize(n);
    double pdet = 1.0;
    for (int i = 0; i < n; i++) {
        double ev = es.eigenvalues()(i);
        eigenvalues_out[i] = ev;
        pdet *= ev;
        if (std::abs(ev) < 1e-8) s.sig_zero++;
        else if (ev > 0) s.sig_pos++;
        else s.sig_neg++;
    }
    s.det = det_cap(pdet);
    return s;
}

// ############################################################################
// PART 4.5: PROPER c_ext (NHC-aware, level-based)
// ############################################################################
//
// For each EXTERNAL gauge-carrying curve i (i ∈ ext_indices and
// curve_gauges[i] non-trivial):
//   k_i = sum of |IF(i,j)| where j ∉ ext_indices (intersection to LST side)
//   c_i = central_charge(curve_gauges[i], k_i) = k_i·dim/(k_i+h)
// c_ext = Σ_i c_i.
//
// We do NOT count LST curves enhanced into a higher gauge by ext attachment —
// only the curves that are themselves external carry the c_ext charge.
// Filter: reject if c_ext > 16.

inline double compute_proper_c_ext(const Eigen::MatrixXi& IF,
                                    const NHCResult& nhc,
                                    const std::set<int>& ext_indices) {
    int n = IF.rows();
    double total = 0.0;
    for (int i : ext_indices) {
        if (i < 0 || i >= n) continue;
        if (i >= (int)nhc.curve_gauges.size()) continue;
        if (nhc.curve_gauges[i].dim == 0) continue;
        int k = 0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            if (ext_indices.count(j) > 0) continue;   // skip intra-ext bonds
            k += IF(i, j);
        }
        if (k > 0) total += central_charge(nhc.curve_gauges[i], k);
    }
    return total;
}

static double g_proper_c_ext_max = 16.0;  // proper NHC-aware c_ext cap

// ----------------------------------------------------------------------------
// Fiber-weighted true c_ext (final SUGRA filter)
// ----------------------------------------------------------------------------
// Find positive integer vector X in ker(A), normalized so gcd(X) = 1.
// Returns empty vector if no such X exists or A has no usable null space.
inline Eigen::VectorXi find_positive_integer_kernel(const Eigen::MatrixXi& A_int) {
    int n = A_int.rows();
    if (n == 0) return Eigen::VectorXi();
    Eigen::MatrixXd A = A_int.cast<double>();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
    auto s = svd.singularValues();
    auto V = svd.matrixV();
    int n_s = (int)s.size();
    if (n_s == 0) return Eigen::VectorXi();
    if (s(n_s - 1) >= 1e-8) return Eigen::VectorXi();  // no null space

    Eigen::VectorXd v = V.col(n_s - 1);
    // Flip sign if needed so all entries are non-negative (within tolerance)
    bool any_pos = false, any_neg = false;
    for (int i = 0; i < n; i++) {
        if (v(i) > 1e-9) any_pos = true;
        if (v(i) < -1e-9) any_neg = true;
    }
    if (any_neg && !any_pos) v = -v;
    else if (any_pos && any_neg) return Eigen::VectorXi();  // mixed signs

    for (int i = 0; i < n; i++) v(i) = std::abs(v(i));
    double min_abs = std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; i++)
        if (v(i) > 1e-9 && v(i) < min_abs) min_abs = v(i);
    if (!std::isfinite(min_abs)) return Eigen::VectorXi();

    Eigen::VectorXd scaled = v / min_abs;
    Eigen::VectorXi X(n);
    for (int i = 0; i < n; i++) {
        double rounded = std::round(scaled(i));
        if (std::abs(scaled(i) - rounded) > 1e-6) return Eigen::VectorXi();
        X(i) = (int)rounded;
        if (X(i) < 0) return Eigen::VectorXi();
    }
    int g = 0;
    for (int i = 0; i < n; i++) if (X(i) > 0) g = (g == 0) ? X(i) : std::gcd(g, X(i));
    if (g > 1) for (int i = 0; i < n; i++) X(i) /= g;
    return X;
}

// True c_ext using fiber-weighted k_eff = Σ IF(ext, i) * X[i] for i in LST.
// Returns -1.0 if kernel can't be determined (caller should not filter in that case).
inline double compute_true_c_ext_fiber(const Eigen::MatrixXi& IF,
                                        const NHCResult& nhc,
                                        const std::set<int>& ext_indices) {
    int n = IF.rows();
    std::vector<int> lst_idx;
    for (int i = 0; i < n; i++)
        if (ext_indices.count(i) == 0) lst_idx.push_back(i);
    int m = (int)lst_idx.size();
    if (m == 0) return 0.0;
    Eigen::MatrixXi A(m, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            A(i, j) = IF(lst_idx[i], lst_idx[j]);
    Eigen::VectorXi X = find_positive_integer_kernel(A);
    if (X.size() == 0) return -1.0;
    double total = 0.0;
    for (int e : ext_indices) {
        if (e < 0 || e >= n) continue;
        if (e >= (int)nhc.curve_gauges.size()) continue;
        if (nhc.curve_gauges[e].dim == 0) continue;
        int k_eff = 0;
        for (int i = 0; i < m; i++) k_eff += IF(e, lst_idx[i]) * X(i);
        if (k_eff <= 0) continue;
        total += central_charge(nhc.curve_gauges[e], k_eff);
    }
    return total;
}

// Cached variant: kernel X and lst_idx are pre-computed once per (LST, ext_set) pair.
// X.size() == 0 means kernel undefined → return -1.0 (caller skips filter).
// This skips the per-leaf SVD which is the dominant cost in enumeration.
inline double compute_true_c_ext_fiber_cached(const Eigen::MatrixXi& IF,
                                               const NHCResult& nhc,
                                               const std::set<int>& ext_indices,
                                               const std::vector<int>& lst_idx,
                                               const Eigen::VectorXi& X) {
    if (X.size() == 0) return -1.0;
    int n = IF.rows();
    int m = (int)lst_idx.size();
    double total = 0.0;
    for (int e : ext_indices) {
        if (e < 0 || e >= n) continue;
        if (e >= (int)nhc.curve_gauges.size()) continue;
        if (nhc.curve_gauges[e].dim == 0) continue;
        int k_eff = 0;
        for (int i = 0; i < m; i++) k_eff += IF(e, lst_idx[i]) * X(i);
        if (k_eff <= 0) continue;
        total += central_charge(nhc.curve_gauges[e], k_eff);
    }
    return total;
}

static double g_true_c_ext_max = 16.0;  // true (fiber-weighted) c_ext cap

// ############################################################################
// PART 4b: FAST INCREMENTAL SIGNATURE (Schur complement / Haynsworth)
// ############################################################################
// A base IF B (n×n, symmetric) is shared across MANY candidate attachments that
// differ only in the appended curves. Precompute B's inertia + pseudo-inverse
// ONCE; then each candidate M = [[B, A],[Aᵀ, D]] (k appended curves) gets its
// signature in O(n²k + (z+k)³) instead of a fresh O(n³) eigensolve/LDLT.
//
// Haynsworth inertia additivity: integrating out the invertible (non-null) part
// of B leaves a (z+k)×(z+k) effective form S:
//     In(M) = (p, q, 0)  +  In(S),
//     S = [[ NᵀBN = 0_z ,  NᵀA ],
//          [   AᵀN      ,  D − AᵀB⁺A ]]
// where p,q,z = #pos/#neg/#null eigenvalues of B, N = (n×z) null eigenvectors,
// B⁺ = pseudo-inverse of B on the non-null subspace. So:
//     sig_pos(M) = p + sig_pos(S),  sig_neg(M) = q + sig_neg(S),  sig_zero(M)=sig_zero(S).
//
// LST base (phase1 / nhc_ext): inertia (0, n−1, 1) → z=1 null = tensionless dir;
//   S is (1+k); a valid SUGRA attach must lift the null to exactly +1.
// Lorentzian SUGRA base (phase2 / phase3): inertia (1, T, 0) → z=0; S is plain
//   k×k (B⁻¹ regular); the appended curves must add NO new positive eigenvalue.
// Same code, only the base B differs.
struct BaseSchur {
    int n = 0;
    int p = 0, q = 0;          // # positive / negative eigenvalues of B
    Eigen::MatrixXd Nnull;     // n×z null eigenvectors
    Eigen::MatrixXd Bpinv;     // n×n pseudo-inverse of B on the non-null subspace
    double zero_tol = 0.0;     // threshold used to classify null eigenvalues of B
};

inline BaseSchur precompute_base_schur(const Eigen::MatrixXi& B) {
    BaseSchur bs;
    int n = B.rows();
    bs.n = n;
    if (n == 0) return bs;
    Eigen::MatrixXd Bd = B.cast<double>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Bd);
    const auto& lam = es.eigenvalues();
    const auto& U = es.eigenvectors();
    double scale = lam.cwiseAbs().maxCoeff();
    double tol = 1e-7 * std::max(scale, 1.0);
    bs.zero_tol = tol;
    bs.Bpinv = Eigen::MatrixXd::Zero(n, n);
    std::vector<int> null_cols;
    for (int i = 0; i < n; i++) {
        double l = lam(i);
        if (std::abs(l) <= tol) { null_cols.push_back(i); }
        else {
            if (l > 0) bs.p++; else bs.q++;
            bs.Bpinv.noalias() += (1.0 / l) * U.col(i) * U.col(i).transpose();
        }
    }
    bs.Nnull.resize(n, (int)null_cols.size());
    for (int j = 0; j < (int)null_cols.size(); j++) bs.Nnull.col(j) = U.col(null_cols[j]);
    return bs;
}

// Full signature of M (size (n+k)×(n+k)), whose top-left n×n block is the base B
// that produced bs. pos_tol classifies the (z+k)×(z+k) Schur eigenvalues.
inline SigInfo schur_sig(const BaseSchur& bs, const Eigen::MatrixXi& M, double pos_tol = 1e-8) {
    SigInfo s = {0, 0, 0, 0};
    int N = M.rows();
    int n = bs.n;
    int k = N - n;
    if (k < 0) return s;
    if (k == 0) { s.sig_pos = bs.p; s.sig_neg = bs.q; s.sig_zero = (int)bs.Nnull.cols(); return s; }
    int z = (int)bs.Nnull.cols();
    Eigen::MatrixXd A(n, k), D(k, k);
    for (int j = 0; j < k; j++) {
        for (int i = 0; i < n; i++) A(i, j) = (double)M(i, n + j);
        for (int i = 0; i < k; i++) D(i, j) = (double)M(n + i, n + j);
    }
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(z + k, z + k);
    S.block(z, z, k, k) = D - A.transpose() * (bs.Bpinv * A);  // cluster Schur complement
    if (z > 0) {
        Eigen::MatrixXd NtA = bs.Nnull.transpose() * A;        // z×k null coupling
        S.block(0, z, z, k) = NtA;
        S.block(z, 0, k, z) = NtA.transpose();
        // S.block(0,0,z,z) stays 0 (= NᵀBN)
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S, Eigen::EigenvaluesOnly);
    int sp = 0, sn = 0, sz = 0;
    for (int i = 0; i < z + k; i++) {
        double e = es.eigenvalues()(i);
        if (e > pos_tol) sp++; else if (e < -pos_tol) sn++; else sz++;
    }
    s.sig_pos = bs.p + sp;
    s.sig_neg = bs.q + sn;
    s.sig_zero = sz;
    return s;
}

// ############################################################################
// PART 5: ANOMALY
// ############################################################################

struct AnomalyResult {
    int T, H_charged, V;
    int H_neutral;      // 273 + V - H_charged - 29T >= 0
};

// LST-T anomaly mode: T = LST tensor count (number of curves in the LST sub-IF),
// set by the caller per input entry before invoking compute_anomaly. This is
// the canonical anomaly convention for the pipeline; the legacy "T = sig_neg"
// path has been removed.
static int g_current_base_T = 0;
#ifdef _OPENMP
// Each OpenMP thread mutates g_current_base_T independently per input entry.
#pragma omp threadprivate(g_current_base_T)
#endif

inline AnomalyResult compute_anomaly(const Eigen::MatrixXi& /*IF*/, const NHCResult& nhc) {
    AnomalyResult a;
    a.T = g_current_base_T;
    a.H_charged = nhc.total_H();
    a.V = nhc.total_V();
    a.H_neutral = 273 + a.V - a.H_charged - 29 * a.T;
    return a;
}

// Overload kept for API symmetry; `sig` is unused under the LST-T convention.
inline AnomalyResult compute_anomaly(const NHCResult& nhc, const SigInfo& /*sig*/) {
    AnomalyResult a;
    a.T = g_current_base_T;
    a.H_charged = nhc.total_H();
    a.V = nhc.total_V();
    a.H_neutral = 273 + a.V - a.H_charged - 29 * a.T;
    return a;
}


// Overload: use LST base_T instead of sig_neg (--use-lst-T mode)
inline AnomalyResult compute_anomaly_lstT(const NHCResult& nhc, const SigInfo& sig, int base_T) {
    AnomalyResult a;
    a.T = base_T;  // LST tensor number, not sig_neg
    a.H_charged = nhc.total_H();
    a.V = nhc.total_V();
    a.H_neutral = 273 + a.V - a.H_charged - 29 * a.T;
    return a;
}

// ############################################################################
// PART 6: ATTACHMENT RULES
// ############################################################################

struct AttachmentRule {
    int ext_self_int;
    int target_self_int;
    int intersection_num;
    std::vector<int> allowed_sidelinks;         // empty = all
    std::vector<std::string> allowed_nobles;    // empty = all
    std::vector<int> allowed_block_params;      // empty = all (for Nk)
};

// ┌───────────────────────────────────────────────────────────────┐
// │ USER: Fill allowed_sidelinks / allowed_nobles per rule        │
// │ Example:                                                      │
// │   {-2, -3, 1, {33, 331}, {"313"}, {}},                      │
// └───────────────────────────────────────────────────────────────┘

inline std::vector<AttachmentRule> build_attachment_rules() {
    std::vector<AttachmentRule> rules;
    
    rules.push_back({-2, -1, 1, {}, {}, {}});
    rules.push_back({-2, -1, 2, {}, {}, {}});
    rules.push_back({-2, -1, 3, {}, {}, {}});
    rules.push_back({-2, -1, 4, {}, {}, {}});
    rules.push_back({-2, -1, 5, {}, {}, {}});
    rules.push_back({-2, -1, 6, {}, {}, {}});
    rules.push_back({-2, -1, 7, {}, {}, {}});
    rules.push_back({-2, -1, 8, {}, {}, {}});
    rules.push_back({-2, -1, 9, {}, {}, {}});
    rules.push_back({-2, -3, 1, {}, {}, {}});
    rules.push_back({-2, -4, 2, {}, {}, {}});
    
    rules.push_back({-3, -1, 1, {}, {}, {}});
    rules.push_back({-3, -1, 2, {}, {}, {}});
    rules.push_back({-3, -1, 3, {}, {}, {}});
    rules.push_back({-3, -1, 4, {}, {}, {}});
    rules.push_back({-3, -1, 5, {}, {}, {}});
    rules.push_back({-3, -1, 6, {}, {}, {}});
    rules.push_back({-3, -1, 7, {}, {}, {}});
    rules.push_back({-3, -1, 8, {}, {}, {}});
    rules.push_back({-3, -1, 9, {}, {}, {}});









    rules.push_back({-3, -2, 1, {}, {}, {}});
    
    rules.push_back({-4, -1, 1, {}, {}, {}});
    rules.push_back({-4, -2, 2, {}, {}, {}});
    
    rules.push_back({-5,  -1, 1, {}, {}, {}});
    rules.push_back({-6,  -1, 1, {}, {}, {}});
    rules.push_back({-7,  -1, 1, {}, {}, {}});
    rules.push_back({-8,  -1, 1, {}, {}, {}});
    rules.push_back({-12, -1, 1, {}, {}, {}});
    
    return rules;
}

inline bool rule_allows_base(const AttachmentRule& rule, const CatalogEntry& entry) {
    if (entry.is_nn()) {
        auto info = nn_detail::parse_nn_topo(entry.topo);
        if (!info.valid) return false;
        if (info.is_sidelink && !rule.allowed_sidelinks.empty()) {
            bool ok = false;
            for (int p : rule.allowed_sidelinks) if (p == info.sidelink_param) { ok = true; break; }
            if (!ok) return false;
        }
        if (!info.is_sidelink && !rule.allowed_nobles.empty()) {
            bool ok = false;
            for (const auto& nm : rule.allowed_nobles) if (nm == info.noble_name) { ok = true; break; }
            if (!ok) return false;
        }
    }
    return true;
}

// ############################################################################
// PART 7: DUMMY LST + IF OPERATIONS
// ############################################################################

struct DummyLST {
    int length;
    Eigen::MatrixXi IF;
    int T;
};

inline DummyLST make_dummy_lst(int length) {
    DummyLST d;
    d.length = length;
    d.T = length - 1;
    d.IF = Eigen::MatrixXi::Zero(length, length);
    d.IF(0, 0) = -1;
    d.IF(length - 1, length - 1) = -1;
    for (int i = 1; i < length - 1; i++) d.IF(i, i) = -2;
    for (int i = 0; i < length - 1; i++) {
        d.IF(i, i + 1) = 1; d.IF(i + 1, i) = 1;
    }
    return d;
}

inline std::vector<int> find_curves(const Eigen::MatrixXi& IF, int si) {
    std::vector<int> r;
    for (int i = 0; i < IF.rows(); i++) if (IF(i, i) == si) r.push_back(i);
    return r;
}

inline Eigen::MatrixXi attach_curve(const Eigen::MatrixXi& IF,
                                     int target, int ext_si, int int_num) {
    int n = IF.rows();
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(n + 1, n + 1);
    nIF.block(0, 0, n, n) = IF;
    nIF(n, n) = ext_si;
    nIF(n, target) = int_num;
    nIF(target, n) = int_num;
    return nIF;
}

inline Eigen::MatrixXi glue_dummy(const Eigen::MatrixXi& base,
                                   int target, const DummyLST& dummy, int ep = 0) {
    int nb = base.rows(), nd = dummy.length;
    Eigen::MatrixXi nIF = Eigen::MatrixXi::Zero(nb + nd - 1, nb + nd - 1);
    nIF.block(0, 0, nb, nb) = base;
    auto map_d = [&](int j) -> int {
        if (j == ep) return target;
        return nb + (j < ep ? j : j - 1);
    };
    for (int j = 0; j < nd; j++) {
        if (j == ep) continue;
        nIF(map_d(j), map_d(j)) = dummy.IF(j, j);
    }
    for (int j = 0; j < nd; j++)
        for (int k = j + 1; k < nd; k++)
            if (dummy.IF(j, k) != 0) {
                nIF(map_d(j), map_d(k)) += dummy.IF(j, k);
                nIF(map_d(k), map_d(j)) += dummy.IF(k, j);
            }
    return nIF;
}

// ############################################################################
// PART 8: SUGRA RESULT + CONFIG
// ############################################################################

struct SUGRAResult {
    int catalog_id;
    std::string catalog_type;
    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;
    
    struct Attachment { int ext_si, target_idx, int_num; };
    std::vector<Attachment> externals;
    
    struct DummyGlue { int length, target_idx, endpoint; };
    std::vector<DummyGlue> dummies;
    
    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;
    
    bool valid;
    std::string fail_reason;
};

struct SUGRAConfig {
    int T_max = 193;
    int max_externals = 10;
    int max_dummy = 0;
    int dummy_min_length = 2;
    int dummy_max_length = 10;
    
    double cc_budget = 8.0;
    
    bool check_nhc = true;
    bool check_anomaly = false;
    bool check_determinant = false;
    int required_det = 0;
    bool check_det_perfect_square = false;  // |det| = n² filter
    
    bool include_nn = true;
    bool include_nk = true;
    bool include_dummy = true;
    int catalog_T_min = 0;
    int catalog_T_max = 9999;
    
    bool verbose = false;

    // Large intersection mode
    bool largeint_single = false;  // skip NHC-forming mixed multi-target
    bool largeint_mixed = false;   // skip single-target, only NHC-forming
    int mixed_int_max = 5;  // default: enumerate int_nums up to 5

    // Save non-SUGRA bases (anomaly fail or reject_nonuni fail)
    bool save_nonsugra = false;
};

// ############################################################################
// PART 9: CANDIDATE FINDING
// ############################################################################

inline double used_cc_on_minus1(const Eigen::MatrixXi& IF, int m1) {
    double used = 0;
    for (int j = 0; j < IF.rows(); j++) {
        if (j == m1 || IF(m1, j) == 0 || IF(j, j) == -1) continue;
        used += central_charge(gauge_from_si(IF(j, j)), std::abs(IF(m1, j)));
    }
    return used;
}

// NHC-aware gauge detection for a curve
inline GaugeInfo gauge_on_curve_nhc(const Eigen::MatrixXi& IF, int idx) {
    int si = IF(idx, idx);
    int n = IF.rows();
    if (si == -2) {
        // Part of 2-3 NHC?
        for (int j = 0; j < n; j++)
            if (j != idx && IF(idx, j) != 0 && IF(j, j) == -3) return GAUGE_SU2;
        // Part of 2-4 NHC? (int=2 to -4)
        for (int j = 0; j < n; j++)
            if (j != idx && std::abs(IF(idx, j)) == 2 && IF(j, j) == -4) return GAUGE_SU8;
        return GAUGE_NONE;  // standalone -2
    }
    if (si == -3) {
        // Part of 2-3 NHC?
        for (int j = 0; j < n; j++)
            if (j != idx && IF(idx, j) != 0 && IF(j, j) == -2) return GAUGE_G2;
        return GAUGE_SU3;  // standalone -3
    }
    return gauge_from_si(si);
}

// NHC-aware cc used on a (-1) curve
inline double used_cc_on_minus1_nhc(const Eigen::MatrixXi& IF, int m1) {
    double used = 0;
    int n = IF.rows();
    for (int j = 0; j < n; j++) {
        if (j == m1 || IF(m1, j) == 0 || IF(j, j) == -1) continue;
        GaugeInfo g = gauge_on_curve_nhc(IF, j);
        if (g.dim > 0)
            used += central_charge(g, std::abs(IF(m1, j)));
    }
    return used;
}

// Quick cc check on ALL (-1) curves in IF. Returns false if any exceeds budget.
inline bool quick_cc_check(const Eigen::MatrixXi& IF, double budget = 8.0) {
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (IF(i, i) != -1) continue;
        if (used_cc_on_minus1(IF, i) > budget + 1e-9) return false;
    }
    return true;
}

struct AttachCandidate {
    int rule_idx, target_idx, ext_si, int_num;
};

inline std::vector<AttachCandidate> find_candidates(
    const Eigen::MatrixXi& IF,
    const std::vector<AttachmentRule>& rules,
    double budget = 8.0)
{
    std::vector<AttachCandidate> cands;
    for (int ri = 0; ri < (int)rules.size(); ri++) {
        const auto& rule = rules[ri];
        for (int tidx : find_curves(IF, rule.target_self_int)) {
            if (rule.target_self_int == -1) {
                auto g = gauge_from_si(rule.ext_self_int);
                if (g.dim > 0) {
                    double rem = budget - used_cc_on_minus1(IF, tidx);
                    if (central_charge(g, rule.intersection_num) > rem + 1e-9) continue;
                }
            }
            cands.push_back({ri, tidx, rule.ext_self_int, rule.intersection_num});
        }
    }
    return cands;
}

// ############################################################################
// PART 10: GENERATOR
// ############################################################################

// Dummy LST IF builders
// A(n): 1-2-2-...-2-1, n curves total (n >= 2)
inline Eigen::MatrixXi build_dummy_A(int n) {
    Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
    IF(0, 0) = -1;
    IF(n-1, n-1) = -1;
    for (int i = 1; i < n-1; i++) IF(i, i) = -2;
    for (int i = 0; i < n-1; i++) { IF(i, i+1) = 1; IF(i+1, i) = 1; }
    return IF;
}

// D(n): D-type dummy LST
// n=3: linear chain [-2, -1, -2] (degenerate, no branch)
// n>=4: 2{2,2}2...21
//   Main chain: [0]=-2, [1]=-2, ..., [n-4]=-2, [n-3]=-1  (n-2 curves)
//   Two branches: [n-2]=-2, [n-1]=-2, both connected to [0]
//   Node 0 has degree: 2 (branch) + 1 (chain) = 3  (for n=4: chain is just -1)

inline Eigen::MatrixXi build_dummy_D(int n) {
    if (n == 3) {
        // Linear chain: -2, -1, -2
        Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(3, 3);
        IF(0, 0) = -2; IF(1, 1) = -1; IF(2, 2) = -2;
        IF(0, 1) = 1; IF(1, 0) = 1;
        IF(1, 2) = 1; IF(2, 1) = 1;
        return IF;
    }
    // n >= 4: total n curves
    // Chain: indices 0..n-3, lengths n-2
    //   [0]=-2 (branch point), [1]=-2, ..., [n-4]=-2, [n-3]=-1
    // Two branches: [n-2]=-2, [n-1]=-2, connected to [0]
    Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
    for (int i = 0; i <= n - 4; i++) IF(i, i) = -2;
    IF(n - 3, n - 3) = -1;
    IF(n - 2, n - 2) = -2;
    IF(n - 1, n - 1) = -2;
    // Chain adjacency
    for (int i = 0; i < n - 3; i++) { IF(i, i + 1) = 1; IF(i + 1, i) = 1; }
    // Two branches at node 0
    IF(0, n - 2) = 1; IF(n - 2, 0) = 1;
    IF(0, n - 1) = 1; IF(n - 1, 0) = 1;
    return IF;
}

// Generate from raw IF (for dummy LSTs)
inline std::vector<SUGRAResult> generate_from_IF(
    const Eigen::MatrixXi& base_IF,
    const std::string& label,
    int label_id,
    const std::vector<AttachmentRule>& rules,
    const SUGRAConfig& config)
{
    std::vector<SUGRAResult> results;
    
    auto base_sig = compute_sig_fast(base_IF);
    if (base_sig.sig_neg > config.T_max) return results;
    
    auto cands = find_candidates(base_IF, rules, config.cc_budget);
    
    for (const auto& c : cands) {
        // No rule_allows_base check for dummy LSTs (all rules allowed)
        
        Eigen::MatrixXi new_IF = attach_curve(base_IF, c.target_idx, c.ext_si, c.int_num);
        auto new_sig = compute_sig_fast(new_IF);
        
        if (new_sig.sig_neg > config.T_max) continue;
        if (config.check_determinant && std::abs(new_sig.det) != 1) continue;
        if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) continue;
        
        NHCResult nhc = {true, "", {}, {}, {}};
        if (config.check_nhc) {
            nhc = check_nhc(new_IF, config.cc_budget);
            if (!nhc.passes) continue;
        }
        
        // External -2 → -1: su(2) gauge enhancement
        enhance_external_m2_gauge(nhc, new_IF, new_IF.rows() - 1,
                                  base_IF(c.target_idx, c.target_idx), config.cc_budget);
        if (!nhc.passes) continue;
        
        AnomalyResult anom = {0, 0, 0, 0};
        if (nhc.passes) {
            anom = compute_anomaly(new_IF, nhc);
            if (config.check_anomaly && anom.H_neutral < 0) continue;
        }
        
        SUGRAResult r;
        r.catalog_id = label_id;
        r.catalog_type = label;
        r.base_IF = base_IF;
        r.final_IF = new_IF;
        r.externals.push_back({c.ext_si, c.target_idx, c.int_num});
        r.nhc = nhc;
        r.anomaly = anom;
        r.sig = new_sig;
        r.valid = true;
        results.push_back(r);
    }
    
    return results;
}

inline std::vector<SUGRAResult> generate_from_entry(
    const CatalogEntry& entry,
    const std::vector<AttachmentRule>& rules,
    const SUGRAConfig& config)
{
    std::vector<SUGRAResult> results;
    
    Eigen::MatrixXi base_IF = reconstruct_IF(entry);
    if (base_IF.rows() == 0) return results;
    
    auto base_sig = compute_sig_fast(base_IF);
    if (base_sig.sig_neg > config.T_max) return results;
    
    auto cands = find_candidates(base_IF, rules, config.cc_budget);
    
    for (const auto& c : cands) {
        if (!rule_allows_base(rules[c.rule_idx], entry)) continue;
        
        Eigen::MatrixXi new_IF = attach_curve(base_IF, c.target_idx, c.ext_si, c.int_num);
        auto new_sig = compute_sig_fast(new_IF);
        
        if (new_sig.sig_neg > config.T_max) continue;
        if (config.check_determinant && std::abs(new_sig.det) != 1) continue;
        if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) continue;
        
        NHCResult nhc = {true, "", {}, {}, {}};
        if (config.check_nhc) {
            nhc = check_nhc(new_IF, config.cc_budget);
            if (!nhc.passes) continue;
        }
        
        // External -2 → -1: su(2) gauge enhancement
        enhance_external_m2_gauge(nhc, new_IF, new_IF.rows() - 1,
                                  base_IF(c.target_idx, c.target_idx), config.cc_budget);
        if (!nhc.passes) continue;
        
        AnomalyResult anom = {0, 0, 0, 0};
        if (nhc.passes) {
            anom = compute_anomaly(new_IF, nhc);
            if (config.check_anomaly && anom.H_neutral < 0) continue;
        }
        
        SUGRAResult r;
        r.catalog_id = entry.id;
        r.catalog_type = entry.type;
        r.base_IF = base_IF;
        r.final_IF = new_IF;
        r.externals.push_back({c.ext_si, c.target_idx, c.int_num});
        r.nhc = nhc;
        r.anomaly = anom;
        r.sig = new_sig;
        r.valid = true;
        results.push_back(r);
    }
    
    return results;
}

inline std::vector<SUGRAResult> generate_sugra(
    const std::vector<CatalogEntry>& catalog,
    const SUGRAConfig& config)
{
    auto rules = build_attachment_rules();
    std::vector<SUGRAResult> all;
    
    // ── Catalog entries ──
    int processed = 0;
    for (const auto& entry : catalog) {
        if (entry.is_nn() && !config.include_nn) continue;
        if (!entry.is_nn() && !config.include_nk) continue;
        if (entry.T < config.catalog_T_min || entry.T > config.catalog_T_max) continue;
        
        auto res = generate_from_entry(entry, rules, config);
        for (auto& r : res) all.push_back(std::move(r));
        
        processed++;
        if (config.verbose && processed % 1000 == 0)
            std::cout << processed << " → " << all.size() << " bases\n";
    }
    if (config.verbose)
        std::cout << "Catalog: " << processed << " → " << all.size() << " bases\n";
    
    // ── Dummy LSTs ──
    if (config.include_dummy) {
        int dummy_count = 0;
        
        // A-type: 1-2-...-2-1, n curves → T_base = n-1
        // Cap: catalog_T_min ≤ T_base ≤ catalog_T_max
        int a_min = std::max(2, config.catalog_T_min + 1);  // n = T_base + 1
        int a_max = config.catalog_T_max + 1;
        for (int n = a_min; n <= a_max; n++) {
            auto IF = build_dummy_A(n);
            auto res = generate_from_IF(IF, "DM:A(" + std::to_string(n) + ")", -n, rules, config);
            for (auto& r : res) all.push_back(std::move(r));
            dummy_count += res.size();
        }
        
        // D-type: D(3)=212, D(4)=2{2,2}1, D(5)=2{2,2}21, ...
        int d_min = std::max(3, config.catalog_T_min + 1);
        int d_max = config.catalog_T_max + 1;
        for (int n = d_min; n <= d_max; n++) {
            auto IF = build_dummy_D(n);
            auto res = generate_from_IF(IF, "DM:D(" + std::to_string(n) + ")", -(1000+n), rules, config);
            for (auto& r : res) all.push_back(std::move(r));
            dummy_count += res.size();
        }
        
        if (config.verbose)
            std::cout << "Dummy LSTs: " << dummy_count << " bases\n";
    }
    
    if (config.verbose)
        std::cout << "Total: " << all.size() << " bases\n";
    return all;
}

// ############################################################################
// PART 11: OUTPUT
// ############################################################################

inline void print_sugra_summary(const std::vector<SUGRAResult>& results) {
    std::map<int, int> T_dist, ext_dist;
    int anom_pass = 0;
    for (const auto& r : results) {
        T_dist[r.anomaly.T]++;
        for (const auto& e : r.externals) ext_dist[e.ext_si]++;
        if (r.anomaly.H_neutral >= 0) anom_pass++;
    }
    std::cout << "SUGRA: " << results.size() << " bases | H_neutral>=0: " << anom_pass << "\n";
    std::cout << "T distribution:\n";
    for (auto& [T, cnt] : T_dist) std::cout << "  T=" << T << ": " << cnt << "\n";
    std::cout << "Externals:\n";
    for (auto& [si, cnt] : ext_dist) std::cout << "  " << si << ": " << cnt << "\n";
}

inline void print_gauge_assignment(const SUGRAResult& r) {
    std::cout << "ID=" << r.catalog_id << " " << r.catalog_type
              << " T=" << r.anomaly.T << " Hc=" << r.anomaly.H_charged
              << " V=" << r.anomaly.V << " Hn=" << r.anomaly.H_neutral << "\n";
    std::cout << "  IF " << r.final_IF.rows() << "x" << r.final_IF.cols()
              << " det=" << r.sig.det << "\n";
    for (int ci = 0; ci < (int)r.nhc.clusters.size(); ci++) {
        auto& cl = r.nhc.clusters[ci];
        std::cout << "  cluster " << ci << ": ";
        if (cl.nhc) std::cout << cl.nhc->name;
        else std::cout << "pure-2";
        std::cout << " H=" << cl.H << " V=" << cl.V << " curves={";
        for (int i = 0; i < (int)cl.curve_indices.size(); i++) {
            if (i) std::cout << ",";
            int c = cl.curve_indices[i];
            std::cout << c << "(" << r.final_IF(c, c) << ")";
        }
        std::cout << "}\n";
    }
}

// ############################################################################
// PART 8: MULTI-TARGET UTILITIES (v6 + hat1 multi)
// ############################################################################

// Max level k such that cc(g, k) ≤ budget
inline int max_level_for_cc(const GaugeInfo& g, double budget = 20.0) {
    if (g.dim == 0) return 999;
    if ((double)g.dim <= budget) return 999;  // cc(g,k) < dim always
    // cc(g,k) = k*dim/(k+h) ≤ budget  →  k ≤ budget*h/(dim-budget)
    return std::max(1, (int)(budget * g.h_dual / (g.dim - budget)));
}

// Blowdown for level extraction:
// Build extended IF (with b₀), blowdown all exceptional curves
// (C²=-1 AND b₀·C=1, i.e. non-hat1 (-1) curves),
// then extract level = intersection between 0-curve and positive-curve.
inline int extract_level_from_blowdown(const Eigen::MatrixXi& IF,
                                        const std::set<int>& hat1_set) {
    int n = IF.rows();
    // Build extended IF with b₀
    int T = 0;
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
        for (int i = 0; i < n; i++)
            if (es.eigenvalues()(i) < -1e-8) T++;
    }
    Eigen::MatrixXi cur = Eigen::MatrixXi::Zero(n+1, n+1);
    cur.block(0, 0, n, n) = IF;
    for (int i = 0; i < n; i++) {
        int b = hat1_set.count(i) ? -1 : (IF(i,i) + 2);
        cur(i, n) = b; cur(n, i) = b;
    }
    cur(n, n) = 9 - T;

    // Track hat1 indices through blowdowns
    std::set<int> h1 = hat1_set;

    // Blowdown loop: remove non-hat1 (-1) curves with b₀·C = 1
    while (true) {
        int nn = cur.rows(), m = -1;
        for (int i = 0; i < nn - 1; i++) {  // skip b₀ (last)
            if (h1.count(i)) continue;        // skip hat1
            if (cur(i, i) == -1 && cur(i, nn-1) == 1) { m = i; break; }
        }
        if (m < 0) break;
        // Blowdown curve m
        Eigen::MatrixXi next = Eigen::MatrixXi::Zero(nn-1, nn-1);
        for (int i = 0, ni = 0; i < nn; i++) {
            if (i == m) continue;
            for (int j = 0, nj = 0; j < nn; j++) {
                if (j == m) continue;
                next(ni, nj) = cur(i, j) + cur(i, m) * cur(m, j);
                nj++;
            }
            ni++;
        }
        cur = next;
        std::set<int> nh;
        for (int idx : h1) nh.insert(idx > m ? idx - 1 : idx);
        h1 = nh;
    }

    // Extract level: find 0-curve intersecting positive-self-int geometric curve
    // b₀ is the last index (nn-1) — exclude it from search
    int nn = cur.rows();
    for (int i = 0; i < nn - 1; i++) {          // geometric curves only
        if (cur(i, i) != 0) continue;
        // Found 0-curve, find intersection with positive-self-int geometric curve
        for (int j = 0; j < nn - 1; j++) {      // exclude b₀
            if (j == i) continue;
            if (cur(j, j) > 0)
                return std::abs(cur(i, j));
        }
        // Only one 0-curve, no positive geometric curve → standard Hirzebruch, level 1
        return 1;
    }
    return 1;  // no 0-curve found → level 1 (default)
}

// Check if multi-target attachment passes level-based cc bound
// Returns true if passes, false if rejected
inline bool check_multi_target_cc(const Eigen::MatrixXi& new_IF,
                                   const GaugeInfo& ext_gauge,
                                   const std::set<int>& hat1_set,
                                   double cc_bound = 20.0) {
    if (ext_gauge.dim == 0) return true;
    int level = extract_level_from_blowdown(new_IF, hat1_set);
    double cc = (double)(level * ext_gauge.dim) / (double)(level + ext_gauge.h_dual);
    return cc <= cc_bound + 1e-9;
}

// Enumerate all int_num assignments for n_targets, each from 1 to max_k
template<typename Func>
inline void enumerate_int_nums(int n_targets, int max_k, Func&& func) {
    std::vector<int> nums(n_targets, 1);
    while (true) {
        func(nums);
        int i = n_targets - 1;
        while (i >= 0) {
            nums[i]++;
            if (nums[i] <= max_k) break;
            nums[i] = 1;
            i--;
        }
        if (i < 0) break;
    }
}

// Enumerate k-element subsets of pool, calling func(subset)
template<typename Func>
inline void enumerate_combinations_v2(const std::vector<int>& pool, int k, Func&& func) {
    int n = pool.size();
    if (k > n || k < 1) return;
    std::vector<int> idx(k);
    for (int i = 0; i < k; i++) idx[i] = i;
    while (true) {
        std::vector<int> targets(k);
        for (int i = 0; i < k; i++) targets[i] = pool[idx[i]];
        func(targets);
        int i = k - 1;
        while (i >= 0 && idx[i] == n - k + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < k; j++) idx[j] = idx[j-1] + 1;
    }
}

template<typename Func>
inline void enumerate_subsets_v2(const std::vector<int>& pool, int min_k, int max_k, Func&& func) {
    for (int k = min_k; k <= std::min(max_k, (int)pool.size()); k++)
        enumerate_combinations_v2(pool, k, func);
}

// ############################################################################
// Non-unimodular T_min/T_max rejection
// ############################################################################

// f_multi[T] = min H-V+29T achievable by any LST partition of T tensors
static const int f_multi[] = {
    0, 29, 58, 79, 88, 93, 96, 98, 99, 110, 110,  // T=0..10
    112, 100, 111, 112, 113, 114, 115, 116, 117, 116,  // T=11..20
    117, 117, 117, 118, 118, 118, 118, 119, 119, 119,  // T=21..30
    119, 119, 119, 120, 121, 120, 119, 120, 121, 120,  // T=31..40
    121, 120, 122, 120, 125, 123, 123, 120, 121, 126,  // T=41..50
    121, 122, 121, 132, 153, 158, 124, 125, 121, 150,  // T=51..60
    128, 123, 132, 122, 151, 180, 201, 125, 126, 123,  // T=61..70
    152, 129, 124, 133, 125, 154, 183, 204, 126, 127,  // T=71..80
    128, 157, 130, 132, 134, 135, 164, 193, 214, 136,  // T=81..90
    148, 157, 162, 183, 211, 232, 234, 235, 241, 242,  // T=91..100
    241, 236, 241, 242, 243, 242, 241, 242, 244, 242,  // T=101..110
    243, 242, 243, 243, 243, 244, 243, 242, 244, 247,  // T=111..120
    244, 245, 243, 245, 246, 245, 246, 244, 244, 247,  // T=121..130
    248, 245, 248, 245, 247, 250, 246, 247, 247, 246,  // T=131..140
    249, 250, 247, 249, 248, 248, 251, 249, 249, 250,  // T=141..150
    251, 250, 251, 251, 252, 253, 256, 252, 253, 254,  // T=151..160
    255, 256, 257, 258, 260, 260, 262, 264, 262, 263,  // T=161..170
    264, 270, 266, 268, 270, 271, 283, 292, 297, 272,  // T=171..180
    284, 293, 298, 310, 319, 324, 345, 365, 367, 367,  // T=181..190
    366, 366, 99999  // T=191..193
};

// max T_add such that f_multi[T_add] <= budget
inline int max_t_add_from_fmulti(int budget) {
    if (budget < 0) return 0;
    int tmax = 0;
    for (int t = 1; t <= 192; t++) {
        if (f_multi[t] <= budget) tmax = t;
    }
    return tmax;
}

// Max possible Δ(V-H) recovery from remaining cc budget (--use-lst-T mode, no 29T penalty)
// Precomputed greedy table: best externals by Δ(V-H)/cc ratio
// e8(248/8.0=31), e7(133/7.0=19), e6(78/6.0=13), f4(52/5.2=10), so8(28/4.0=7), su3(8/2.0=4)
inline int max_recovery_from_cc(double remaining_cc) {
    if (remaining_cc < 2.0 - 1e-9) return 0;
    // Greedy: fill with highest ratio externals
    struct ExtOpt { double cc; int dvh; };
    static const ExtOpt opts[] = {
        {8.0, 248},  // e8
        {7.0, 133},  // e7
        {6.0, 78},   // e6
        {5.2, 52},   // f4
        {4.0, 28},   // so8
        {2.0, 8},    // su3
    };
    double rem = remaining_cc;
    int total = 0;
    for (auto& o : opts) {
        while (rem >= o.cc - 1e-9) {
            rem -= o.cc;
            total += o.dvh;
        }
    }
    return total;
}

// Check if an anomaly-failing entry can possibly recover with more externals.
// Returns true if recovery is impossible (should skip, not save to non-SUGRA).
inline bool recovery_impossible(int H_neutral, double remaining_total_cc) {
    return H_neutral + max_recovery_from_cc(remaining_total_cc) < 0;
}

inline bool reject_nonuni(const Eigen::MatrixXi& IF, const SigInfo& sig,
                           const AnomalyResult& anom, const std::set<int>& hat1_set) {
    // H_neutral at T_sugra = current tensor count
    int T_sugra = anom.T;
    int H_neutral = 273 + anom.V - anom.H_charged - 29 * T_sugra;
    // T_max = T_sugra + max additional tensors from LST content
    int T_max_anom = (H_neutral >= 0) ?
        std::min(T_sugra + max_t_add_from_fmulti(H_neutral), 193) : -1;

    if (std::abs(sig.det) <= 1) {
        // Unimodular: check extended sig at T = sig_neg
        // (T <= T_max is guaranteed by anomaly check upstream)
        if (check_extended_sig(IF, sig.sig_neg, hat1_set)) return false;  // sig OK → keep
        // sig_pos != 1 at T=sig_neg: fall through to T_min calculation like non-uni
    }

    // Compute T_min from det=0 crossing of extended IF
    int n = IF.rows();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n+1, n+1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        double b = hat1_set.count(i) ? -1.0 : (IF(i,i) + 2.0);
        M(i, n) = b; M(n, i) = b;
    }
    M(n, n) = 0; double B = M.determinant();
    M(n, n) = 1; double A = M.determinant() - B;
    if (std::abs(A) < 1e-10) {
        // A=0: det is constant in T (no crossing). Fallback: scan T for sig_pos=1.
        for (int T = sig.sig_neg; T <= T_max_anom; T++) {
            if (check_extended_sig(IF, T, hat1_set)) return false;  // found valid T → keep
        }
        return true;  // no valid T found → reject
    }
    double T_crit = 9.0 - (-B / A);
    int T_min = (int)std::ceil(T_crit - 1e-9);
    if (T_min < 0 || T_min > T_max_anom) return true;

    // Check extended sig at T_min
    return !check_extended_sig(IF, T_min, hat1_set);
}

// end of sugra_generator.h
