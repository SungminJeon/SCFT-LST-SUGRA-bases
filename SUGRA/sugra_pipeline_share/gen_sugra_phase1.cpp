// gen_sugra_phase1.cpp
// Phase 1: Single-external SUGRA base catalog
//
// Follows v2 algorithm exactly for standard externals,
// v5 algorithm for hat{1} externals.
//
// Output: per-spec .cat files with stacking metadata
//   → Phase 2 reads these and stacks additional externals
//
// For each ExternalSpec:
//   for each catalog entry:
//     find valid target curves → attach single external → filter → record
//   + dummy LSTs (A-type, D-type)
//   deduplicate (v2-style: eigenvalue + VDesc + edge)
//   write spec .cat file
//
// Compile:
//   g++ -std=c++17 -O2 -o gen_sugra_phase1 gen_sugra_phase1.cpp Tensor.C \
//       TopoLineCompact_enhanced.cpp Topology_enhanced.cpp \
//       -I/usr/include/eigen3 -I.
//
// Usage:
//   ./gen_sugra_phase1 unified.cat [T_max=193] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>
#include <map>
#include <algorithm>
#include <functional>
#include <filesystem>
#include <mutex>
#include <atomic>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

static double g_max_ext_cc = 16.0;

// ============================================================================
// PART 1: ExternalSpec — reused from enhanced, trimmed
// ============================================================================

struct ExternalSpec {
    int id;
    int ext_si;
    int target_si;
    int int_num;
    bool is_hat1;
    double cc;
    std::string label;
    std::string tag;
};

inline double compute_ext_cc(int ext_si, int int_num, bool is_hat1, int target_si) {
    if (is_hat1) {
        if (target_si == -1) return 7.0;
        if (target_si == -2) return 15.0;
        return 10.0;
    }
    GaugeInfo g = gauge_from_si(ext_si);
    if (g.dim > 0)
        return central_charge(g, int_num);
    if (ext_si == -2 && target_si == -1)
        return central_charge(GAUGE_SU2, int_num);
    if (ext_si == -2 && target_si == -3) return central_charge(GAUGE_SU2, int_num);   // su2
    if (ext_si == -2 && target_si == -4) return central_charge(GAUGE_SU8, int_num);   // su8
    if (ext_si == -3 && target_si == -2) return central_charge(GAUGE_G2, int_num);    // g2
    if (ext_si == -4 && target_si == -2) return central_charge(GAUGE_SO16, int_num);  // so16
    return 0.0;
}

// Total ext cc for multi-target: sum cc over all int_nums
inline double total_ext_cc_multi(int ext_si, const std::vector<int>& int_nums, bool is_hat1, int target_si) {
    double total = 0.0;
    for (int k : int_nums)
        total += compute_ext_cc(ext_si, k, is_hat1, target_si);
    return total;
}

inline std::vector<ExternalSpec> build_all_specs() {
    std::vector<ExternalSpec> specs;
    int id = 0;
    auto add = [&](int ext_si, int tgt_si, int inum, bool hat1,
                    const std::string& lbl, const std::string& tag) {
        double cc = compute_ext_cc(ext_si, inum, hat1, tgt_si);
        specs.push_back({id++, ext_si, tgt_si, inum, hat1, cc, lbl, tag});
    };

    // Standard externals
    for (int k = 1; k <= 2; k++)
        add(-2, -1, k, false, "su(2) k=" + std::to_string(k), "su2");
    add(-2,  -3, 1, false, "(-2)->(-3)",  "su2n3");
    add(-2,  -4, 2, false, "su8->so16",   "su8");
    for (int k = 1; k <= 2; k++)
        add(-3, -1, k, false, "su(3) k=" + std::to_string(k), "su3");
    add(-3,  -2, 1, false, "(-3)->(-2)",  "su3n2");
    add(-4,  -1, 1, false, "so(8)",       "so8");
    add(-4,  -2, 2, false, "(-4)->(-2)",  "so16n2");
    add(-5,  -1, 1, false, "f4",          "f4");
    add(-6,  -1, 1, false, "e6",          "e6");
    add(-7,  -1, 1, false, "e7'",         "e7p");
    add(-8,  -1, 1, false, "e7",          "e7");
    add(-12, -1, 1, false, "e8",          "e8");

    // hat{1}
    add(-1, -1, 1, true, "hat1->(-1)",  "hat1m1");
    add(-1, -2, 1, true, "hat1->(-2)",  "hat1m2");

    return specs;
}

// ============================================================================
// PART 2: Phase1Result — single-external result with stacking metadata
// ============================================================================

struct Phase1Result {
    // Identity
    int spec_id;
    int catalog_id;
    std::string catalog_type;
    int base_T;

    // Attachment
    int target_idx;            // curve in base_IF where external attaches
    int ext_curve_idx;         // index of external in final_IF (= base_IF.rows())
    std::vector<int> hat1_multi_targets;  // for hat1 multi-target (v5 style)

    // Matrices
    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;

    // Physics
    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;

    // Cached eigenvalues for dedup (avoids redundant eigensolve)
    std::vector<double> eigenvalues;

    // Stacking metadata: available targets in final_IF
    // For each curve in final_IF, record (curve_idx, self_int, remaining_cc_on_adjacent_m1)
    struct AvailTarget {
        int curve_idx;
        int self_int;
        double available_cc;   // remaining cc budget on this -1 curve (if self_int == -1)
    };
    std::vector<AvailTarget> remaining_targets;
};

// ============================================================================
// PART 3: Dedup key — v2 style (eigenvalue + VDesc + edge)
// ============================================================================

// Use pre-computed eigenvalues to avoid redundant eigensolve
inline std::string dedup_key_v2_cached(const Eigen::MatrixXi& IF,
                                        const std::vector<double>& eigenvalues,
                                        const std::string& extra = "") {
    int n = IF.rows();

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8);

    // 1. Sorted eigenvalues (already sorted from SelfAdjointEigenSolver)
    ss << "E:";
    for (int i = 0; i < n; i++) {
        double v = eigenvalues[i];
        if (std::abs(v) < 1e-10) v = 0.0;
        if (i) ss << ",";
        ss << v;
    }

    // 2. Sorted vertex descriptors (v2 style)
    struct VDesc {
        int si, deg;
        std::vector<int> nbr_si;
        bool operator<(const VDesc& o) const {
            if (si != o.si) return si < o.si;
            if (deg != o.deg) return deg < o.deg;
            return nbr_si < o.nbr_si;
        }
    };
    std::vector<VDesc> vds(n);
    for (int i = 0; i < n; i++) {
        vds[i].si = IF(i, i);
        vds[i].deg = 0;
        for (int j = 0; j < n; j++) {
            if (j != i && IF(i, j) != 0) {
                vds[i].deg++;
                vds[i].nbr_si.push_back(IF(j, j));
            }
        }
        std::sort(vds[i].nbr_si.begin(), vds[i].nbr_si.end());
    }
    std::sort(vds.begin(), vds.end());
    ss << "|V:";
    for (auto& vd : vds) {
        ss << "(" << vd.si << "," << vd.deg << ",{";
        for (int k = 0; k < (int)vd.nbr_si.size(); k++) {
            if (k) ss << ",";
            ss << vd.nbr_si[k];
        }
        ss << "})";
    }

    // 3. Sorted edge multiset
    std::vector<std::tuple<int, int, int>> edges;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (IF(i, j) != 0) {
                int a = IF(i, i), b = IF(j, j);
                if (a > b) std::swap(a, b);
                edges.push_back({a, b, IF(i, j)});
            }
    std::sort(edges.begin(), edges.end());
    ss << "|E:";
    for (auto& [a, b, w] : edges) ss << "(" << a << "," << b << "," << w << ")";

    // 4. Extra discriminator (spec tag for Phase1, combo_key for Phase2)
    if (!extra.empty()) ss << "|X:" << extra;

    return ss.str();
}

inline std::string dedup_key_v2(const Eigen::MatrixXi& IF, const std::string& extra = "") {
    int n = IF.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    std::vector<double> ev_vec(n);
    for (int i = 0; i < n; i++) ev_vec[i] = es.eigenvalues()(i);
    return dedup_key_v2_cached(IF, ev_vec, extra);
}

// ============================================================================
// PART 4: base_spec_key — v2 style (eigenvalue spectrum for grouping)
// ============================================================================

inline std::string base_spec_key_v2(const Eigen::MatrixXi& IF) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    auto ev = es.eigenvalues();
    std::vector<double> vals(ev.data(), ev.data() + ev.size());
    std::sort(vals.begin(), vals.end());
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < vals.size(); i++) {
        if (i) ss << ",";
        ss << vals[i];
    }
    return ss.str();
}

// ============================================================================
// PART 5: Compute remaining targets in final_IF (stacking metadata)
// ============================================================================

inline std::vector<Phase1Result::AvailTarget> compute_remaining_targets(
    const Eigen::MatrixXi& final_IF,
    const NHCResult& nhc,
    double cc_budget = 8.0)
{
    std::vector<Phase1Result::AvailTarget> targets;
    int n = final_IF.rows();

    for (int i = 0; i < n; i++) {
        int si = final_IF(i, i);
        Phase1Result::AvailTarget t;
        t.curve_idx = i;
        t.self_int = si;
        t.available_cc = -1.0;

        if (si == -1) {
            // Compute remaining cc budget
            // Include (-1) neighbors that have gauge (e.g. hat1 with su(8))
            double used = 0.0;
            for (int j = 0; j < n; j++) {
                if (j == i || final_IF(i, j) == 0) continue;
                GaugeInfo g = (j < (int)nhc.curve_gauges.size() && nhc.curve_gauges[j].dim > 0)
                    ? nhc.curve_gauges[j]
                    : gauge_from_si(final_IF(j, j));
                if (final_IF(j, j) == -1 && g.dim == 0) continue;
                int k = std::abs(final_IF(i, j));
                if (g.dim > 0)
                    used += central_charge(g, k);
            }
            t.available_cc = cc_budget - used;
        }
        targets.push_back(t);
    }
    return targets;
}

// reject_nonuni is now in sugra_generator.h

// ============================================================================
// PART 5c: Combination enumeration (replaces bitmask for large nm1)
// ============================================================================

// Calls func(targets) for each k-element subset of pool
template<typename Func>
inline void enumerate_combinations(const std::vector<int>& pool, int k, Func&& func) {
    int n = pool.size();
    if (k > n || k < 1) return;
    std::vector<int> idx(k);
    for (int i = 0; i < k; i++) idx[i] = i;
    while (true) {
        std::vector<int> targets(k);
        for (int i = 0; i < k; i++) targets[i] = pool[idx[i]];
        func(targets);
        // Next combination
        int i = k - 1;
        while (i >= 0 && idx[i] == n - k + i) i--;
        if (i < 0) break;
        idx[i]++;
        for (int j = i + 1; j < k; j++) idx[j] = idx[j-1] + 1;
    }
}

// Calls func(targets) for each subset of pool with size in [min_k, max_k]
template<typename Func>
inline void enumerate_subsets(const std::vector<int>& pool, int min_k, int max_k, Func&& func) {
    for (int k = min_k; k <= max_k; k++)
        enumerate_combinations(pool, k, func);
}

// ============================================================================
// PART 6: Generate single-external results for one spec
// ============================================================================

// Forward declare MultiTargetKey for skip parameter
struct MultiTargetKey;

struct CachedIF {
    Eigen::MatrixXi IF;
    SigInfo sig;
    bool valid;
};

struct Phase1SpecResult {
    std::vector<Phase1Result> sugra;
    std::vector<Phase1Result> nonsugra;
};

inline Phase1SpecResult generate_single_spec(
    const ExternalSpec& spec,
    const std::vector<CatalogEntry>& catalog,
    const SUGRAConfig& config,
    std::set<std::string>* early_dedup = nullptr,
    const void* skip_multi_target = nullptr,
    const std::vector<CachedIF>* cached_ifs = nullptr)
{
    std::vector<Phase1Result> results;
    std::vector<Phase1Result> nonsugra_results;
    // Helper: build dedup extra string for this spec
    auto make_extra = [&](const Phase1Result& r) -> std::string {
        std::string extra = spec.tag;
        if (spec.is_hat1)
            extra += "|T:" + std::to_string(spec.target_si)
                   + "|N:" + std::to_string(r.hat1_multi_targets.size());
        // Multi-target: use actual int pattern for dedup (not spec tag)
        if (!r.hat1_multi_targets.empty() && (int)r.hat1_multi_targets.size() >= 2 && !spec.is_hat1) {
            // Build pattern from IF: ext curve's connections
            int ext = r.ext_curve_idx;
            int n = r.final_IF.rows();
            std::vector<int> ints;
            for (int j = 0; j < n; j++) {
                if (j == ext) continue;
                int v = r.final_IF(ext, j);
                if (v != 0) ints.push_back(v);
            }
            std::sort(ints.begin(), ints.end());
            extra = "mt" + std::to_string(spec.ext_si) + "_";
            for (int k = 0; k < (int)ints.size(); k++) {
                if (k) extra += "_";
                extra += std::to_string(ints[k]);
            }
        }
        // Include catalog_id (LST source) so different build paths are not merged
        extra += "|cid:" + std::to_string(r.catalog_id);
        return extra;
    };
    // Helper: push result with optional early dedup (thread-safe)
    std::mutex results_mutex;
    auto push_result = [&](Phase1Result& r, const Eigen::MatrixXi& final_IF, const NHCResult& nhc_ref) {
        r.remaining_targets = compute_remaining_targets(final_IF, nhc_ref, config.cc_budget);
        std::lock_guard<std::mutex> lock(results_mutex);
        if (early_dedup) {
            std::string extra = make_extra(r);
            std::string key = r.eigenvalues.empty()
                ? dedup_key_v2(final_IF, extra)
                : dedup_key_v2_cached(final_IF, r.eigenvalues, extra);
            if (!early_dedup->insert(key).second) return; // duplicate, skip
        }
        results.push_back(std::move(r));
    };
    auto push_nonsugra = [&](Phase1Result& r, const Eigen::MatrixXi& final_IF, const NHCResult& nhc_ref) {
        if (!config.save_nonsugra) return;
        // Skip if recovery is impossible
        double cur_cc = spec.cc;  // Phase 1: single external
        if (recovery_impossible(r.anomaly.H_neutral, g_max_ext_cc - cur_cc)) return;
        r.remaining_targets = compute_remaining_targets(final_IF, nhc_ref, config.cc_budget);
        std::lock_guard<std::mutex> lock(results_mutex);
        nonsugra_results.push_back(std::move(r));
    };

    // --- Catalog entries (OpenMP parallelized) ---

    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t ei = 0; ei < catalog.size(); ei++) {
        const auto& entry = catalog[ei];
        if (entry.is_nn() && !config.include_nn) continue;
        if (!entry.is_nn() && !config.include_nk) continue;
        if (entry.T < config.catalog_T_min || entry.T > config.catalog_T_max) continue;

        Eigen::MatrixXi base_IF;
        SigInfo base_sig;
        if (cached_ifs && ei < cached_ifs->size() && (*cached_ifs)[ei].valid) {
            base_IF = (*cached_ifs)[ei].IF;
            base_sig = (*cached_ifs)[ei].sig;
        } else {
            base_IF = reconstruct_IF(entry);
            if (base_IF.rows() == 0) continue;
            base_sig = compute_sig(base_IF);
        }
        if (base_IF.rows() == 0) continue;
        if (base_sig.sig_neg + 1 > config.T_max) continue;  // +1 for the external

        int n_base = base_IF.rows();
        g_current_base_T = n_base;  // for --use-lst-T mode: LST curve count

        // Pre-compute cc used on each (-1) curve (once per entry, reuse across all specs)
        std::vector<int> all_m1_curves;
        std::vector<double> m1_used_cc;  // parallel to all_m1_curves
        for (int i = 0; i < n_base; i++) {
            if (base_IF(i, i) == -1) {
                all_m1_curves.push_back(i);
                m1_used_cc.push_back(used_cc_on_minus1(base_IF, i));
            }
        }

        // Section timers (thread-safe accumulation)
        auto t_start = [](){ return std::chrono::steady_clock::now(); };
        auto t_elapsed_ms = [](std::chrono::steady_clock::time_point s){
            return std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - s).count();
        };
        static std::atomic<long long> time_A{0}, time_B{0}, time_D{0}, time_E{0}, time_F{0};
        static std::atomic<int> count_A{0}, count_B{0}, count_D{0}, count_E{0}, count_F{0};
        static std::atomic<int> entry_count{0};
        int ec = ++entry_count;
        if (ec % 100 == 0) {
            std::lock_guard<std::mutex> lock(results_mutex);
            std::cerr << "[entry " << ec << "/" << catalog.size() << "] "
                      << "A=" << time_A/1000 << "s(" << count_A << ") "
                      << "B=" << time_B/1000 << "s(" << count_B << ") "
                      << "D=" << time_D/1000 << "s(" << count_D << ") "
                      << "E=" << time_E/1000 << "s(" << count_E << ") "
                      << "F=" << time_F/1000 << "s(" << count_F << ")\n";
        }

        auto tA = t_start();
        // Find target curves matching spec.target_si
        // For (-1) targets: use cached cc to pre-filter
        std::vector<int> target_curves;
        if (!config.largeint_mixed) {
            if (spec.target_si == -1) {
                GaugeInfo ext_g = gauge_from_si(spec.ext_si);
                if (ext_g.dim == 0 && spec.ext_si == -2) ext_g = GAUGE_SU2;
                // hat1 → (-1): su(8) gauge, hat1 → (-2): su(16) gauge
                if (spec.is_hat1 && spec.target_si == -1) ext_g = GAUGE_SU8;
                if (spec.is_hat1 && spec.target_si == -2) ext_g = GAUGE_SU16;
                double ext_cc = (ext_g.dim > 0) ? central_charge(ext_g, spec.int_num) : 0;
                for (int mi = 0; mi < (int)all_m1_curves.size(); mi++) {
                    if (m1_used_cc[mi] + ext_cc <= config.cc_budget + 1e-9)
                        target_curves.push_back(all_m1_curves[mi]);
                }
            } else {
                for (int i = 0; i < n_base; i++) {
                    if (base_IF(i, i) != spec.target_si) continue;
                    // hat1→(-2): target must be standalone (not part of NHC chain)
                    if (spec.is_hat1 && spec.target_si == -2) {
                        bool standalone = true;
                        for (int j = 0; j < n_base; j++) {
                            if (j == i || base_IF(i, j) == 0) continue;
                            if (base_IF(j, j) != -1) { standalone = false; break; }
                        }
                        if (!standalone) continue;
                    }
                    target_curves.push_back(i);
                }
            }
        }

        if (target_curves.empty() && !config.largeint_mixed) continue;

        for (int tidx : target_curves) {
            // Attach
            Eigen::MatrixXi new_IF = attach_curve(base_IF, tidx, spec.ext_si, spec.int_num);
            int ext_idx = new_IF.rows() - 1;

            // 1. LDLT sig_pos prefilter (cheapest)
            if (!spec.is_hat1 && sig_pos_exceeds_one_fast(new_IF, config.T_max)) continue;

            // 2. Full eigensolve
            std::vector<double> ev;
            auto new_sig = compute_sig(new_IF, ev);
            if (!spec.is_hat1) {
                if (new_sig.sig_pos != 1) continue;
            }

            // 3. det check
            if (spec.is_hat1 && std::abs(new_sig.det) == 1) continue;
            if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) continue;

            // 4. NHC check
            NHCResult nhc = {true, "", {}, {}, {}};
            if (config.check_nhc) {
                nhc = check_nhc(new_IF, config.cc_budget);
                if (!nhc.passes) continue;
            }

            // 5. Gauge enhancement
            if (spec.is_hat1) {
                enhance_hat1_gauge(nhc, new_IF, ext_idx, spec.target_si, config.cc_budget);
            } else if (spec.ext_si == -2 && spec.target_si == -1) {
                enhance_external_m2_gauge(nhc, new_IF, ext_idx, spec.target_si, config.cc_budget);
            }
            if (!nhc.passes) continue;

            // 5.5 Proper c_ext (NHC-aware bridge-level) ≤ 16 — prune over-budget
            std::set<int> ext_set = {ext_idx};
            if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
                continue;

            // 5.6 True (fiber-weighted) c_ext ≤ 16 — REJECT (no recovery)
            {
                double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
                if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) continue;
            }

            // 6. Exact anomaly (with real T)
            AnomalyResult anom = compute_anomaly(nhc, new_sig);

            // Anomaly + reject_nonuni check
            bool is_sugra = true;
            if (config.check_anomaly && anom.H_neutral < 0) is_sugra = false;
            if (is_sugra) {
                std::set<int> h1;
                if (spec.is_hat1) h1.insert(ext_idx);
                if (reject_nonuni(new_IF, new_sig, anom, h1)) is_sugra = false;
            }

            // Record
            Phase1Result r;
            r.spec_id = spec.id;
            r.catalog_id = entry.id;
            r.catalog_type = entry.type;
            r.base_T = entry.T;
            r.target_idx = tidx;
            r.ext_curve_idx = ext_idx;
            r.base_IF = base_IF;
            r.final_IF = new_IF;
            r.nhc = nhc;
            r.anomaly = anom;
            r.sig = new_sig;
            r.eigenvalues = std::move(ev);
            if (is_sugra) {
                push_result(r, new_IF, nhc);
            } else {
                push_nonsugra(r, new_IF, nhc);
            }
        }

        // hat1 multi-target: single hat1 curve → multiple -1 curves (v5 style)
        if (spec.is_hat1 && spec.target_si == -1 && !skip_multi_target) {
            // hat1→(-1): su(8) cc=7.0. Only (-1) with used_cc ≤ 1.0
            double hat1_cc = central_charge(GAUGE_SU8, 1);
            std::vector<int> m1_curves;
            for (int mi = 0; mi < (int)all_m1_curves.size(); mi++)
                if (m1_used_cc[mi] + hat1_cc <= config.cc_budget + 1e-9)
                    m1_curves.push_back(all_m1_curves[mi]);

            int nm1 = (int)m1_curves.size();
            if (nm1 >= 2) {
                int hat1_max = std::min(nm1, 6);
                enumerate_subsets(m1_curves, 2, hat1_max, [&](const std::vector<int>& targets) {
                    Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                    new_IF.block(0, 0, n_base, n_base) = base_IF;
                    new_IF(n_base, n_base) = -1;
                    for (int t : targets) { new_IF(n_base, t) = 1; new_IF(t, n_base) = 1; }
                    int ext_idx = n_base;

                    std::vector<double> ev;
                    auto new_sig = compute_sig(new_IF, ev);
                    if (new_sig.sig_neg > config.T_max) return;

                    // hat1→(-1): unimodular cannot reach Hirzebruch
                    if (std::abs(new_sig.det) == 1) return;
                    if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) return;

                    NHCResult nhc = {true, "", {}, {}, {}};
                    if (config.check_nhc) {
                        nhc = check_nhc(new_IF, config.cc_budget);
                        if (!nhc.passes) return;
                    }
                    enhance_hat1_gauge(nhc, new_IF, ext_idx, -1, config.cc_budget);
                    if (!nhc.passes) return;

                    AnomalyResult anom = compute_anomaly(nhc, new_sig);
                    bool _is_sugra = true;
                    if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                    if (_is_sugra) { std::set<int> h1; h1.insert(ext_idx);
                      if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                    Phase1Result r;
                    r.spec_id = spec.id;
                    r.catalog_id = entry.id;
                    r.catalog_type = entry.type;
                    r.base_T = entry.T;
                    r.target_idx = targets[0];
                    r.ext_curve_idx = ext_idx;
                    r.hat1_multi_targets = targets;
                    r.base_IF = base_IF;
                    r.final_IF = new_IF;
                    r.nhc = nhc;
                    r.anomaly = anom;
                    r.sig = new_sig;
                    r.eigenvalues = std::move(ev);
                    if (_is_sugra) push_result(r, new_IF, nhc);
                    else push_nonsugra(r, new_IF, nhc);
                });
            }
        }

        time_A += t_elapsed_ms(tA); count_A++;
        auto tB = t_start();
        // v6 multi-target: standard external → multiple -1 curves
        // Blowdown method: attach, build extended IF, blowdown, extract level, cc check
        if (!spec.is_hat1 && spec.target_si == -1 && !skip_multi_target) {
            GaugeInfo ext_gauge = gauge_from_si(spec.ext_si);
            if (ext_gauge.dim == 0 && spec.ext_si == -2)
                ext_gauge = GAUGE_SU2;

            // Pre-filter: only (-1) curves with enough cc budget (NHC-aware)
            std::vector<int> m1_curves;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i, i) != -1) continue;
                double used = used_cc_on_minus1_nhc(base_IF, i);
                double min_cc = (ext_gauge.dim > 0) ? central_charge(ext_gauge, 1) : 0;
                if (used + min_cc <= config.cc_budget + 1e-9)
                    m1_curves.push_back(i);
            }

            int nm1 = (int)m1_curves.size();
            // Physically motivated limits per gauge algebra:
            // su2 (-2): max 3 targets, int_num=1 only
            // su3 (-3): max 4 targets
            int max_k, max_tgt;
            if (spec.ext_si == -2) {
                max_k = config.mixed_int_max;
                max_tgt = std::min(nm1, 193);  // NHC pre-filter already limits nm1
            } else if (spec.ext_si == -3) {
                max_k = config.mixed_int_max;
                max_tgt = std::min(nm1, 193);
            } else {
                // so8 and higher: k=1 only (k=2 would be so16 etc, handled separately)
                max_k = 1;
                max_tgt = std::min(nm1, 193);
            }

            if (max_tgt >= 2) {
                enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](const std::vector<int>& targets) {
                    enumerate_int_nums((int)targets.size(), max_k, [&](const std::vector<int>& int_nums) {
                        // Total ext cc check
                        if (total_ext_cc_multi(spec.ext_si, int_nums, spec.is_hat1, spec.target_si) > g_max_ext_cc + 1e-9) return;

                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                        new_IF.block(0, 0, n_base, n_base) = base_IF;
                        new_IF(n_base, n_base) = spec.ext_si;
                        for (int i = 0; i < (int)targets.size(); i++) {
                            new_IF(n_base, targets[i]) = int_nums[i];
                            new_IF(targets[i], n_base) = int_nums[i];
                        }
                        int ext_idx = n_base;

                        // Quick cc check BEFORE expensive LDLT
                        if (!quick_cc_check(new_IF, config.cc_budget)) return;

                        // Fast signature pre-filter (LDLT, no eigenvalues)
                        if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;


                        NHCResult nhc = {true, "", {}, {}, {}};
                        if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }
                        if (spec.ext_si == -2) enhance_external_m2_gauge(nhc, new_IF, ext_idx, -1, config.cc_budget);
                        if (!nhc.passes) return;

                        // Full eigensolve for passing candidates (needed for anomaly, det, dedup)
                        std::vector<double> ev;
                        auto new_sig = compute_sig(new_IF, ev);
                        if (new_sig.sig_pos != 1 || new_sig.sig_neg > config.T_max) return;

                        AnomalyResult anom = compute_anomaly(nhc, new_sig);
                        bool _is_sugra = true;
                        if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                        if (_is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                        Phase1Result r;
                        r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                        r.base_T = entry.T; r.target_idx = targets[0]; r.ext_curve_idx = ext_idx;
                        r.hat1_multi_targets = targets; r.base_IF = base_IF; r.final_IF = new_IF;
                        r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                        if (_is_sugra) push_result(r, new_IF, nhc);
                        else push_nonsugra(r, new_IF, nhc);
                    });
                });
            }
        }

        time_B += t_elapsed_ms(tB); count_B++;
        auto tD = t_start();
        // --- Mixed multi-target: su2 → (-1) and (-3) simultaneously ---
        // Skip in largeint mode (single-target only); always run in largeint_mixed mode
        if (spec.ext_si == -2 && spec.target_si == -1 && !spec.is_hat1
            && (!skip_multi_target || config.largeint_mixed)
            && !config.largeint_single) {
            // Collect (-1) with cc budget (NHC-aware) and (-3) curves
            std::vector<int> m1_curves, m3_curves;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i, i) == -1) {
                    double used = used_cc_on_minus1_nhc(base_IF, i);
                    if (used + central_charge(GAUGE_SU2, 1) <= config.cc_budget + 1e-9)
                        m1_curves.push_back(i);
                }
                else if (base_IF(i, i) == -3) m3_curves.push_back(i);
            }

            if (!m1_curves.empty() && !m3_curves.empty()) {
                // For each (-3) target, combine with subsets of (-1) targets
                // (-1) targets can have int=1~mixed_int_max, (-3) always int=1
                int int_max = config.mixed_int_max;
                for (int m3 : m3_curves) {
                    int max_m1 = std::min((int)m1_curves.size(), 193);
                    enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                        enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                        new_IF.block(0, 0, n_base, n_base) = base_IF;
                        new_IF(n_base, n_base) = -2;  // su2 external
                        // Connect to (-3) with int=1
                        new_IF(n_base, m3) = 1; new_IF(m3, n_base) = 1;
                        // Connect to (-1)s with varying int
                        for (int i = 0; i < (int)m1_targets.size(); i++) {
                            new_IF(n_base, m1_targets[i]) = int_nums[i];
                            new_IF(m1_targets[i], n_base) = int_nums[i];
                        }
                        int ext_idx = n_base;

                        if (!quick_cc_check(new_IF, config.cc_budget)) return;

                        if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;

                        NHCResult nhc = {true, "", {}, {}, {}};
                        if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }
                        // No su2 gauge enhancement here — su2 is forming NHC with (-3)
                        if (!nhc.passes) return;

                        std::vector<double> ev;
                        auto new_sig = compute_sig(new_IF, ev);
                        if (new_sig.sig_pos != 1 || new_sig.sig_neg > config.T_max) return;
                        if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) return;

                        AnomalyResult anom = compute_anomaly(nhc, new_sig);
                        bool _is_sugra = true;
                        if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                        if (_is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                        Phase1Result r;
                        r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                        r.base_T = entry.T; r.target_idx = m3; r.ext_curve_idx = ext_idx;
                        std::vector<int> all_targets = {m3};
                        all_targets.insert(all_targets.end(), m1_targets.begin(), m1_targets.end());
                        r.hat1_multi_targets = all_targets;
                        r.base_IF = base_IF; r.final_IF = new_IF;
                        r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                        if (_is_sugra) push_result(r, new_IF, nhc);
                        else push_nonsugra(r, new_IF, nhc);
                        });
                    });
                }
            }
        }

        time_D += t_elapsed_ms(tD); count_D++;
        auto tE = t_start();
        // --- Mixed multi-target: su2 → (-4) with int=2, and (-1) simultaneously ---
        if (spec.ext_si == -2 && spec.target_si == -1 && !spec.is_hat1 && !skip_multi_target
            && !config.largeint_single) {
            std::vector<int> m1_curves, m4_curves;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i, i) == -1) {
                    double used = used_cc_on_minus1_nhc(base_IF, i);
                    if (used + central_charge(GAUGE_SU2, 1) <= config.cc_budget + 1e-9)
                        m1_curves.push_back(i);
                }
                else if (base_IF(i, i) == -4) m4_curves.push_back(i);
            }

            if (!m1_curves.empty() && !m4_curves.empty()) {
                for (int m4 : m4_curves) {
                    int max_m1 = std::min((int)m1_curves.size(), 193);
                    enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                        enumerate_int_nums((int)m1_targets.size(), config.mixed_int_max, [&](const std::vector<int>& int_nums) {
                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                        new_IF.block(0, 0, n_base, n_base) = base_IF;
                        new_IF(n_base, n_base) = -2;
                        // Connect to (-4) with int=2
                        new_IF(n_base, m4) = 2; new_IF(m4, n_base) = 2;
                        // Connect to (-1)s with varying int
                        for (int i = 0; i < (int)m1_targets.size(); i++) {
                            new_IF(n_base, m1_targets[i]) = int_nums[i];
                            new_IF(m1_targets[i], n_base) = int_nums[i];
                        }
                        int ext_idx = n_base;

                        if (!quick_cc_check(new_IF, config.cc_budget)) return;

                        if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;

                        NHCResult nhc = {true, "", {}, {}, {}};
                        if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }
                        if (!nhc.passes) return;

                        std::vector<double> ev;
                        auto new_sig = compute_sig(new_IF, ev);
                        if (new_sig.sig_pos != 1 || new_sig.sig_neg > config.T_max) return;
                        if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) return;

                        AnomalyResult anom = compute_anomaly(nhc, new_sig);
                        bool _is_sugra = true;
                        if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                        if (_is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                        Phase1Result r;
                        r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                        r.base_T = entry.T; r.target_idx = m4; r.ext_curve_idx = ext_idx;
                        std::vector<int> all_targets = {m4};
                        all_targets.insert(all_targets.end(), m1_targets.begin(), m1_targets.end());
                        r.hat1_multi_targets = all_targets;
                        r.base_IF = base_IF; r.final_IF = new_IF;
                        r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                        if (_is_sugra) push_result(r, new_IF, nhc);
                        else push_nonsugra(r, new_IF, nhc);
                        });
                    });
                }
            }
        }

        time_E += t_elapsed_ms(tE); count_E++;
        auto tF = t_start();
        // --- Mixed multi-target: so16 → (-2) with int=2, and (-1) simultaneously ---
        // --- Generic mixed helper lambda ---
        auto try_mixed = [&](int ext_si, int special_curve, int special_int,
                             const std::vector<int>& m1c, const std::string& mtag) {
            if (m1c.empty()) return;
            int max_m1 = std::min((int)m1c.size(), 4);
            enumerate_subsets_v2(m1c, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                new_IF.block(0, 0, n_base, n_base) = base_IF;
                new_IF(n_base, n_base) = ext_si;
                new_IF(n_base, special_curve) = special_int;
                new_IF(special_curve, n_base) = special_int;
                for (int t : m1_targets) {
                    new_IF(n_base, t) = 1; new_IF(t, n_base) = 1;
                }
                int ext_idx = n_base;

                if (!quick_cc_check(new_IF, config.cc_budget)) return;
                if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;

                NHCResult nhc = {true, "", {}, {}, {}};
                if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }

                std::vector<double> ev;
                auto new_sig = compute_sig(new_IF, ev);
                if (new_sig.sig_pos != 1) return;
                if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) return;

                AnomalyResult anom = compute_anomaly(nhc, new_sig);
                bool is_sugra = true;
                if (config.check_anomaly && anom.H_neutral < 0) is_sugra = false;
                if (is_sugra) {
                    std::set<int> h1;
                    if (reject_nonuni(new_IF, new_sig, anom, h1)) is_sugra = false;
                }

                Phase1Result r;
                r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                r.base_T = entry.T; r.target_idx = special_curve; r.ext_curve_idx = ext_idx;
                std::vector<int> all_tgts = {special_curve};
                all_tgts.insert(all_tgts.end(), m1_targets.begin(), m1_targets.end());
                r.hat1_multi_targets = all_tgts;
                r.base_IF = base_IF; r.final_IF = new_IF;
                r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                if (is_sugra) push_result(r, new_IF, nhc);
                else push_nonsugra(r, new_IF, nhc);
            });
        };

        // so16n2mix: (-4) → (-2) int=2 + (-1)s
        if (spec.ext_si == -4 && spec.target_si == -1 && !spec.is_hat1 && !skip_multi_target) {
            std::vector<int> m1c, m2c;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i,i) == -1) {
                    double used = used_cc_on_minus1_nhc(base_IF, i);
                    if (used + central_charge(gauge_from_si(-4), 1) <= config.cc_budget + 1e-9)
                        m1c.push_back(i);
                } else if (base_IF(i,i) == -2) m2c.push_back(i);
            }
            for (int m2 : m2c) try_mixed(-4, m2, 2, m1c, "so16n2mix");
        }

        // su3n2mix: (-3) → (-2) int=1 + (-1)s
        if (spec.ext_si == -3 && spec.target_si == -1 && !spec.is_hat1 && !skip_multi_target) {
            std::vector<int> m1c, m2c;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i,i) == -1) {
                    double used = used_cc_on_minus1_nhc(base_IF, i);
                    if (used + central_charge(GAUGE_G2, 1) <= config.cc_budget + 1e-9)
                        m1c.push_back(i);
                } else if (base_IF(i,i) == -2) m2c.push_back(i);
            }
            for (int m2 : m2c) try_mixed(-3, m2, 1, m1c, "su3n2mix");
        }

        // so7: (-3) → (-2)+(-2) (no (-1) connection)
        if (spec.ext_si == -3 && spec.target_si == -1 && !spec.is_hat1 && !skip_multi_target) {
            std::vector<int> m2c;
            for (int i = 0; i < n_base; i++)
                if (base_IF(i,i) == -2) m2c.push_back(i);
            if ((int)m2c.size() >= 2) {
                for (int a = 0; a < (int)m2c.size(); a++) {
                    for (int b = a+1; b < (int)m2c.size(); b++) {
                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base+1, n_base+1);
                        new_IF.block(0,0,n_base,n_base) = base_IF;
                        new_IF(n_base, n_base) = -3;
                        new_IF(n_base, m2c[a]) = 1; new_IF(m2c[a], n_base) = 1;
                        new_IF(n_base, m2c[b]) = 1; new_IF(m2c[b], n_base) = 1;
                        int ext_idx = n_base;

                        if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) continue;
                        NHCResult nhc = {true, "", {}, {}, {}};
                        if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) continue; }

                        std::vector<double> ev;
                        auto new_sig = compute_sig(new_IF, ev);
                        if (new_sig.sig_pos != 1) continue;

                        AnomalyResult anom = compute_anomaly(nhc, new_sig);
                        bool is_sugra = true;
                        if (config.check_anomaly && anom.H_neutral < 0) is_sugra = false;
                        if (is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) is_sugra = false; }

                        Phase1Result r;
                        r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                        r.base_T = entry.T; r.target_idx = m2c[a]; r.ext_curve_idx = ext_idx;
                        r.hat1_multi_targets = {m2c[a], m2c[b]};
                        r.base_IF = base_IF; r.final_IF = new_IF;
                        r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                        if (is_sugra) push_result(r, new_IF, nhc);
                        else push_nonsugra(r, new_IF, nhc);
                    }
                }
            }
        }

        // so7mix: (-3) → (-2)+(-2) + (-1)s
        if (spec.ext_si == -3 && spec.target_si == -1 && !spec.is_hat1 && !skip_multi_target) {
            std::vector<int> m1c, m2c;
            for (int i = 0; i < n_base; i++) {
                if (base_IF(i,i) == -1) {
                    double used = used_cc_on_minus1_nhc(base_IF, i);
                    if (used + central_charge(GAUGE_SO7, 1) <= config.cc_budget + 1e-9)
                        m1c.push_back(i);
                } else if (base_IF(i,i) == -2) m2c.push_back(i);
            }
            if ((int)m2c.size() >= 2 && !m1c.empty()) {
                for (int a = 0; a < (int)m2c.size(); a++) {
                    for (int b = a+1; b < (int)m2c.size(); b++) {
                        int max_m1 = std::min((int)m1c.size(), 4);
                        enumerate_subsets_v2(m1c, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base+1, n_base+1);
                            new_IF.block(0,0,n_base,n_base) = base_IF;
                            new_IF(n_base, n_base) = -3;
                            new_IF(n_base, m2c[a]) = 1; new_IF(m2c[a], n_base) = 1;
                            new_IF(n_base, m2c[b]) = 1; new_IF(m2c[b], n_base) = 1;
                            for (int t : m1_targets) { new_IF(n_base, t) = 1; new_IF(t, n_base) = 1; }
                            int ext_idx = n_base;

                            if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;
                            NHCResult nhc = {true, "", {}, {}, {}};
                            if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }

                            std::vector<double> ev;
                            auto new_sig = compute_sig(new_IF, ev);
                            if (new_sig.sig_pos != 1) return;

                            AnomalyResult anom = compute_anomaly(nhc, new_sig);
                            bool is_sugra = true;
                            if (config.check_anomaly && anom.H_neutral < 0) is_sugra = false;
                            if (is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) is_sugra = false; }

                            Phase1Result r;
                            r.spec_id = spec.id; r.catalog_id = entry.id; r.catalog_type = entry.type;
                            r.base_T = entry.T; r.target_idx = m2c[a]; r.ext_curve_idx = ext_idx;
                            r.hat1_multi_targets = {m2c[a], m2c[b]};
                            r.hat1_multi_targets.insert(r.hat1_multi_targets.end(), m1_targets.begin(), m1_targets.end());
                            r.base_IF = base_IF; r.final_IF = new_IF;
                            r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                            if (is_sugra) push_result(r, new_IF, nhc);
                            else push_nonsugra(r, new_IF, nhc);
                        });
                    }
                }
            }
        }
        time_F += t_elapsed_ms(tF); count_F++;
    }

    // --- Dummy LSTs ---
    if (config.include_dummy) {
        auto try_dummy = [&](const Eigen::MatrixXi& base_IF, int id, const std::string& type) {
            int n_base = base_IF.rows();
            auto base_sig = compute_sig(base_IF);
            if (base_sig.sig_neg + 1 > config.T_max) return;
            g_current_base_T = n_base;  // for --use-lst-T mode: LST curve count

            std::vector<int> target_curves;
            for (int i = 0; i < n_base; i++)
                if (base_IF(i, i) == spec.target_si)
                    target_curves.push_back(i);

            for (int tidx : target_curves) {
                if (spec.target_si == -1 && !spec.is_hat1) {
                    GaugeInfo g = gauge_from_si(spec.ext_si);
                    if (g.dim > 0) {
                        double rem = config.cc_budget - used_cc_on_minus1(base_IF, tidx);
                        if (central_charge(g, spec.int_num) > rem + 1e-9) continue;
                    }
                    if (spec.ext_si == -2) {
                        double rem = config.cc_budget - used_cc_on_minus1(base_IF, tidx);
                        if (central_charge(GAUGE_SU2, spec.int_num) > rem + 1e-9) continue;
                    }
                }

                Eigen::MatrixXi new_IF = attach_curve(base_IF, tidx, spec.ext_si, spec.int_num);
                int ext_idx = new_IF.rows() - 1;

                std::vector<double> ev;
                auto new_sig = compute_sig(new_IF, ev);
                if (spec.is_hat1) {
                    if (new_sig.sig_neg > config.T_max) continue;
                } else {
                    if (new_sig.sig_pos != 1 || new_sig.sig_neg > config.T_max) continue;
                }

                // hat1: unimodular cannot reach Hirzebruch
                if (spec.is_hat1 && std::abs(new_sig.det) == 1) continue;
                if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) continue;

                NHCResult nhc = {true, "", {}, {}, {}};
                if (config.check_nhc) {
                    nhc = check_nhc(new_IF, config.cc_budget);
                    if (!nhc.passes) continue;
                }
                if (spec.is_hat1)
                    enhance_hat1_gauge(nhc, new_IF, ext_idx, spec.target_si, config.cc_budget);
                else if (spec.ext_si == -2 && spec.target_si == -1)
                    enhance_external_m2_gauge(nhc, new_IF, ext_idx, spec.target_si, config.cc_budget);
                if (!nhc.passes) continue;

                AnomalyResult anom = compute_anomaly(nhc, new_sig);
                bool _is_sugra = true;
                if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                if (_is_sugra) {
                    std::set<int> h1;
                    if (spec.is_hat1) h1.insert(ext_idx);
                    if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false;
                }

                Phase1Result r;
                r.spec_id = spec.id;
                r.catalog_id = id;
                r.catalog_type = type;
                r.base_T = base_sig.sig_neg;
                r.target_idx = tidx;
                r.ext_curve_idx = ext_idx;
                r.base_IF = base_IF;
                r.final_IF = new_IF;
                r.nhc = nhc;
                r.anomaly = anom;
                r.sig = new_sig;
                r.eigenvalues = std::move(ev);
                if (_is_sugra) push_result(r, new_IF, nhc);
                else push_nonsugra(r, new_IF, nhc);
            }

            // hat1 multi-target on dummies
            if (spec.is_hat1 && spec.target_si == -1 && !skip_multi_target) {
                std::vector<int> m1_curves;
                for (int i = 0; i < n_base; i++)
                    if (base_IF(i, i) == -1) m1_curves.push_back(i);
                int nm1 = (int)m1_curves.size();
                if (nm1 >= 2) {
                    int hat1_max = std::min(nm1, 6);
                    enumerate_subsets(m1_curves, 2, hat1_max, [&](const std::vector<int>& targets) {
                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                        new_IF.block(0, 0, n_base, n_base) = base_IF;
                        new_IF(n_base, n_base) = -1;
                        for (int t : targets) { new_IF(n_base, t) = 1; new_IF(t, n_base) = 1; }
                        int ext_idx = n_base;

                        std::vector<double> ev;
                        auto new_sig = compute_sig(new_IF, ev);
                        if (new_sig.sig_neg > config.T_max) return;
                        if (std::abs(new_sig.det) == 1) return; // unimodular hat1→(-1) rejected
                        if (config.check_det_perfect_square && !is_perfect_square(std::abs(new_sig.det))) return;
                        NHCResult nhc = {true, "", {}, {}, {}};
                        if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }
                        enhance_hat1_gauge(nhc, new_IF, ext_idx, -1, config.cc_budget);
                        if (!nhc.passes) return;
                        AnomalyResult anom = compute_anomaly(nhc, new_sig);
                        bool _is_sugra = true;
                        if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                        if (_is_sugra) { std::set<int> h1; h1.insert(ext_idx); if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                        Phase1Result r;
                        r.spec_id = spec.id; r.catalog_id = id; r.catalog_type = type;
                        r.base_T = base_sig.sig_neg; r.target_idx = targets[0]; r.ext_curve_idx = ext_idx;
                        r.hat1_multi_targets = targets; r.base_IF = base_IF; r.final_IF = new_IF;
                        r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                        if (_is_sugra) push_result(r, new_IF, nhc);
                        else push_nonsugra(r, new_IF, nhc);
                    });
                }
            }

            // v6 multi-target on dummies: standard external → multiple -1 curves
            // Blowdown method
            if (!spec.is_hat1 && spec.target_si == -1 && !skip_multi_target) {
                GaugeInfo ext_gauge = gauge_from_si(spec.ext_si);
                if (ext_gauge.dim == 0 && spec.ext_si == -2)
                    ext_gauge = GAUGE_SU2;

                std::vector<int> m1_curves;
                for (int i = 0; i < n_base; i++)
                    if (base_IF(i, i) == -1) m1_curves.push_back(i);

                int nm1 = (int)m1_curves.size();
                int max_k, max_tgt;
                if (spec.ext_si == -2) {
                    max_k = 5;
                    max_tgt = std::min(nm1, 3);
                } else if (spec.ext_si == -3) {
                    max_k = 4;
                    max_tgt = std::min(nm1, 3);
                } else {
                    max_k = std::min(max_level_for_cc(ext_gauge, 20.0), 9);
                    max_tgt = std::min(nm1, 6);
                }

                if (max_tgt >= 2) {
                    enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](const std::vector<int>& targets) {
                        enumerate_int_nums((int)targets.size(), max_k, [&](const std::vector<int>& int_nums) {
                            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n_base + 1, n_base + 1);
                            new_IF.block(0, 0, n_base, n_base) = base_IF;
                            new_IF(n_base, n_base) = spec.ext_si;
                            for (int i = 0; i < (int)targets.size(); i++) {
                                new_IF(n_base, targets[i]) = int_nums[i];
                                new_IF(targets[i], n_base) = int_nums[i];
                            }
                            int ext_idx = n_base;

                            // Fast signature pre-filter (LDLT, no eigenvalues)
                            if (!quick_cc_check(new_IF, config.cc_budget)) return;

                            if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;

                            NHCResult nhc = {true, "", {}, {}, {}};
                            if (config.check_nhc) { nhc = check_nhc(new_IF, config.cc_budget); if (!nhc.passes) return; }
                            if (spec.ext_si == -2) enhance_external_m2_gauge(nhc, new_IF, ext_idx, -1, config.cc_budget);
                            if (!nhc.passes) return;

                            // Full eigensolve for passing candidates
                            std::vector<double> ev;
                            auto new_sig = compute_sig(new_IF, ev);
                            if (new_sig.sig_pos != 1 || new_sig.sig_neg > config.T_max) return;
                            AnomalyResult anom = compute_anomaly(nhc, new_sig);
                            bool _is_sugra = true;
                            if (config.check_anomaly && anom.H_neutral < 0) _is_sugra = false;
                            if (_is_sugra) { std::set<int> h1; if (reject_nonuni(new_IF, new_sig, anom, h1)) _is_sugra = false; }

                            Phase1Result r;
                            r.spec_id = spec.id; r.catalog_id = id; r.catalog_type = type;
                            r.base_T = base_sig.sig_neg; r.target_idx = targets[0]; r.ext_curve_idx = ext_idx;
                            r.hat1_multi_targets = targets; r.base_IF = base_IF; r.final_IF = new_IF;
                            r.nhc = nhc; r.anomaly = anom; r.sig = new_sig; r.eigenvalues = std::move(ev);
                            if (_is_sugra) push_result(r, new_IF, nhc);
                            else push_nonsugra(r, new_IF, nhc);
                        });
                    });
                }
            }
        };

        int a_min = std::max(2, config.catalog_T_min + 1);
        int a_max = config.catalog_T_max + 1;
        for (int nc = a_min; nc <= a_max; nc++)
            try_dummy(build_dummy_A(nc), -nc, "DM:A(" + std::to_string(nc) + ")");

        int d_min = std::max(3, config.catalog_T_min + 1);
        int d_max = config.catalog_T_max + 1;
        for (int nc = d_min; nc <= d_max; nc++)
            try_dummy(build_dummy_D(nc), -(1000 + nc), "DM:D(" + std::to_string(nc) + ")");
    }

    return {std::move(results), std::move(nonsugra_results)};
}

// ============================================================================
// PART 7: Write .cat file for one spec (stacking-ready format)
// ============================================================================
//
// Format:
//   # header
//   ENTRY <result_id>
//   SPEC <spec_id> <tag> <ext_si> <target_si> <int_num> <hat1>
//   BASE <catalog_id> <type> <base_T> <n_base>
//   ATTACH <target_idx> <ext_curve_idx> [multi:<idx1,idx2,...>]
//   PHYSICS <T> <H_c> <V> <H_n> <det> <sig+> <sig-> <sig0>
//   IF <n>
//   <row0>
//   <row1>
//   ...
//   REMAINING <n_targets>
//   <curve_idx> <self_int> <available_cc>
//   ...
//   END

inline void write_phase1_cat(const std::string& filename,
                              const ExternalSpec& spec,
                              const std::vector<Phase1Result>& results) {
    std::ofstream out(filename);

    out << "# Phase1 SUGRA catalog: " << spec.tag << "\n";
    out << "# External: " << spec.label
        << " (ext_si=" << spec.ext_si
        << ", target_si=" << spec.target_si
        << ", int_num=" << spec.int_num
        << ", hat1=" << spec.is_hat1
        << ", cc=" << std::fixed << std::setprecision(2) << spec.cc << ")\n";
    out << "# Total: " << results.size() << " bases\n";
    out << "#\n";

    for (int ri = 0; ri < (int)results.size(); ri++) {
        auto& r = results[ri];
        int n = r.final_IF.rows();

        out << "ENTRY " << ri << "\n";
        out << "SPEC " << r.spec_id << " " << spec.tag
            << " " << spec.ext_si << " " << spec.target_si
            << " " << spec.int_num << " " << (spec.is_hat1 ? 1 : 0) << "\n";
        out << "BASE " << r.catalog_id << " " << r.catalog_type
            << " " << r.base_T << " " << r.base_IF.rows() << "\n";

        out << "ATTACH " << r.target_idx << " " << r.ext_curve_idx;
        if (!r.hat1_multi_targets.empty() && (int)r.hat1_multi_targets.size() > 1) {
            out << " multi:";
            for (int i = 0; i < (int)r.hat1_multi_targets.size(); i++) {
                if (i) out << ",";
                out << r.hat1_multi_targets[i];
            }
        }
        out << "\n";

        out << "PHYSICS " << r.anomaly.T << " " << r.anomaly.H_charged
            << " " << r.anomaly.V << " " << r.anomaly.H_neutral
            << " " << r.sig.det << " " << r.sig.sig_pos
            << " " << r.sig.sig_neg << " " << r.sig.sig_zero << "\n";

        out << "IF " << n << "\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j) out << " ";
                out << r.final_IF(i, j);
            }
            out << "\n";
        }

        // Remaining targets for stacking
        auto& rem = r.remaining_targets;
        // Filter to only useful targets (exclude the external itself)
        std::vector<Phase1Result::AvailTarget> useful;
        for (auto& t : rem) {
            if (t.curve_idx == r.ext_curve_idx) continue;
            // Keep -1 curves with remaining cc > 0, and -2/-3/-4 curves
            if (t.self_int == -1 && t.available_cc > 0.1)
                useful.push_back(t);
            else if (t.self_int <= -2)
                useful.push_back(t);
        }
        out << "REMAINING " << useful.size() << "\n";
        for (auto& t : useful) {
            out << t.curve_idx << " " << t.self_int << " "
                << std::fixed << std::setprecision(2) << t.available_cc << "\n";
        }
        out << "END\n\n";
    }
    out.close();
}

// ============================================================================
// PART 8: Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <catalog> [T_max=193] [T_min=0] [--det-sq]\n"
                  << "       " << argv[0]
                  << " <catalog> T_max T_min --largeint k_min k_max\n";
        return 1;
    }

    std::string catalog_file = argv[1];
    int T_max = (argc > 2) ? std::atoi(argv[2]) : 193;
    int T_min = (argc > 3) ? std::atoi(argv[3]) : 0;
    bool det_sq_mode = false;
    bool no_su2 = false;
    bool largeint_mode = false;
    bool largeint_mixed_mode = false;
    int largeint_kmin = 1, largeint_kmax = 30;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--det-sq") det_sq_mode = true;
        if (std::string(argv[i]) == "--largeint") {
            largeint_mode = true;
            if (i + 1 < argc) largeint_kmin = std::atoi(argv[i + 1]);
            if (i + 2 < argc) largeint_kmax = std::atoi(argv[i + 2]);
        }
        if (std::string(argv[i]) == "--largeint-mixed") {
            largeint_mixed_mode = true;
            if (i + 1 < argc) largeint_kmin = std::atoi(argv[i + 1]);
            if (i + 2 < argc) largeint_kmax = std::atoi(argv[i + 2]);
        }
        if (std::string(argv[i]) == "--mix-k" && i + 1 < argc) {
            largeint_kmax = std::atoi(argv[i + 1]);  // reuse for mixed_int_max override
        }
        if (std::string(argv[i]) == "--use-lst-T") g_use_lst_T = true;
        if (std::string(argv[i]) == "--no-su2") no_su2 = true;
        if (std::string(argv[i]) == "--save-nonsugra") {}  // handled below in config
    }
    std::string suffix = "_T" + std::to_string(T_min) + "_" + std::to_string(T_max);

    std::cout << "=== gen_sugra_phase1 ===\n";
    std::cout << "Catalog: " << catalog_file << "\n";
    std::cout << "T range: [" << T_min << ", " << T_max << "]"
              << (det_sq_mode ? "  [|det|=n² filter ON]" : "")
              << (largeint_mode ? "  [LARGEINT su2 k=" + std::to_string(largeint_kmin)
                                  + "~" + std::to_string(largeint_kmax) + "]" : "")
              << (largeint_mixed_mode ? "  [LARGEINT-MIXED su2→(-3)+(-1) k=" + std::to_string(largeint_kmin)
                                  + "~" + std::to_string(largeint_kmax) + "]" : "")
              << "\n\n";

    // Load catalog
    auto catalog = load_catalog(catalog_file);
    std::cout << "Catalog loaded: " << catalog.size() << " entries\n";

    // Build spec table
    std::vector<ExternalSpec> all_specs;
    if (largeint_mode) {
        int id = 0;
        for (int k = largeint_kmin; k <= largeint_kmax; k++) {
            double cc = compute_ext_cc(-2, k, false, -1);
            all_specs.push_back({id++, -2, -1, k, false, cc,
                "su(2) k=" + std::to_string(k), "su2"});
        }
    } else if (largeint_mixed_mode) {
        // Single spec: su2→(-1) k=1 as placeholder; mixed multi-target does the real work
        double cc = compute_ext_cc(-2, 1, false, -1);
        all_specs.push_back({0, -2, -1, 1, false, cc, "su2 mixed", "su2"});
    } else {
        all_specs = build_all_specs();
    }
    if (no_su2) {
        all_specs.erase(std::remove_if(all_specs.begin(), all_specs.end(),
            [](const ExternalSpec& s) { return s.ext_si == -2 && s.target_si == -1 && !s.is_hat1; }),
            all_specs.end());
    }
    std::cout << "External specs: " << all_specs.size() << " types"
              << (no_su2 ? " [no-su2]" : "") << "\n\n";

    // Config — following v2: T_max = 193 for SUGRA bound, catalog range separate
    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;
    config.check_determinant = false;
    config.check_det_perfect_square = det_sq_mode;
    config.catalog_T_min = T_min;
    config.catalog_T_max = T_max;
    config.include_nn = true;
    config.include_nk = true;
    config.include_dummy = true;
    config.verbose = false;
    config.largeint_single = largeint_mode;
    config.largeint_mixed = largeint_mixed_mode;
    // --mix-k overrides mixed_int_max; default 5 for normal mode
    int default_mix_k = 5;
    for (int i = 1; i < argc; i++)
        if (std::string(argv[i]) == "--mix-k" && i + 1 < argc)
            default_mix_k = std::atoi(argv[i + 1]);
    config.mixed_int_max = (largeint_mode || largeint_mixed_mode) ? largeint_kmax : default_mix_k;
    for (int i = 1; i < argc; i++)
        if (std::string(argv[i]) == "--save-nonsugra") config.save_nonsugra = true;

    // Precompute IF and signature for all catalog entries (shared across all specs)
    std::vector<CachedIF> cached_ifs(catalog.size());
    int n_cached = 0;
    for (size_t i = 0; i < catalog.size(); i++) {
        auto& entry = catalog[i];
        cached_ifs[i].valid = false;
        if (entry.is_nn() && !config.include_nn) continue;
        if (!entry.is_nn() && !config.include_nk) continue;
        if (entry.T < config.catalog_T_min || entry.T > config.catalog_T_max) continue;
        cached_ifs[i].IF = reconstruct_IF(entry);
        if (cached_ifs[i].IF.rows() == 0) continue;
        cached_ifs[i].sig = compute_sig(cached_ifs[i].IF);
        cached_ifs[i].valid = true;
        n_cached++;
    }
    std::cout << "Precomputed IFs: " << n_cached << " entries\n\n";

    // Output directory (--out <dir> to override)
    std::string dir = largeint_mode ? "cat_largeint" :
                      largeint_mixed_mode ? "cat_largeint_mixed" : "cat_1ext";
    for (int i = 1; i < argc; i++)
        if (std::string(argv[i]) == "--out" && i + 1 < argc) dir = argv[i + 1];
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    // Group specs by (ext_si, target_si) to cache multi-target results
    // Multi-target enumeration depends only on (ext_si, target_si), not int_num,
    // so specs like su2_k1..k9 produce identical multi-target IFs.
    // We cache multi-target results from the first spec and reuse for subsequent ones.
    struct MultiTargetKey {
        int ext_si, target_si;
        bool is_hat1;
        bool operator<(const MultiTargetKey& o) const {
            if (ext_si != o.ext_si) return ext_si < o.ext_si;
            if (target_si != o.target_si) return target_si < o.target_si;
            return is_hat1 < o.is_hat1;
        }
    };
    // Cache: for each (ext_si, target_si, is_hat1), store the multi-target results
    // from the first spec that computed them
    std::map<MultiTargetKey, std::vector<Phase1Result>> multi_target_cache;
    std::set<MultiTargetKey> multi_target_computed;

    // Group specs by tag → accumulate results, write once per tag
    std::map<std::string, std::vector<Phase1Result>> results_by_tag;
    std::map<std::string, std::vector<Phase1Result>> nonsugra_by_tag;
    std::map<std::string, const ExternalSpec*> tag_to_spec;  // representative spec per tag
    std::map<std::string, std::set<std::string>> seen_by_tag;  // shared dedup per tag

    int total_bases = 0;

    for (auto& spec : all_specs) {
        std::cout << "[" << spec.id << "] " << spec.label
                  << " (ext=" << spec.ext_si << "→" << spec.target_si
                  << ", #=" << spec.int_num << ")";
        std::cout.flush();

        MultiTargetKey mtk{spec.ext_si, spec.target_si, spec.is_hat1};
        bool is_first_in_group = (multi_target_computed.find(mtk) == multi_target_computed.end());

        // Generate with early dedup — shared seen set per tag for cross-k dedup
        auto& seen = seen_by_tag[spec.tag];
        auto spec_result = generate_single_spec(spec, catalog, config, &seen,
                                         is_first_in_group ? nullptr : (const void*)&mtk,
                                         &cached_ifs);
        auto& raw = spec_result.sugra;

        std::cout << " → " << raw.size() << " raw, total "
                  << (results_by_tag[spec.tag].size() + raw.size()) << "\n";

        // Accumulate non-SUGRA results
        if (config.save_nonsugra) {
            auto& ns_group = nonsugra_by_tag[spec.tag];
            for (auto& r : spec_result.nonsugra) ns_group.push_back(std::move(r));
        }

        // If this is the first spec in the group, cache its multi-target results
        if (is_first_in_group && spec.target_si == -1) {
            multi_target_cache[mtk].clear();
            for (auto& r : raw) {
                if (!r.hat1_multi_targets.empty() && (int)r.hat1_multi_targets.size() >= 2) {
                    multi_target_cache[mtk].push_back(r); // copy for reuse
                }
            }
            multi_target_computed.insert(mtk);
        } else if (!is_first_in_group && spec.target_si == -1) {
            // Reuse cached multi-target results with updated spec_id
            auto it = multi_target_cache.find(mtk);
            if (it != multi_target_cache.end()) {
                for (auto cached_r : it->second) {
                    cached_r.spec_id = spec.id;
                    // Check early dedup
                    std::string extra = spec.tag;
                    if (spec.is_hat1)
                        extra += "|T:" + std::to_string(spec.target_si)
                               + "|N:" + std::to_string(cached_r.hat1_multi_targets.size());
                    // Include catalog_id (LST source) so different build paths are not merged
                    extra += "|cid:" + std::to_string(cached_r.catalog_id);
                    std::string key = cached_r.eigenvalues.empty()
                        ? dedup_key_v2(cached_r.final_IF, extra)
                        : dedup_key_v2_cached(cached_r.final_IF, cached_r.eigenvalues, extra);
                    if (seen.insert(key).second)
                        raw.push_back(std::move(cached_r));
                }
            }
            multi_target_computed.insert(mtk);
        }

        std::cout << " → " << raw.size() << " raw";

        // Accumulate into tag group
        auto& group = results_by_tag[spec.tag];
        for (auto& r : raw) group.push_back(std::move(r));
        if (tag_to_spec.find(spec.tag) == tag_to_spec.end())
            tag_to_spec[spec.tag] = &spec;

        std::cout << ", total " << group.size() << "\n";
    }

    // Deduplicate within each tag group and write
    std::ofstream index(dir + "/INDEX.txt");
    index << "# Phase1 SUGRA Catalog Index\n";
    index << "# T range: [" << T_min << ", " << T_max << "]\n";
    index << "#\n";
    index << "# tag | n_dedup | filename\n\n";

    for (auto& [tag, group] : results_by_tag) {
        // No additional dedup here — early dedup in generate_single_spec
        // already uses shared 'seen' set across specs with same tag.
        auto& deduped = group;

        // Sort
        std::sort(deduped.begin(), deduped.end(),
                  [](const Phase1Result& a, const Phase1Result& b) {
            if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
            if (a.sig.det != b.sig.det) return std::abs(a.sig.det) < std::abs(b.sig.det);
            return a.catalog_id < b.catalog_id;
        });

        // Write
        const ExternalSpec* rep = tag_to_spec[tag];
        std::string fname = dir + "/" + tag + ".cat";
        write_phase1_cat(fname, *rep, deduped);

        index << tag << " | " << deduped.size() << " | " << tag << ".cat\n";
        total_bases += (int)deduped.size();
    }

    index << "\n# Total: " << total_bases << " bases across "
          << all_specs.size() << " specs\n";
    index.close();

    // Non-SUGRA output
    int total_nonsugra = 0;
    if (config.save_nonsugra) {
        std::string ns_dir = dir + "_nonsugra";
        std::filesystem::create_directories(ns_dir);
        for (auto& [tag, group] : nonsugra_by_tag) {
            if (group.empty()) continue;
            // Dedup non-SUGRA
            std::set<std::string> ns_seen;
            std::vector<Phase1Result> ns_deduped;
            for (auto& r : group) {
                std::string extra = tag;
                // Include catalog_id (LST source) so different build paths are not merged
                extra += "|cid:" + std::to_string(r.catalog_id);
                std::string key = r.eigenvalues.empty()
                    ? dedup_key_v2(r.final_IF, extra)
                    : dedup_key_v2_cached(r.final_IF, r.eigenvalues, extra);
                if (ns_seen.insert(key).second) ns_deduped.push_back(std::move(r));
            }
            std::sort(ns_deduped.begin(), ns_deduped.end(),
                      [](const Phase1Result& a, const Phase1Result& b) {
                if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
                return a.catalog_id < b.catalog_id;
            });
            const ExternalSpec* rep = tag_to_spec[tag];
            std::string fname = ns_dir + "/" + tag + ".cat";
            write_phase1_cat(fname, *rep, ns_deduped);
            total_nonsugra += (int)ns_deduped.size();
        }
        std::cout << "Non-SUGRA: " << total_nonsugra << " bases to " << ns_dir << "/\n";
    }

    // Summary
    std::cout << "\n=== Phase1 Summary ===\n";
    std::cout << "Total: " << total_bases << " bases\n";
    std::cout << "Output: " << dir << "/\n";
    std::cout << "Index:  " << dir << "/INDEX.txt\n";


    // ── LaTeX ──
    {
        std::string texfile = dir + "/results.tex";
        std::ofstream tex(texfile);
        tex << "\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n"
            << "\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n\\usepackage[dvipsnames]{xcolor}\n"
            << "\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n"
            << "\\section*{1-External SUGRA Bases}\n\\textcolor{red}{Red} = external.\n\n"
            << "Total: " << total_bases << " bases.\n\n\\end{document}\n";
        tex.close();
        std::cout << "LaTeX: " << texfile << "\n";
    }

    return 0;
}
