// gen_sugra_nhc_ext_phase2.cpp — Phase 2 for NHC external catalogs
//
// Reads cat_nhc_ext/ or previous NHC ext Phase 2 output.
// Attaches single-curve externals + NHC cluster externals.
// Output: cat_nhc_ext_Next/ (N = round number)
//
// Usage: ./gen_sugra_nhc_ext_phase2 <input_dir> [round=2]

#include "sugra_generator.h"

static double g_max_ext_cc = 16.0;
static bool g_det_sq_mode = false;
static bool g_save_nonsugra = false;
static bool g_no_su2 = false;
static bool g_no_223 = false;
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <numeric>
#include <map>
#include <algorithm>
#include <functional>
#include <filesystem>
#include <cassert>
#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// Types
// ============================================================================

struct ExternalSpec {
    int id, ext_si, target_si, int_num;
    bool is_hat1;
    double cc;
    std::string label, tag;
};

inline double compute_ext_cc(int es, int in, bool h, int ts) {
    if (h) return (ts == -1) ? 7.0 : (ts == -2) ? 15.0 : 10.0;
    GaugeInfo g = gauge_from_si(es);
    if (g.dim > 0) return central_charge(g, in);
    if (es == -2 && ts == -1) return central_charge(GAUGE_SU2, in);
    if (es == -2 && ts == -3) return central_charge(GAUGE_SU2, in);  // su2: 3/3=1.0
    if (es == -2 && ts == -4) return central_charge(GAUGE_SU8, in);   // su8: 63*2/(2+8)=12.6
    if (es == -3 && ts == -2) return central_charge(GAUGE_G2, in);   // g2: 14/5=2.8
    if (es == -4 && ts == -2) return central_charge(GAUGE_SO16, in); // so16: 120*2/(2+14)=15.0
    return 0.0;
}

inline std::vector<ExternalSpec> build_all_specs() {
    std::vector<ExternalSpec> sp;
    int id = 0;
    auto a = [&](int es, int ts, int in, bool h, const std::string& l, const std::string& t) {
        sp.push_back({id++, es, ts, in, h, compute_ext_cc(es, in, h, ts), l, t});
    };
    for (int k = 1; k <= 2; k++)
        a(-2, -1, k, false, "su(2) k=" + std::to_string(k), "su2");
    a(-2, -3, 1, false, "su2->su3", "su2n3");
    a(-2, -4, 2, false, "su8->so16", "su8");
    for (int k = 1; k <= 2; k++)
        a(-3, -1, k, false, "su(3) k=" + std::to_string(k), "su3");
    a(-3, -2, 1, false, "g2->su2", "su3n2");
    a(-4, -1, 1, false, "so(8)", "so8");
    a(-4, -2, 2, false, "so8->su2", "so16n2");
    a(-5, -1, 1, false, "$\\mathfrak{f}_4$", "f4");
    a(-6, -1, 1, false, "$\\mathfrak{e}_6$", "e6");
    a(-7, -1, 1, false, "$\\mathfrak{e}_7'$", "e7p");
    a(-8, -1, 1, false, "$\\mathfrak{e}_7$", "e7");
    a(-12, -1, 1, false, "$\\mathfrak{e}_8$", "e8");
    a(-1, -1, 1, true, "$\\hat{1}{\\to}1$", "hat1m1");
    a(-1, -2, 1, true, "$\\hat{1}{\\to}2$", "hat1m2");
    return sp;
}

// NHC cluster spec (for attaching whole NHC as external)
struct NHCClusterSpec {
    int id;
    std::string tag, label;
    std::vector<int> self_ints;
    struct InternalInt { int a, b, k; };
    std::vector<InternalInt> internal_ints;
    int H, V;
};

inline std::vector<NHCClusterSpec> build_nhc_cluster_specs() {
    std::vector<NHCClusterSpec> specs;
    int id = 100;  // offset to not clash with single-curve spec ids
    specs.push_back({id++, "nhc_2_3", "(-2)-(-3) NHC",
                     {-2, -3}, {{0, 1, 1}}, 8, 17});
    specs.push_back({id++, "nhc_2_3_2", "(-2)-(-3)-(-2) NHC",
                     {-2, -3, -2}, {{0, 1, 1}, {1, 2, 1}}, 16, 27});
    // Curve order: [-3, -2, -2] (process -3 first → cc-prune early)
    specs.push_back({id++, "nhc_2_2_3", "(-2)-(-2)-(-3) NHC",
                     {-3, -2, -2}, {{0, 1, 1}, {1, 2, 1}}, 8, 17});
    specs.push_back({id++, "nhc_2_4", "(-2)-(-4) NHC",
                     {-2, -4}, {{0, 1, 2}}, 128, 183});
    return specs;
}

struct ExtInfo {
    int curve_idx, spec_id;
    std::string tag;
    int ext_si, target_si, int_num;
    bool is_hat1;
};

struct CatEntry {
    int entry_id;
    std::string combo_key;
    std::vector<ExtInfo> externals;
    int catalog_id;
    std::string catalog_type;
    int base_T;
    int T, H_charged, V, H_neutral, det, sig_pos, sig_neg, sig_zero;
    Eigen::MatrixXi final_IF;
    // Sparse in-memory storage of the (symmetric, ~3.5%-dense) base IF for LOADED input
    // entries: holds upper-triangle nonzeros as flat [i,j,v, i,j,v, ...]. final_IF is left
    // empty until build_dense() reconstructs it just-in-time when the entry is processed,
    // then free_dense() releases it. Cuts the huge chain input from ~17 GB to ~1 GB.
    // Output entries (built during attach) keep a dense final_IF and leave if_n == 0.
    int if_n = 0;
    std::vector<int> if_nz;
    void build_dense() {
        if (if_n == 0) return;  // no sparse store → leave final_IF as-is
        final_IF = Eigen::MatrixXi::Zero(if_n, if_n);
        for (size_t k = 0; k + 2 < if_nz.size(); k += 3) {
            int i = if_nz[k], j = if_nz[k + 1], v = if_nz[k + 2];
            final_IF(i, j) = v;
            if (i != j) final_IF(j, i) = v;
        }
    }
    void free_dense() { if (if_n != 0) final_IF = Eigen::MatrixXi(); }
    // Store a dense matrix sparsely. Used for OUTPUT entries so the accumulating chain
    // results stay small in RAM; dedup/write reconstruct via build_dense() per entry.
    void set_sparse_from(const Eigen::MatrixXi& M) {
        if_n = (int)M.rows();
        if_nz.clear();
        for (int i = 0; i < if_n; i++)
            for (int j = i; j < if_n; j++)
                if (M(i, j) != 0) { if_nz.push_back(i); if_nz.push_back(j); if_nz.push_back(M(i, j)); }
    }
    struct RemTarget { int curve_idx, self_int; double available_cc; };
    std::vector<RemTarget> remaining;
    // Cached for true_c_ext: LST sub-IF kernel (= curves not in any external).
    // Computed once per entry; reused by every attach function and every leaf.
    // Empty vector means kernel undefined → filter is skipped per existing semantics.
    std::vector<int> cached_lst_idx;
    mutable Eigen::VectorXi cached_kernel;  // LST sub-IF null vector (JacobiSVD),
    mutable bool kernel_ready = false;      // lazily computed (entry_kernel) in attach
    mutable BaseSchur schur;           // base inertia + pseudo-inverse (lazily computed),
    mutable bool schur_ready = false;  // skipped for chain entries that attach no candidate
};

// Populate cached_lst_idx and cached_kernel for an entry. Idempotent; safe to call
// multiple times. Called once after entry is loaded / created.
inline void populate_lst_cache(CatEntry& e) {
    // Input entries store the IF sparsely (final_IF empty until build_dense); use if_n.
    int n = e.if_n != 0 ? e.if_n : (int)e.final_IF.rows();
    std::set<int> ext_set;
    for (auto& ei : e.externals) ext_set.insert(ei.curve_idx);
    e.cached_lst_idx.clear();
    for (int i = 0; i < n; i++)
        if (ext_set.count(i) == 0) e.cached_lst_idx.push_back(i);
    // cached_kernel (JacobiSVD) and schur (eigendecomp) are deferred to entry_kernel /
    // entry_schur on the first attach candidate. This moves the per-entry heavy linear
    // algebra out of the SERIAL load loop into the PARALLEL attach region (≈Ncore speedup)
    // and skips it entirely for entries that attach nothing — the dominant cost on the
    // huge non-SUGRA chain input once parsing was made fast. Verified byte-identical &
    // deterministic: each entry is owned by one OpenMP thread, and Eigen disables nested
    // threading inside the parallel region, so concurrent on-demand decompositions are safe.
}

// Lazy per-entry LST sub-IF null vector (JacobiSVD), computed on first attach candidate.
inline const Eigen::VectorXi& entry_kernel(const CatEntry& e) {
    if (!e.kernel_ready) {
        int m = (int)e.cached_lst_idx.size();
        if (m == 0) { e.cached_kernel = Eigen::VectorXi(); }
        else {
            Eigen::MatrixXi A(m, m);
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++)
                    A(i, j) = e.final_IF(e.cached_lst_idx[i], e.cached_lst_idx[j]);
            e.cached_kernel = find_positive_integer_kernel(A);
        }
        e.kernel_ready = true;
    }
    return e.cached_kernel;
}

// Lazy per-entry base inertia + pseudo-inverse for schur_sig.
inline const BaseSchur& entry_schur(const CatEntry& e) {
    if (!e.schur_ready) { e.schur = precompute_base_schur(e.final_IF); e.schur_ready = true; }
    return e.schur;
}

// Fast signature pre-filter via Schur complement against the cached base inertia
// (O((z+k)³) vs O(n³) LDLT). Conservative threshold (1e-6) so a true sig_pos==1 is
// never rejected; compute_sig_fast downstream stays authoritative. Mirrors the
// reject condition of sig_pos_exceeds_one_fast (sig_pos>1 OR sig_neg>T_max).
inline bool schur_reject(const BaseSchur& bs, const Eigen::MatrixXi& M, int T_max) {
    SigInfo s = schur_sig(bs, M, 1e-6);
    return s.sig_pos > 1 || s.sig_neg > T_max;
}

static std::vector<CatEntry>* g_nonsugra_results = nullptr;
#ifdef _OPENMP
// Each OpenMP thread points g_nonsugra_results to its own thread-local nonsugra
// vector; main thread merges after the parallel region.
#pragma omp threadprivate(g_nonsugra_results)
#endif

// NHC cluster curve → actual gauge cc
inline double nhc_curve_cc(const std::string& tag, int ext_si, int int_num, int target_si = -1) {
    if (tag == "nhc_2_3") {
        if (ext_si == -3 && target_si == -2)
            return central_charge(GAUGE_SO7, int_num) + central_charge(GAUGE_SU2, int_num);
        if (ext_si == -2) return central_charge(GAUGE_SU2, int_num);
        if (ext_si == -3) return central_charge(GAUGE_G2, int_num);
    } else if (tag == "nhc_2_3_2") {
        if (ext_si == -2) return central_charge(GAUGE_SU2, int_num);
        if (ext_si == -3) return central_charge(GAUGE_SO7, int_num);
    } else if (tag == "nhc_2_2_3") {
        if (ext_si == -3) return central_charge(GAUGE_G2, int_num);
        if (ext_si == -2) return central_charge(GAUGE_SP1, int_num);
    } else if (tag == "nhc_2_4") {
        if (ext_si == -2) return central_charge(GAUGE_SU8, int_num);
        if (ext_si == -4) return central_charge(GAUGE_SO16, int_num);
    }
    return 0.0;
}

inline double total_ext_cc(const std::vector<ExtInfo>& exts,
                           int new_ext_si = 0, int new_target_si = 0, int new_int_num = 0, bool new_is_hat1 = false) {
    double total = 0.0;
    for (auto& e : exts) {
        if (e.tag.substr(0, 4) == "nhc_")
            total += nhc_curve_cc(e.tag, e.ext_si, e.int_num, e.target_si);
        else
            total += compute_ext_cc(e.ext_si, e.int_num, e.is_hat1, e.target_si);
    }
    if (new_ext_si != 0 || new_is_hat1)
        total += compute_ext_cc(new_ext_si, new_int_num, new_is_hat1, new_target_si);
    return total;
}

inline bool should_save_nonsugra(const AnomalyResult& anom, const std::vector<ExtInfo>& all_exts) {
    double cur_cc = total_ext_cc(all_exts);
    return !recovery_impossible(anom.H_neutral, g_max_ext_cc - cur_cc);
}

// ============================================================================
// I/O
// ============================================================================

// Fast whitespace-separated field parsing for the hot .cat load path. Per-line
// std::istringstream (formatted >>) dominated load time when parsing the large
// non-SUGRA chain input (5.5 GB / 234k entries at T=108: sample showed ~all time
// in load_cat). These parse identical values, so catalog output stays byte-identical.
inline const char* skip_ws(const char* p, const char* end) {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
    return p;
}
inline const char* fast_int(const char* p, const char* end, int& out) {
    p = skip_ws(p, end);
    bool neg = false;
    if (p < end && (*p == '+' || *p == '-')) { neg = (*p == '-'); ++p; }
    long v = 0;
    while (p < end && *p >= '0' && *p <= '9') v = v * 10 + (*p++ - '0');
    out = neg ? (int)(-v) : (int)v;
    return p;
}

inline std::vector<CatEntry> load_cat(const std::string& path) {
    std::vector<CatEntry> entries;
    std::ifstream fin(path);
    if (!fin) return entries;
    std::string line;
    while (std::getline(fin, line)) {
        if (line.substr(0, 5) != "ENTRY") continue;
        CatEntry e;
        e.entry_id = std::stoi(line.substr(6));

        std::getline(fin, line);
        if (line.substr(0, 4) == "SPEC") {
            std::istringstream ss(line.substr(5));
            ExtInfo ei; int h;
            ss >> ei.spec_id >> ei.tag >> ei.ext_si >> ei.target_si >> ei.int_num >> h;
            ei.is_hat1 = (h != 0);
            e.combo_key = ei.tag;
            std::getline(fin, line);
            { std::istringstream ss2(line.substr(5)); ss2 >> e.catalog_id >> e.catalog_type >> e.base_T; }
            std::getline(fin, line);
            { std::istringstream ss2(line.substr(7));
              int target_idx, ext_curve_idx;
              ss2 >> target_idx >> ext_curve_idx;
              ei.curve_idx = ext_curve_idx;
              e.externals.push_back(ei);
            }
            std::getline(fin, line);
            { std::istringstream ss2(line.substr(8));
              ss2 >> e.T >> e.H_charged >> e.V >> e.H_neutral >> e.det >> e.sig_pos >> e.sig_neg >> e.sig_zero; }
        } else if (line.substr(0, 5) == "COMBO") {
            e.combo_key = line.substr(6);
            std::getline(fin, line);
            int n_ext = std::stoi(line.substr(10));
            for (int i = 0; i < n_ext; i++) {
                std::getline(fin, line);
                std::istringstream ss(line);
                ExtInfo ei; int h;
                ss >> ei.curve_idx >> ei.spec_id >> ei.tag >> ei.ext_si >> ei.target_si >> ei.int_num >> h;
                ei.is_hat1 = (h != 0);
                e.externals.push_back(ei);
            }
            std::getline(fin, line);
            { std::istringstream ss(line.substr(5)); ss >> e.catalog_id >> e.catalog_type >> e.base_T; }
            std::getline(fin, line);
            { std::istringstream ss(line.substr(8));
              ss >> e.T >> e.H_charged >> e.V >> e.H_neutral >> e.det >> e.sig_pos >> e.sig_neg >> e.sig_zero; }
        } else {
            continue;
        }

        std::getline(fin, line);
        int n = std::stoi(line.substr(3));
        // Parse the dense IF text but store only upper-triangle nonzeros (sparse). final_IF
        // is reconstructed on demand by build_dense() when the entry is processed.
        e.if_n = n;
        for (int i = 0; i < n; i++) {
            std::getline(fin, line);
            const char* p = line.c_str(); const char* end = p + line.size();
            for (int j = 0; j < n; j++) {
                int v; p = fast_int(p, end, v);
                if (v != 0 && i <= j) { e.if_nz.push_back(i); e.if_nz.push_back(j); e.if_nz.push_back(v); }
            }
        }
        std::getline(fin, line);
        int nr = std::stoi(line.substr(10));
        for (int i = 0; i < nr; i++) {
            std::getline(fin, line);
            const char* p = line.c_str(); const char* end = p + line.size();
            CatEntry::RemTarget rt;
            p = fast_int(p, end, rt.curve_idx);
            p = fast_int(p, end, rt.self_int);
            rt.available_cc = std::strtod(skip_ws(p, end), nullptr);
            e.remaining.push_back(rt);
        }
        std::getline(fin, line); // END
        populate_lst_cache(e);  // pre-compute LST kernel once per entry
        entries.push_back(std::move(e));
    }
    return entries;
}

inline std::vector<CatEntry> load_all_cats(const std::string& dir) {
    std::vector<CatEntry> all;
    for (auto& p : std::filesystem::directory_iterator(dir)) {
        if (p.path().extension() != ".cat") continue;
        for (auto& x : load_cat(p.path().string())) all.push_back(std::move(x));
    }
    return all;
}

// ============================================================================
// Dedup + remaining
// ============================================================================

// Binary-fingerprint variant of dedup_key_v2.
// Produces a 128-bit key (two-hash) over the same canonical content as
// dedup_key_v2 (sorted eigenvalues + sorted vertex descriptors + sorted edge
// multiset + extra) but packs bytes directly instead of formatting decimal
// strings. Used with std::unordered_set<DedupKey128, DedupKey128Hash>.
inline DedupKey128 dedup_key_v2_binary(const Eigen::MatrixXi& IF,
                                        const std::string& extra = "") {
    int n = IF.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    auto ev = es.eigenvalues();

    std::vector<uint8_t> buf;
    buf.reserve(64 + n * 32 + extra.size());

    // 1. Sorted eigenvalues quantized to int64 at 1e-8 resolution (matches the
    //    8-decimal text format used by dedup_key_v2).
    for (int i = 0; i < n; i++) {
        double v = ev(i);
        if (std::abs(v) < 1e-10) v = 0.0;
        int64_t qv = (int64_t)std::llround(v * 1e8);
        buf_append_i64(buf, qv);
    }
    buf.push_back(0xFF);  // section separator

    // 2. Sorted vertex descriptors
    struct VDesc {
        int si, deg;
        std::vector<int> nbr_si;
        bool operator<(const VDesc& o) const {
            return std::tie(si, deg, nbr_si) < std::tie(o.si, o.deg, o.nbr_si);
        }
    };
    std::vector<VDesc> vds(n);
    for (int i = 0; i < n; i++) {
        vds[i].si = IF(i, i);
        vds[i].deg = 0;
        for (int j = 0; j < n; j++)
            if (j != i && IF(i, j) != 0) {
                vds[i].deg++;
                vds[i].nbr_si.push_back(IF(j, j));
            }
        std::sort(vds[i].nbr_si.begin(), vds[i].nbr_si.end());
    }
    std::sort(vds.begin(), vds.end());
    for (auto& vd : vds) {
        buf_append_i32(buf, vd.si);
        buf_append_i32(buf, vd.deg);
        buf_append_i32(buf, (int32_t)vd.nbr_si.size());
        for (int s : vd.nbr_si) buf_append_i32(buf, s);
    }
    buf.push_back(0xFE);

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
    for (auto& [a, b, w] : edges) {
        buf_append_i32(buf, a);
        buf_append_i32(buf, b);
        buf_append_i32(buf, w);
    }
    buf.push_back(0xFD);

    // 4. Extra discriminator
    buf_append_bytes(buf, extra.data(), extra.size());

    return hash128(buf);
}

inline std::string dedup_key_v2(const Eigen::MatrixXi& IF, const std::string& extra = "") {
    int n = IF.rows();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(IF.cast<double>());
    auto ev = es.eigenvalues();
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(8) << "E:";
    for (int i = 0; i < ev.size(); i++) {
        double v = ev(i); if (std::abs(v) < 1e-10) v = 0.0;
        if (i) ss << ","; ss << v;
    }
    struct VDesc { int si, deg; std::vector<int> nbr_si;
        bool operator<(const VDesc& o) const {
            return std::tie(si, deg, nbr_si) < std::tie(o.si, o.deg, o.nbr_si); }};
    std::vector<VDesc> vds(n);
    for (int i = 0; i < n; i++) {
        vds[i].si = IF(i,i); vds[i].deg = 0;
        for (int j = 0; j < n; j++)
            if (j != i && IF(i,j) != 0) { vds[i].deg++; vds[i].nbr_si.push_back(IF(j,j)); }
        std::sort(vds[i].nbr_si.begin(), vds[i].nbr_si.end());
    }
    std::sort(vds.begin(), vds.end());
    ss << "|V:";
    for (auto& vd : vds) {
        ss << "(" << vd.si << "," << vd.deg << ",{";
        for (int k = 0; k < (int)vd.nbr_si.size(); k++) { if (k) ss << ","; ss << vd.nbr_si[k]; }
        ss << "})";
    }
    std::vector<std::tuple<int,int,int>> edges;
    for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
            if (IF(i,j) != 0) { int a=IF(i,i),b=IF(j,j); if(a>b)std::swap(a,b); edges.push_back({a,b,IF(i,j)}); }
    std::sort(edges.begin(), edges.end());
    ss << "|E:";
    for (auto& [a,b,w] : edges) ss << "(" << a << "," << b << "," << w << ")";
    if (!extra.empty()) ss << "|X:" << extra;
    return ss.str();
}

inline std::vector<CatEntry::RemTarget> compute_remaining(
    const Eigen::MatrixXi& IF, const NHCResult& nhc,
    const std::set<int>& ext_idx, double cc = 8.0)
{
    std::vector<CatEntry::RemTarget> ts;
    int n = IF.rows();
    for (int i = 0; i < n; i++) {
        if (ext_idx.count(i)) continue;
        CatEntry::RemTarget t;
        t.curve_idx = i; t.self_int = IF(i,i); t.available_cc = -1.0;
        if (t.self_int == -1) {
            double used = 0;
            for (int j = 0; j < n; j++) {
                if (j == i || IF(i,j) == 0) continue;
                GaugeInfo g = (j < (int)nhc.curve_gauges.size() && nhc.curve_gauges[j].dim > 0)
                    ? nhc.curve_gauges[j] : gauge_from_si(IF(j,j));
                if (IF(j,j) == -1 && g.dim == 0) continue;
                int k = std::abs(IF(i,j));
                if (g.dim > 0) used += central_charge(g, k);
            }
            t.available_cc = cc - used;
        }
        if (t.self_int == -1 && t.available_cc > 0.1) ts.push_back(t);
        else if (t.self_int <= -2) ts.push_back(t);
    }
    return ts;
}

// ============================================================================
// Write
// ============================================================================

inline void write_cat_entry(std::ofstream& out, int eid, const CatEntry& r) {
    int n = r.final_IF.rows();
    out << "ENTRY " << eid << "\n"
        << "COMBO " << r.combo_key << "\n"
        << "EXTERNALS " << r.externals.size() << "\n";
    for (auto& ei : r.externals)
        out << ei.curve_idx << " " << ei.spec_id << " " << ei.tag << " "
            << ei.ext_si << " " << ei.target_si << " " << ei.int_num
            << " " << (ei.is_hat1 ? 1 : 0) << "\n";
    out << "BASE " << r.catalog_id << " " << r.catalog_type << " " << r.base_T << " " << n << "\n"
        << "PHYSICS " << r.T << " " << r.H_charged << " " << r.V << " " << r.H_neutral
        << " " << r.det << " " << r.sig_pos << " " << r.sig_neg << " " << r.sig_zero << "\n"
        << "IF " << n << "\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) { if (j) out << " "; out << r.final_IF(i,j); }
        out << "\n";
    }
    out << "REMAINING " << r.remaining.size() << "\n";
    for (auto& t : r.remaining)
        out << t.curve_idx << " " << t.self_int << " "
            << std::fixed << std::setprecision(2) << t.available_cc << "\n";
    out << "END\n\n";
}

// ============================================================================
// Helpers
// ============================================================================

inline bool enhance_all_externals(
    NHCResult& nhc, const Eigen::MatrixXi& new_IF,
    const std::vector<ExtInfo>& all_exts, double cc_budget)
{
    for (auto& ei : all_exts) {
        if (ei.is_hat1)
            enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
        else if (ei.ext_si == -2 && ei.target_si == -1)
            enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
        if (!nhc.passes) return false;
    }
    return true;
}

inline CatEntry make_result(
    const CatEntry& parent, const std::string& combo_tag,
    const std::vector<ExtInfo>& new_exts,
    const Eigen::MatrixXi& new_IF,
    const NHCResult& nhc, const AnomalyResult& anom,
    const SigInfo& sig, double cc_budget)
{
    CatEntry r;
    r.externals = parent.externals;
    for (auto& ne : new_exts) r.externals.push_back(ne);
    // Build sorted combo_key for dedup (order-independent). An NHC cluster attaches
    // one ExtInfo PER CURVE (needed downstream), but in the combo it must count as
    // ONE cluster, not one tag per curve: a -2-2-3 cluster is "nhc_2_2_3", NOT
    // "nhc_2_2_3+nhc_2_2_3+nhc_2_2_3". Collapse each NHC tag's curve count to its
    // cluster count (curves / curves-per-cluster). This matches the seed (gen_sugra
    // _nhc_ext writes COMBO = single cluster tag) so seed and chain are consistent.
    {
        auto nhc_curves_per = [](const std::string& t) -> int {
            if (t == "nhc_2_3")   return 2;
            if (t == "nhc_2_3_2") return 3;
            if (t == "nhc_2_2_3") return 3;
            if (t == "nhc_2_4")   return 2;
            return 1;   // single-curve external
        };
        std::map<std::string,int> cnt;
        for (auto& ei : r.externals) cnt[ei.tag]++;
        std::vector<std::string> tags;
        for (auto& kv : cnt) {
            int per = nhc_curves_per(kv.first);
            int n = kv.second / per;   // cluster instances (curve count is a multiple of `per`)
            for (int i = 0; i < n; i++) tags.push_back(kv.first);
        }
        std::sort(tags.begin(), tags.end());
        r.combo_key = "";
        for (size_t i = 0; i < tags.size(); i++) {
            if (i) r.combo_key += "+";
            r.combo_key += tags[i];
        }
    }
    r.catalog_id = parent.catalog_id;
    r.catalog_type = parent.catalog_type;
    r.base_T = parent.base_T;
    r.set_sparse_from(new_IF);   // store sparsely; reconstructed at dedup/write
    r.T = anom.T; r.H_charged = anom.H_charged; r.V = anom.V; r.H_neutral = anom.H_neutral;
    r.det = sig.det; r.sig_pos = sig.sig_pos; r.sig_neg = sig.sig_neg; r.sig_zero = sig.sig_zero;
    std::set<int> ei_set;
    for (auto& ei : r.externals) ei_set.insert(ei.curve_idx);
    r.remaining = compute_remaining(new_IF, nhc, ei_set, cc_budget);
    return r;
}

inline bool has_hat1(const CatEntry& entry, bool new_is_hat1 = false) {
    if (new_is_hat1) return true;
    for (auto& ei : entry.externals) if (ei.is_hat1) return true;
    return false;
}

// Canonical-order skip: non-hat1 externals must be tag-non-decreasing relative
// to the entry's last non-hat1 block. hat1 attach always allowed.
// Output-preserving: dedup_key_v2 (graph-iso) absorbs duplicate paths to same
// final IF; canonical just picks one representative.
// --no-canonical: disable the tag-non-decreasing skip so EVERY attachment order
// is explored (dedup_key_v2 still removes duplicate paths). Required when the
// seed is restricted to a late-canonical tag (e.g. a 223-only seed cannot attach
// any tag < "nhc_2_2_3" under the canonical rule, losing all such combos).
inline bool g_no_canonical = false;

inline bool canonical_skip(const CatEntry& entry,
                            const std::string& new_tag, bool new_is_hat1) {
    if (g_no_canonical) return false;
    if (new_is_hat1) return false;
    for (auto it = entry.externals.rbegin(); it != entry.externals.rend(); ++it) {
        if (!it->is_hat1) return new_tag < it->tag;
    }
    return false;
}

// ============================================================================
// Mode A: Single-curve external attachment
// ============================================================================

inline void attach_single(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    if (canonical_skip(entry, spec.tag, spec.is_hat1)) return;
    std::set<int> ext_set;
    for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);

    for (auto& rt : entry.remaining) {
        if (rt.self_int != spec.target_si) continue;
        // hat1→(-2): target must be standalone (not part of NHC chain)
        if (spec.is_hat1 && spec.target_si == -2) {
            bool standalone = true;
            int n = entry.final_IF.rows();
            for (int j = 0; j < n; j++) {
                if (j == rt.curve_idx || entry.final_IF(rt.curve_idx, j) == 0) continue;
                if (ext_set.count(j)) continue;
                if (entry.final_IF(j, j) != -1) { standalone = false; break; }
            }
            if (!standalone) continue;
        }
        if (spec.target_si == -1) {
            GaugeInfo g = gauge_from_si(spec.ext_si);
            double need = 0;
            if (g.dim > 0) need = central_charge(g, spec.int_num);
            else if (spec.ext_si == -2) need = central_charge(GAUGE_SU2, spec.int_num);
            if (spec.is_hat1 && spec.target_si == -1) need = central_charge(GAUGE_SU8, spec.int_num);
            if (need > rt.available_cc + 1e-9) continue;
        }
        // c_ext check
        if (total_ext_cc(entry.externals, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1) > g_max_ext_cc + 1e-9) continue;

        tried++;
        int n = entry.final_IF.rows();
        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
        int ext_idx = n;
        new_IF(ext_idx, ext_idx) = spec.ext_si;
        new_IF(ext_idx, rt.curve_idx) = spec.int_num;
        new_IF(rt.curve_idx, ext_idx) = spec.int_num;

        if (!spec.is_hat1 && schur_reject(entry_schur(entry), new_IF, 193)) continue;
        new_IF.block(0,0,n,n) = entry.final_IF;  // base copied only after schur passes (schur reads only appended cols)
        auto sig = compute_sig_fast(new_IF);
        if (!spec.is_hat1) {
            if (sig.sig_pos != 1) continue;
        }
        if (g_det_sq_mode && !has_hat1(entry, spec.is_hat1) && !is_perfect_square(std::abs(sig.det))) continue;

        NHCResult nhc = check_nhc(new_IF, cc_budget);
        if (!nhc.passes) continue;

        // Build full ext list for enhancement
        std::vector<ExtInfo> all_exts = entry.externals;
        ExtInfo new_ei = {ext_idx, spec.id, spec.tag, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1};
        all_exts.push_back(new_ei);
        if (!enhance_all_externals(nhc, new_IF, all_exts, cc_budget)) continue;

        // Proper c_ext ≤ 16 prune (NHC-aware bridge-level)
        std::set<int> ext_set;
        for (auto& ei : all_exts) ext_set.insert(ei.curve_idx);
        if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
            continue;
        // True c_ext ≤ 16 — REJECT (no recovery)
        {
            double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
            if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) continue;
        }

        AnomalyResult anom = compute_anomaly(nhc, sig);
        bool is_sugra = true;
        if (anom.H_neutral < 0) is_sugra = false;
        if (is_sugra) {
            std::set<int> h1;
            for (auto& ei : entry.externals) if (ei.is_hat1) h1.insert(ei.curve_idx);
            if (spec.is_hat1) h1.insert(ext_idx);
            if (reject_nonuni(new_IF, sig, anom, h1)) is_sugra = false;
        }
        auto res = make_result(entry, spec.tag, {new_ei}, new_IF, nhc, anom, sig, cc_budget);
        if (is_sugra) { results.push_back(std::move(res)); passed++; }
        else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, res.externals)) g_nonsugra_results->push_back(std::move(res));
    }
}

// ============================================================================
// Mode B: v6 multi-target (standard ext → multiple -1 curves)
// ============================================================================

inline void attach_v6_multi(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    if (spec.is_hat1 || spec.target_si != -1) return;
    if (canonical_skip(entry, spec.tag, false)) return;
    // c_ext pre-check (at least 1 ext added)
    if (total_ext_cc(entry.externals, spec.ext_si, spec.target_si, spec.int_num, false) > g_max_ext_cc + 1e-9) return;

    GaugeInfo ext_gauge = gauge_from_si(spec.ext_si);
    if (ext_gauge.dim == 0 && spec.ext_si == -2) ext_gauge = GAUGE_SU2;

    std::vector<int> m1_curves;
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
    if ((int)m1_curves.size() < 2) return;

    int max_k, max_tgt;
    if (spec.ext_si == -2 || spec.ext_si == -3) {
        max_k = 5;  // match phase1 config.mixed_int_max default
        max_tgt = (int)m1_curves.size();
    } else {
        max_k = std::min(max_level_for_cc(ext_gauge, 8.0), 9);
        int max_level = std::min(max_level_for_cc(ext_gauge, 20.0), 9);
        max_tgt = std::min({(int)m1_curves.size(), max_level});
    }

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    // Per-(-1)-target available_cc lookup (built once, used in pre-check below).
    // Mirrors attach_single's line 545 check across all targets.
    int n_total = entry.final_IF.rows();
    std::vector<double> avail_cc(n_total, -1.0);
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1) avail_cc[rt.curve_idx] = rt.available_cc;

    enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](const std::vector<int>& targets) {
        enumerate_int_nums((int)targets.size(), max_k, [&](const std::vector<int>& int_nums) {
            tried++;
            // Pre-check: each (-1) target must have room for the new ext's cc
            // contribution. Single attach (-1 target) does not change existing
            // non-(-1) cluster gauges, so available_cc upper-bounds the cc
            // added on each target. Skips most L2 fails before new_IF / NHC.
            bool fits = true;
            for (int i = 0; i < (int)targets.size(); i++) {
                double need = central_charge(ext_gauge, int_nums[i]);
                double avail = avail_cc[targets[i]];
                if (avail < 0 || need > avail + 1e-9) { fits = false; break; }
            }
            if (!fits) return;

            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
            int ext_idx = n;
            new_IF(ext_idx, ext_idx) = spec.ext_si;
            for (int i = 0; i < (int)targets.size(); i++) {
                new_IF(ext_idx, targets[i]) = int_nums[i];
                new_IF(targets[i], ext_idx) = int_nums[i];
            }


            if (schur_reject(entry_schur(entry), new_IF, 193)) return;
            new_IF.block(0,0,n,n) = entry.final_IF;  // base copied only after schur passes (schur reads only appended cols)
            auto sig = compute_sig_fast(new_IF);
            if (sig.sig_pos != 1) return;
            if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) return;

            std::vector<ExtInfo> all_exts = entry.externals;
            ExtInfo new_ei = {ext_idx, spec.id, spec.tag, spec.ext_si, -1, int_nums[0], false};
            all_exts.push_back(new_ei);
            if (!enhance_all_externals(nhc, new_IF, all_exts, cc_budget)) return;

            // Proper c_ext ≤ 16 prune (NHC-aware bridge-level)
            std::set<int> ext_set;
            for (auto& ei : all_exts) ext_set.insert(ei.curve_idx);
            if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
                return;
            // True c_ext ≤ 16 — REJECT
            {
                double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
                if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
            }

            AnomalyResult anom = compute_anomaly(nhc, sig);
            bool v6_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto v6_res = make_result(entry, spec.tag, {new_ei}, new_IF, nhc, anom, sig, cc_budget);
            if (v6_sugra) { results.push_back(std::move(v6_res)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, v6_res.externals)) g_nonsugra_results->push_back(std::move(v6_res));
        });
    });
}

// ============================================================================
// Mode C: NHC cluster attachment
// ============================================================================

inline void attach_nhc_cluster(
    const CatEntry& entry, const NHCClusterSpec& nhc_spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    // No canonical_skip: NHC cluster attach on a non-NHC input entry (e.g.,
    // [su3] + attach nhc_2_3) produces final IFs that the canonical path
    // (input nhc_2_3 + attach su3 via single/v6_multi) cannot reach due to
    // intermediate cc-budget divergence. Cluster attach enumerates both
    // directions; dedup absorbs duplicate canonical paths.
    // c_ext pre-check: estimate minimum cc for NHC cluster (using actual NHC gauge)
    {
        double nhc_cc = 0;
        for (int si : nhc_spec.self_ints) nhc_cc += nhc_curve_cc(nhc_spec.tag, si, 1);
        if (total_ext_cc(entry.externals) + nhc_cc > g_max_ext_cc + 1e-9) return;
    }

    // Collect (-1) curves with sufficient cc budget + (-2) curves (only for nhc_2_3)
    std::vector<int> m1_curves;
    std::vector<int> all_targets;
    for (auto& rt : entry.remaining) {
        if (rt.self_int == -1 && rt.available_cc >= 1.0 - 1e-9) {
            m1_curves.push_back(rt.curve_idx);
            all_targets.push_back(rt.curve_idx);
        } else if (rt.self_int == -2 && nhc_spec.tag == "nhc_2_3") {
            all_targets.push_back(rt.curve_idx);
        }
    }
    if (m1_curves.empty()) return;  // still need at least one (-1)

    int ne = (int)nhc_spec.self_ints.size();

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    // Enumerate assignments: each NHC curve → one (-1) base curve
    // Pre-compute available cc per (-1) remaining curve
    std::map<int, double> avail_cc;
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1) avail_cc[rt.curve_idx] = rt.available_cc;

    // Track incremental cc usage from NHC assignments
    std::map<int, double> cc_used;

    std::function<void(int, std::vector<std::vector<std::pair<int,int>>>&)> enumerate;
    enumerate = [&](int nhc_idx, std::vector<std::vector<std::pair<int,int>>>& assignment) {
        if (nhc_idx == ne) {
            tried++;
            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n + ne, n + ne);

            // NHC internal structure
            for (int i = 0; i < ne; i++)
                new_IF(n + i, n + i) = nhc_spec.self_ints[i];
            for (auto& ii : nhc_spec.internal_ints) {
                new_IF(n + ii.a, n + ii.b) = ii.k;
                new_IF(n + ii.b, n + ii.a) = ii.k;
            }
            // NHC → base connections (multi-target)
            for (int i = 0; i < ne; i++) {
                for (auto& [tgt, inum] : assignment[i]) {
                    new_IF(n + i, tgt) += inum;
                    new_IF(tgt, n + i) += inum;
                }
            }

            // Schur signature first: it reads only the appended cluster columns, NOT the
            // base block, so we skip the O(n²) base copy + quick_cc_check for the ~93% it
            // rejects. Reordering vs quick_cc_check is output-preserving (both only reject).
            if (schur_reject(entry_schur(entry), new_IF, 193)) return;
            new_IF.block(0, 0, n, n) = entry.final_IF;  // base copied only after schur passes
            if (!quick_cc_check(new_IF, cc_budget)) return;
            auto sig = compute_sig_fast(new_IF);
            if (sig.sig_pos != 1) return;
            if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) return;

            // Enhance existing externals
            std::vector<ExtInfo> all_exts = entry.externals;
            std::vector<ExtInfo> new_exts;
            for (int i = 0; i < ne; i++) {
                ExtInfo ei;
                ei.curve_idx = n + i;
                ei.spec_id = nhc_spec.id;
                ei.tag = nhc_spec.tag;
                ei.ext_si = nhc_spec.self_ints[i];
                ei.target_si = entry.final_IF(assignment[i][0].first, assignment[i][0].first);
                ei.int_num = assignment[i][0].second;
                ei.is_hat1 = false;
                all_exts.push_back(ei);
                new_exts.push_back(ei);
            }
            if (!enhance_all_externals(nhc, new_IF, all_exts, cc_budget)) return;

            // Proper c_ext ≤ 16 prune (NHC-aware bridge-level)
            std::set<int> ext_set;
            for (auto& ei : all_exts) ext_set.insert(ei.curve_idx);
            if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
                return;
            // True c_ext ≤ 16 — REJECT
            {
                double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
                if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
            }

            AnomalyResult anom = compute_anomaly(nhc, sig);
            bool nc_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto nc_res = make_result(entry, nhc_spec.tag, new_exts, new_IF, nhc, anom, sig, cc_budget);
            if (nc_sugra) { results.push_back(std::move(nc_res)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, nc_res.externals)) g_nonsugra_results->push_back(std::move(nc_res));
            return;
        }

        // CC-aware assignment: check budget before recursing
        int nhc_si = nhc_spec.self_ints[nhc_idx];
        // NHC-aware gauge for cc accounting. A (-2) inside an NHC chain carries SU2
        // (NOT none), so it consumes cc on its (-1) target. gauge_from_si(-2)=NONE -> add=0
        // -> incremental prune never fires -> multi-target enumeration explodes. Charging the
        // true NHC gauge makes the prune fire and matches check_nhc (output-preserving).
        GaugeInfo nhc_g;
        if (nhc_spec.tag == "nhc_2_3_2") {            // su2 x so7 x su2
            static const GaugeInfo* g[] = {&GAUGE_SU2, &GAUGE_SO7, &GAUGE_SU2};
            nhc_g = (nhc_idx < 3) ? *g[nhc_idx] : gauge_from_si(nhc_si);
        } else if (nhc_spec.tag == "nhc_2_2_3") {     // g2 x su2 x none (spec order -3,-2,-2)
            static const GaugeInfo* g[] = {&GAUGE_G2, &GAUGE_SU2, &GAUGE_NONE};
            nhc_g = (nhc_idx < 3) ? *g[nhc_idx] : gauge_from_si(nhc_si);
        } else if (nhc_spec.tag == "nhc_2_4") {       // su8 x so16
            static const GaugeInfo* g[] = {&GAUGE_SU8, &GAUGE_SO16};
            nhc_g = (nhc_idx < 2) ? *g[nhc_idx] : gauge_from_si(nhc_si);
        } else if (nhc_spec.tag == "nhc_2_3") {       // su2 x g2
            nhc_g = (nhc_si == -2) ? GAUGE_SU2 : (nhc_si == -3 ? GAUGE_G2 : gauge_from_si(nhc_si));
        } else {
            nhc_g = gauge_from_si(nhc_si);
        }

        auto try_assignment = [&](const std::vector<std::pair<int,int>>& tgts) {
            for (auto& [t, k] : tgts) {
                if (avail_cc.find(t) == avail_cc.end()) continue; // non-(-1) target
                double add = (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
                if (cc_used[t] + add > avail_cc[t] + 1e-9) return;
            }
            for (auto& [t, k] : tgts) {
                if (avail_cc.find(t) == avail_cc.end()) continue;
                cc_used[t] += (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
            }
            assignment.push_back(tgts);
            enumerate(nhc_idx + 1, assignment);
            assignment.pop_back();
            for (auto& [t, k] : tgts) {
                if (avail_cc.find(t) == avail_cc.end()) continue;
                cc_used[t] -= (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
            }
        };

        // Per-curve minimum cc (NHC-aware gauge at k=1) for target pre-filtering
        double curve_min_cc = (nhc_g.dim > 0) ? central_charge(nhc_g, 1) : 0.0;
        if (nhc_spec.tag == "nhc_2_3_2") {
            static const double tbl[] = {1.0, 3.5, 1.0};
            if (nhc_idx < 3) curve_min_cc = tbl[nhc_idx];
        } else if (nhc_spec.tag == "nhc_2_2_3") {
            static const double tbl[] = {2.8, 1.0, 0.0};
            if (nhc_idx < 3) curve_min_cc = tbl[nhc_idx];
        } else if (nhc_spec.tag == "nhc_2_4") {
            static const double tbl[] = {6.3, 15.0};
            if (nhc_idx < 2) curve_min_cc = tbl[nhc_idx];
        } else if (nhc_spec.tag == "nhc_2_3") {
            int si = nhc_spec.self_ints[nhc_idx];
            if (si == -2) curve_min_cc = 1.0;
            else if (si == -3) curve_min_cc = 2.8;
        }

        // Per-curve target list filtered by avail_cc ≥ curve_min_cc
        std::vector<int> curve_targets, curve_m1_targets;
        for (int t : all_targets) {
            if (entry.final_IF(t, t) == -1) {
                auto it = avail_cc.find(t);
                if (it != avail_cc.end() && it->second >= curve_min_cc - 1e-9) {
                    curve_targets.push_back(t);
                    curve_m1_targets.push_back(t);
                }
            } else {
                curve_targets.push_back(t);
            }
        }

        // nhc_2_2_3: restrict to single-target with k=1 only
        bool restrict_223 = (nhc_spec.tag == "nhc_2_2_3");

        // Single-target options
        for (int t : curve_targets) {
            int kmax = restrict_223 ? 1 : (entry.final_IF(t,t) == -1 ? 2 : 1);
            for (int k = 1; k <= kmax; k++) {
                try_assignment({{t, k}});
            }
        }
        // Multi-target from per-curve filtered list (skipped for nhc_2_2_3)
        int max_multi = restrict_223 ? 0 : std::min((int)curve_m1_targets.size(), 3);
        if (max_multi >= 2) {
            enumerate_subsets_v2(curve_m1_targets, 2, max_multi, [&](const std::vector<int>& targets) {
                enumerate_int_nums((int)targets.size(), 2, [&](const std::vector<int>& knums) {
                    std::vector<std::pair<int,int>> mt;
                    for (size_t i = 0; i < targets.size(); i++)
                        mt.push_back({targets[i], knums[i]});
                    try_assignment(mt);
                });
            });
        }
    };

    std::vector<std::vector<std::pair<int,int>>> assignment;
    enumerate(0, assignment);
}

// ============================================================================
// Mode D: Mixed multi-target
//   su2n3mix:  su2(-2) → (-3) int=1 + (-1) int=k
//   su8mix:    su2(-2) → (-4) int=2 + (-1) int=k
//   su3n2mix:  su3(-3) → (-2) int=1 + (-1) int=k
//   so16n2mix: so8(-4) → (-2) int=2 + (-1) int=k
// ============================================================================

// Generic helper: attach ext_curve (ext_si) to one special target (tgt_curve, tgt_int)
// + multiple (-1) curves simultaneously
inline void attach_mixed_generic(
    const CatEntry& entry, double cc_budget,
    int ext_si, int tgt_curve, int tgt_int,
    const std::vector<int>& m1_curves,
    const std::string& mix_tag, int spec_id,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    // Canonical only for mix_tags that have a phase1 base (alternate path
    // reachable). su2n3mix / su8mix / so16n2mix are phase2-emergent (no phase1
    // cat), so applying canonical here loses final IFs.
    if (mix_tag == "su3n2mix") {
        if (canonical_skip(entry, mix_tag, false)) return;
    }
    // c_ext pre-check
    if (total_ext_cc(entry.externals, ext_si, -1, 1, false) > g_max_ext_cc + 1e-9) return;

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    int int_max = 5;  // match phase1 config.mixed_int_max
    int max_m1 = std::min((int)m1_curves.size(), 193);  // match phase1 (no effective cap)

    // Per-(-1)-target available_cc lookup + ext_gauge for pre-check.
    int n_total = entry.final_IF.rows();
    std::vector<double> avail_cc(n_total, -1.0);
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1) avail_cc[rt.curve_idx] = rt.available_cc;
    GaugeInfo ext_gauge_for_check = gauge_from_si(ext_si);
    if (ext_gauge_for_check.dim == 0 && ext_si == -2) ext_gauge_for_check = GAUGE_SU2;

    enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
        enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
            tried++;
            // Pre-check: each (-1) target must have room for ext_gauge cc.
            // Conservative (uses standalone gauge); cluster joining via
            // tgt_curve may make actual gauge larger (still safe — only
            // under-prunes).
            bool fits = true;
            for (int i = 0; i < (int)m1_targets.size(); i++) {
                double need = central_charge(ext_gauge_for_check, int_nums[i]);
                double avail = avail_cc[m1_targets[i]];
                if (avail < 0 || need > avail + 1e-9) { fits = false; break; }
            }
            if (!fits) return;

            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
            int ext_idx = n;
            new_IF(ext_idx, ext_idx) = ext_si;
            // Connect to special target
            new_IF(ext_idx, tgt_curve) = tgt_int;
            new_IF(tgt_curve, ext_idx) = tgt_int;
            // Connect to (-1)s
            for (int i = 0; i < (int)m1_targets.size(); i++) {
                new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                new_IF(m1_targets[i], ext_idx) = int_nums[i];
            }

            if (schur_reject(entry_schur(entry), new_IF, 193)) return;
            new_IF.block(0,0,n,n) = entry.final_IF;  // base copied only after schur passes
            auto sig = compute_sig_fast(new_IF);
            if (sig.sig_pos != 1) return;
            if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) return;

            std::vector<ExtInfo> all_exts = entry.externals;
            int tgt_si = entry.final_IF(tgt_curve, tgt_curve);
            ExtInfo new_ei = {ext_idx, spec_id, mix_tag, ext_si, tgt_si, tgt_int, false};
            all_exts.push_back(new_ei);
            if (!enhance_all_externals(nhc, new_IF, all_exts, cc_budget)) return;

            // True c_ext ≤ 16 — REJECT
            {
                std::set<int> ext_set;
                for (auto& ei : all_exts) ext_set.insert(ei.curve_idx);
                double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
                if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
            }
            AnomalyResult anom = compute_anomaly(nhc, sig);
            bool mx_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto mx_res = make_result(entry, mix_tag, {new_ei}, new_IF, nhc, anom, sig, cc_budget);
            if (mx_sugra) { results.push_back(std::move(mx_res)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, mx_res.externals)) g_nonsugra_results->push_back(std::move(mx_res));
        });
    });
}

inline void attach_mixed_multi(
    const CatEntry& entry, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    std::vector<int> m1_curves, m2_curves, m3_curves, m4_curves;
    for (auto& rt : entry.remaining) {
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -2) m2_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -3) m3_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -4) m4_curves.push_back(rt.curve_idx);
    }
    // so7 = (-2)-(-3)-(-2): the new (-3) bridges two (-2) curves into su2+so7+su2.
    // A (-2) can be an so7 leg ONLY if it is "isolated" — all its base-IF neighbors are
    // (-1). A (-2) already adjacent to a deeper curve (≤ -2, fiber, etc.) would make a
    // different/invalid cluster, so any pair containing it produces nothing downstream.
    // Pre-filtering to isolated (-2)s avoids the C(m2,2) blow-up: at high T almost every
    // (-2) sits in a chain, so there are ~0 valid pairs (so7/so7mix become near free).
    std::vector<int> iso_m2_curves;
    {
        int nn = entry.final_IF.rows();
        for (int c : m2_curves) {
            bool iso = true;
            for (int j = 0; j < nn; j++) {
                if (j == c || entry.final_IF(c, j) == 0) continue;
                if (entry.final_IF(j, j) != -1) { iso = false; break; }
            }
            if (iso) iso_m2_curves.push_back(c);
        }
    }

    // Don't early-return on m1_curves.empty() — the so7 inline block below only
    // uses m2_curves and must still run. Gate m1-dependent blocks individually.

    if (!m1_curves.empty()) {
        // su2n3mix: su2(-2) → (-3) int=1 + (-1)s
        for (int m3 : m3_curves)
            attach_mixed_generic(entry, cc_budget, -2, m3, 1, m1_curves,
                                 "su2n3mix", 7/*su2n3 spec_id*/, results, tried, passed);

        // su8mix: su2(-2) → (-4) int=2 + (-1)s
        for (int m4 : m4_curves)
            attach_mixed_generic(entry, cc_budget, -2, m4, 2, m1_curves,
                                 "su8mix", 8/*su8 spec_id*/, results, tried, passed);

        // su3n2mix: su3(-3) → (-2) int=1 + (-1)s
        for (int m2 : m2_curves)
            attach_mixed_generic(entry, cc_budget, -3, m2, 1, m1_curves,
                                 "su3n2mix", 16/*su3n2 spec_id*/, results, tried, passed);

        // so16n2mix: so8(-4) → (-2) int=2 + (-1)s
        for (int m2 : m2_curves)
            attach_mixed_generic(entry, cc_budget, -4, m2, 2, m1_curves,
                                 "so16n2mix", 18/*so16n2 spec_id*/, results, tried, passed);
    }

    // so7: (-3) → two (-2) curves (forms (-2)-(-3)-(-2) = su2+so7+su2)
    if ((int)iso_m2_curves.size() >= 2 && !canonical_skip(entry, "so7", false)) {
        std::set<int> hat1_set;
        for (auto& ei : entry.externals)
            if (ei.is_hat1) hat1_set.insert(ei.curve_idx);
        for (int a = 0; a < (int)iso_m2_curves.size(); a++) {
            for (int b = a + 1; b < (int)iso_m2_curves.size(); b++) {
                tried++;
                int n = entry.final_IF.rows();
                Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                int ext_idx = n;
                new_IF(ext_idx, ext_idx) = -3;
                new_IF(ext_idx, iso_m2_curves[a]) = 1; new_IF(iso_m2_curves[a], ext_idx) = 1;
                new_IF(ext_idx, iso_m2_curves[b]) = 1; new_IF(iso_m2_curves[b], ext_idx) = 1;

                if (schur_reject(entry_schur(entry), new_IF, 193)) continue;
                new_IF.block(0,0,n,n) = entry.final_IF;  // base copied only after schur passes
                auto sig = compute_sig_fast(new_IF);
                if (sig.sig_pos != 1) continue;

                NHCResult nhc = check_nhc(new_IF, cc_budget);
                if (!nhc.passes) continue;
                if (!enhance_all_externals(nhc, new_IF, entry.externals, cc_budget)) continue;

                // True c_ext ≤ 16 — REJECT
                {
                    std::set<int> ext_set;
                    for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);
                    ext_set.insert(ext_idx);
                    double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
                    if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) continue;
                }
                AnomalyResult anom = compute_anomaly(nhc, sig);
                bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                ExtInfo new_ei = {ext_idx, 0, "so7", -3, -2, 1, false};
                auto res = make_result(entry, "so7", {new_ei}, new_IF, nhc, anom, sig, cc_budget);
                if (is_sugra) { results.push_back(std::move(res)); passed++; }
                else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, res.externals)) g_nonsugra_results->push_back(std::move(res));
            }
        }
    }

    // so7mix: (-3) → two (-2) + (-1)s
    if ((int)iso_m2_curves.size() >= 2 && !m1_curves.empty() && !canonical_skip(entry, "so7mix", false)) {
        std::set<int> hat1_set;
        for (auto& ei : entry.externals)
            if (ei.is_hat1) hat1_set.insert(ei.curve_idx);
        // Per-(-1) available_cc lookup for the fits pre-check below.
        std::vector<double> avail_cc_m(entry.final_IF.rows(), -1.0);
        for (auto& rt : entry.remaining)
            if (rt.self_int == -1) avail_cc_m[rt.curve_idx] = rt.available_cc;
        for (int a = 0; a < (int)iso_m2_curves.size(); a++) {
            for (int b = a + 1; b < (int)iso_m2_curves.size(); b++) {
                int max_m1 = std::min((int)m1_curves.size(), 193);  // match phase1
                enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                    enumerate_int_nums((int)m1_targets.size(), 5, [&](const std::vector<int>& int_nums) {
                        tried++;
                        // Per-curve cc pre-check (mirrors attach_mixed_generic's `fits`):
                        // each (-1) target must have budget for the external's cc load. Uses
                        // the standalone su3 gauge; the cluster gauge (so7) is ≥ su3, so this
                        // conservatively UNDER-prunes → output-preserving (overloaded targets
                        // fail downstream NHC/cc anyway). Skips ~99% of leaves before new_IF.
                        {
                          bool fits = true;
                          for (int i = 0; i < (int)m1_targets.size(); i++) {
                              double need = central_charge(GAUGE_SU3, int_nums[i]);
                              double av = avail_cc_m[m1_targets[i]];
                              if (av < 0 || need > av + 1e-9) { fits = false; break; }
                          }
                          if (!fits) return;
                        }
                        int n = entry.final_IF.rows();
                        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                        int ext_idx = n;
                        new_IF(ext_idx, ext_idx) = -3;
                        new_IF(ext_idx, iso_m2_curves[a]) = 1; new_IF(iso_m2_curves[a], ext_idx) = 1;
                        new_IF(ext_idx, iso_m2_curves[b]) = 1; new_IF(iso_m2_curves[b], ext_idx) = 1;
                        for (int i = 0; i < (int)m1_targets.size(); i++) {
                            new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                            new_IF(m1_targets[i], ext_idx) = int_nums[i];
                        }

                        if (schur_reject(entry_schur(entry), new_IF, 193)) return;
                        new_IF.block(0,0,n,n) = entry.final_IF;  // base copied only after schur passes
                        auto sig = compute_sig_fast(new_IF);
                        if (sig.sig_pos != 1) return;

                        NHCResult nhc = check_nhc(new_IF, cc_budget);
                        if (!nhc.passes) return;
                        if (!enhance_all_externals(nhc, new_IF, entry.externals, cc_budget)) return;

                        // True c_ext ≤ 16 — REJECT
                        {
                            std::set<int> ext_set;
                            for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);
                            ext_set.insert(ext_idx);
                            double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set, entry.cached_lst_idx, entry_kernel(entry));
                            if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
                        }
                        AnomalyResult anom = compute_anomaly(nhc, sig);
                        bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                        ExtInfo new_ei = {ext_idx, 0, "so7mix", -3, -2, 1, false};
                        auto res = make_result(entry, "so7mix", {new_ei}, new_IF, nhc, anom, sig, cc_budget);
                        if (is_sugra) { results.push_back(std::move(res)); passed++; }
                        else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, res.externals)) g_nonsugra_results->push_back(std::move(res));
                    });
                });
            }
        }
    }

    // 2 / 2mix: new (-2, gauge NONE) attaches to a LST-native (-2) that is
    // currently part of an isolated nhc_2_3 cluster (the partner (-3) must
    // also be LST-native and have no other non-(-1) neighbors). Forms
    // nhc_2_2_3 after attach. cc/anomaly impact: zero (new (-2) has no
    // gauge; existing (-2)'s SU2→SP1 has identical cc; H,V same between
    // nhc_2_3 and nhc_2_2_3).
    {
        // Identify externals and valid target (-2)s
        std::set<int> ext_set_2mix;
        for (auto& ei : entry.externals) ext_set_2mix.insert(ei.curve_idx);

        int n = entry.final_IF.rows();
        std::vector<int> targets_2;
        for (int i = 0; i < n; i++) {
            if (ext_set_2mix.count(i)) continue;  // i must be LST native
            if (entry.final_IF(i, i) != -2) continue;
            // i's non-(-1) neighbors must be exactly one LST-native (-3)
            int j_partner = -1;
            int i_non_m1_count = 0;
            for (int j = 0; j < n; j++) {
                if (j == i || entry.final_IF(i, j) == 0) continue;
                if (entry.final_IF(j, j) == -1) continue;
                i_non_m1_count++;
                j_partner = j;
            }
            if (i_non_m1_count != 1) continue;
            if (ext_set_2mix.count(j_partner)) continue;  // partner must be LST native
            if (entry.final_IF(j_partner, j_partner) != -3) continue;
            // j_partner's non-(-1) neighbors must be exactly {i}
            int j_non_m1_count = 0;
            for (int k = 0; k < n; k++) {
                if (k == j_partner || entry.final_IF(j_partner, k) == 0) continue;
                if (entry.final_IF(k, k) == -1) continue;
                j_non_m1_count++;
            }
            if (j_non_m1_count != 1) continue;
            targets_2.push_back(i);
        }

        if (!targets_2.empty()) {
            std::set<int> hat1_set;
            for (auto& ei : entry.externals)
                if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

            auto run_2_attach = [&](int target_2,
                                    const std::vector<int>& m1_targets,
                                    const std::vector<int>& int_nums,
                                    const std::string& tag) {
                tried++;
                Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n + 1, n + 1);
                int ext_idx = n;
                new_IF(ext_idx, ext_idx) = -2;
                new_IF(ext_idx, target_2) = 1;
                new_IF(target_2, ext_idx) = 1;
                for (int i = 0; i < (int)m1_targets.size(); i++) {
                    new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                    new_IF(m1_targets[i], ext_idx) = int_nums[i];
                }

                if (schur_reject(entry_schur(entry), new_IF, 193)) return;
                new_IF.block(0, 0, n, n) = entry.final_IF;  // base copied only after schur passes
                auto sig = compute_sig_fast(new_IF);
                if (sig.sig_pos != 1) return;
                if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

                NHCResult nhc = check_nhc(new_IF, cc_budget);
                if (!nhc.passes) return;

                std::vector<ExtInfo> all_exts = entry.externals;
                ExtInfo new_ei = {ext_idx, 0, tag, -2, -2, 1, false};
                all_exts.push_back(new_ei);
                if (!enhance_all_externals(nhc, new_IF, all_exts, cc_budget)) return;

                std::set<int> ext_set_new;
                for (auto& ei : all_exts) ext_set_new.insert(ei.curve_idx);
                if (compute_proper_c_ext(new_IF, nhc, ext_set_new) > g_proper_c_ext_max + 1e-9) return;
                {
                    double tc = compute_true_c_ext_fiber_cached(new_IF, nhc, ext_set_new,
                                                                 entry.cached_lst_idx, entry_kernel(entry));
                    if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
                }

                AnomalyResult anom = compute_anomaly(nhc, sig);
                bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                auto res = make_result(entry, tag, {new_ei}, new_IF, nhc, anom, sig, cc_budget);
                if (is_sugra) { results.push_back(std::move(res)); passed++; }
                else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, res.externals))
                    g_nonsugra_results->push_back(std::move(res));
            };

            // Pure "2": new (-2) connects only to one target (-2)
            if (!canonical_skip(entry, "2", false)) {
                for (int target_2 : targets_2) {
                    run_2_attach(target_2, {}, {}, "2");
                }
            }

            // "2mix": + 1 or 2 (-1) targets with int_num in 1..g_twomix_int_max (default 3)
            if (!g_no_2mix && !m1_curves.empty() && !canonical_skip(entry, "2mix", false)) {
                for (int target_2 : targets_2) {
                    int max_m1 = std::min((int)m1_curves.size(), 2);  // user's spec: at most 2 (-1) targets
                    enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                        enumerate_int_nums((int)m1_targets.size(), g_twomix_int_max, [&](const std::vector<int>& int_nums) {
                            run_2_attach(target_2, m1_targets, int_nums, "2mix");
                        });
                    });
                }
            }
        }
    }
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_dir> [round=2] [--det-sq]\n";
        return 1;
    }
    std::string input_dir = argv[1];
    std::string round_str = "2";
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--det-sq") g_det_sq_mode = true;
        else if (std::string(argv[i]) == "--use-lst-T") { /* now canonical, no-op */ }
        else if (std::string(argv[i]) == "--save-nonsugra") g_save_nonsugra = true;
        else if (std::string(argv[i]) == "--no-su2") g_no_su2 = true;
        else if (std::string(argv[i]) == "--no-223") g_no_223 = true;
        else if (std::string(argv[i]) == "--2mix-intmax" && i + 1 < argc) g_twomix_int_max = std::atoi(argv[++i]);
        else if (std::string(argv[i]) == "--no-2mix") g_no_2mix = true;
        else if (std::string(argv[i]) == "--no-canonical") g_no_canonical = true;
        else if (std::string(argv[i]).rfind("--", 0) == 0) { /* unknown flag: ignore, don't eat round arg */ }
        else round_str = argv[i];
    }

    std::cout << "=== gen_sugra_nhc_ext_phase2 ===\n"
              << "Input: " << input_dir << "\n"
              << "Round: " << round_str
              << (g_det_sq_mode ? "  [|det|=n² filter ON]" : "") << "\n\n";

    auto all_specs = build_all_specs();
    auto nhc_cluster_specs = build_nhc_cluster_specs();
    double cc_budget = 8.0;

    auto input_entries = load_all_cats(input_dir);
    std::cout << "Loaded: " << input_entries.size() << " entries\n\n";

    std::string out_dir = "cat_nhc_ext_" + round_str + "ext";
    std::filesystem::remove_all(out_dir);
    std::filesystem::create_directories(out_dir);

    std::vector<CatEntry> results;
    int tried_single = 0, passed_single = 0;
    int tried_v6 = 0, passed_v6 = 0;
    int tried_nhc = 0, passed_nhc = 0;
    int tried_mix = 0, passed_mix = 0;

    std::vector<CatEntry> nonsugra;

    const int n_in = (int)input_entries.size();
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    // Per-thread result vectors; merged after the parallel region.
    std::vector<std::vector<CatEntry>> tl_results(nthreads);
    std::vector<std::vector<CatEntry>> tl_nonsugra(nthreads);

#ifdef _OPENMP
    omp_set_dynamic(0);  // pin thread count so threadprivate stays consistent
#endif
#pragma omp parallel num_threads(nthreads) \
    reduction(+:tried_single,passed_single,tried_v6,passed_v6,tried_nhc,passed_nhc,tried_mix,passed_mix)
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        // Each thread routes nonsugra pushes into its own slot.
        if (g_save_nonsugra) g_nonsugra_results = &tl_nonsugra[tid];

#pragma omp for schedule(static)
        for (int ei = 0; ei < n_in; ei++) {
            auto& entry = input_entries[ei];
            entry.build_dense();   // reconstruct dense IF from sparse store (once per entry)
            g_current_base_T = (int)entry.final_IF.rows() - (int)entry.externals.size();
            auto& my_res = tl_results[tid];

            // Single-curve externals (no ordering constraint — try all specs)
            for (auto& spec : all_specs) {
                if (g_no_su2 && spec.ext_si == -2 && spec.target_si == -1 && !spec.is_hat1) continue;
                attach_single(entry, spec, cc_budget, my_res, tried_single, passed_single);
                attach_v6_multi(entry, spec, cc_budget, my_res, tried_v6, passed_v6);
            }

            // Mixed multi-target (all variants: su2n3mix, su8mix, su3n2mix, so16n2mix, so7, so7mix, 2, 2mix)
            attach_mixed_multi(entry, cc_budget, my_res, tried_mix, passed_mix);

            // NHC cluster externals
            for (auto& ns : nhc_cluster_specs) {
                if (g_no_223 && ns.tag == "nhc_2_2_3") continue;
                attach_nhc_cluster(entry, ns, cc_budget, my_res, tried_nhc, passed_nhc);
            }
            entry.free_dense();    // release the dense IF; sparse store retained
        }
    }

    // Merge thread-local vectors into the main results / nonsugra vectors.
    {
        size_t total_res = 0, total_ns = 0;
        for (auto& v : tl_results) total_res += v.size();
        for (auto& v : tl_nonsugra) total_ns += v.size();
        results.reserve(total_res);
        for (auto& v : tl_results) for (auto& r : v) results.push_back(std::move(r));
        tl_results.clear();
        if (g_save_nonsugra) {
            nonsugra.reserve(total_ns);
            for (auto& v : tl_nonsugra) for (auto& r : v) nonsugra.push_back(std::move(r));
            tl_nonsugra.clear();
            // Restore main pointer for any post-loop usage.
            g_nonsugra_results = &nonsugra;
        }
    }

    std::cout << "\nSingle-target:    tried=" << tried_single << ", passed=" << passed_single << "\n"
              << "v6 multi-target:  tried=" << tried_v6 << ", passed=" << passed_v6 << "\n"
              << "NHC cluster:      tried=" << tried_nhc << ", passed=" << passed_nhc << "\n"
              << "Mixed multi:      tried=" << tried_mix << ", passed=" << passed_mix << "\n"
              << "Total raw: " << results.size() << "\n";

    // Dedup
    std::unordered_set<DedupKey128, DedupKey128Hash> seen;
    seen.reserve(results.size());
    std::vector<CatEntry> filtered;
    for (auto& r : results) {
        // Include catalog_id (LST source) so different build paths → different keys
        r.build_dense();
        DedupKey128 key = dedup_key_v2_binary(r.final_IF, r.combo_key + "|cid:" + std::to_string(r.catalog_id));
        r.free_dense();
        if (seen.insert(key).second) filtered.push_back(std::move(r));
    }
    std::cout << "After dedup: " << filtered.size() << "\n";

    // Sort
    std::sort(filtered.begin(), filtered.end(), [](const CatEntry& a, const CatEntry& b) {
        if (a.combo_key != b.combo_key) return a.combo_key < b.combo_key;
        if (a.T != b.T) return a.T < b.T;
        return std::abs(a.det) < std::abs(b.det);
    });

    // Write .cat files
    std::map<std::string, std::vector<CatEntry*>> by_combo;
    for (auto& r : filtered) by_combo[r.combo_key].push_back(&r);
    int tw = 0;
    for (auto& [combo, entries] : by_combo) {
        std::ofstream out(out_dir + "/" + combo + ".cat");
        out << "# NHC ext Phase2 catalog: " << combo << "\n# Round: " << round_str
            << "\n# Total: " << entries.size() << "\n#\n";
        int eid = 0;
        for (auto* rp : entries) { rp->build_dense(); write_cat_entry(out, eid++, *rp); rp->free_dense(); tw++; }
    }

    // INDEX
    { std::ofstream idx(out_dir + "/INDEX.txt");
      idx << "# NHC ext Phase2 | Round: " << round_str << " | Input: " << input_dir
          << " | Total: " << tw << " bases, " << by_combo.size() << " combos\n#\n";
      for (auto& [c, e] : by_combo) idx << c << " | " << e.size() << "\n"; }

    int n_uni = 0, n_nonuni = 0;
    for (auto& r : filtered) { if (std::abs(r.det) == 1) n_uni++; else n_nonuni++; }
    std::cout << "\nWritten " << tw << " bases (" << by_combo.size() << " combos) to "
              << out_dir << "/\nUnimodular: " << n_uni << ", Non-unimodular: " << n_nonuni << "\n";

    // LaTeX
    { std::string texfile = out_dir + "/results.tex";
      std::ofstream tex(texfile);
      tex << "\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n"
          << "\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n"
          << "\\usepackage[dvipsnames]{xcolor}\n\\setlength{\\parindent}{0pt}\n"
          << "\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n";
      tex << "\\section*{NHC ext Phase2 Round " << round_str << "}\n\n"
          << "Total: " << filtered.size() << " (" << n_uni << " uni, " << n_nonuni << " non-uni).\n\n";
      std::map<int,int> Tc; for (auto& r : filtered) Tc[r.T]++;
      tex << "\\begin{center}\\begin{tabular}{c|c}\\toprule $T$ & Count \\\\\\midrule\n";
      for (auto& [T,c] : Tc) tex << T << " & " << c << " \\\\\n";
      tex << "\\bottomrule\\end{tabular}\\end{center}\n\n";
      tex << "\\end{document}\n";
      std::cout << "Written LaTeX to " << texfile << "\n";
    }

    // Write non-SUGRA catalog
    if (g_save_nonsugra && !nonsugra.empty()) {
        std::unordered_set<DedupKey128, DedupKey128Hash> ns_seen;
        ns_seen.reserve(nonsugra.size());
        std::vector<CatEntry> ns_filtered;
        for (auto& r : nonsugra) {
            // Include catalog_id (LST source) so different build paths → different keys
            r.build_dense();
            DedupKey128 key = dedup_key_v2_binary(r.final_IF, r.combo_key + "|cid:" + std::to_string(r.catalog_id));
            r.free_dense();
            if (ns_seen.insert(key).second) ns_filtered.push_back(std::move(r));
        }
        std::string ns_dir = out_dir + "_nonsugra";
        std::filesystem::remove_all(ns_dir);
        std::filesystem::create_directories(ns_dir);
        std::map<std::string, std::vector<CatEntry*>> ns_by_combo;
        for (auto& r : ns_filtered) ns_by_combo[r.combo_key].push_back(&r);
        int ns_tw = 0;
        for (auto& [combo, entries] : ns_by_combo) {
            std::ofstream out2(ns_dir + "/" + combo + ".cat");
            out2 << "# Non-SUGRA: " << combo << "\n#\n";
            int eid = 0;
            for (auto* rp : entries) { rp->build_dense(); write_cat_entry(out2, eid++, *rp); rp->free_dense(); ns_tw++; }
        }
        std::cout << "Non-SUGRA: " << ns_tw << " bases (" << ns_by_combo.size()
                  << " combos) to " << ns_dir << "/\n";
    }

    return 0;
}
