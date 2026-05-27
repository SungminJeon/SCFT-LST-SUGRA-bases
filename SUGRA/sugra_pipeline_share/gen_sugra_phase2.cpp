// gen_sugra_phase2.cpp — iteratable external stacking
//
// Reads cat_Next/ (Phase 1 or previous Phase 2 output)
// Writes cat_(N+1)ext/
//
// Three attachment modes:
//   A. Standard single-target (existing logic)
//   B. v6 multi-target: standard external → multiple (-1) curves
//   C. hat1 multi-target: hat1 → multiple (-1) curves
//
// Multi-target uses blowdown to extract level, cc(gauge, level) ≤ 20
//
// Usage: ./gen_sugra_phase2 <input_cat_dir> [output_suffix] [--det-sq]

#include "sugra_generator.h"

static bool g_det_sq_mode = false;
static double g_max_ext_cc = 16.0;
static bool g_no_su2 = false;
static bool g_save_nonsugra = false;
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

// ============================================================================
// PART 1: Types
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
    if (es == -2 && ts == -3) return central_charge(GAUGE_SU2, in);   // su2: 1.0
    if (es == -2 && ts == -4) return central_charge(GAUGE_SU8, in);   // su8: 12.6
    if (es == -3 && ts == -2) return central_charge(GAUGE_G2, in);    // g2: 2.8
    if (es == -4 && ts == -2) return central_charge(GAUGE_SO16, in);  // so16: 15.0
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
    struct RemTarget { int curve_idx, self_int; double available_cc; };
    std::vector<RemTarget> remaining;
};

static std::vector<CatEntry>* g_nonsugra_results = nullptr;

// ============================================================================
// PART 2: .cat I/O
// ============================================================================

// Compute naive total external cc from all externals + one new external
// NHC cluster curve → actual gauge cc (based on tag + ext_si)
inline double nhc_curve_cc(const std::string& tag, int ext_si, int int_num, int target_si = -1) {
    if (tag == "nhc_2_3") {
        // (-3) on (-2) target → forms (-2-3-2) → so7 gauge
        if (ext_si == -3 && target_si == -2)
            return central_charge(GAUGE_SO7, int_num) + central_charge(GAUGE_SU2, int_num);  // so7 + su2 from target (-2)
        if (ext_si == -2) return central_charge(GAUGE_SU2, int_num);
        if (ext_si == -3) return central_charge(GAUGE_G2, int_num);
    } else if (tag == "nhc_2_3_2") {
        // (-2)=su2, (-3)=so7, (-2)=su2
        if (ext_si == -2) return central_charge(GAUGE_SU2, int_num);
        if (ext_si == -3) return central_charge(GAUGE_SO7, int_num);
    } else if (tag == "nhc_2_2_3") {
        // (-2)=none, (-2)=sp1, (-3)=g2 → approximate: first -2 has no gauge
        if (ext_si == -3) return central_charge(GAUGE_G2, int_num);
        if (ext_si == -2) return central_charge(GAUGE_SP1, int_num);  // sp1 = su2
    } else if (tag == "nhc_2_4") {
        if (ext_si == -2) return central_charge(GAUGE_SU8, int_num);
        if (ext_si == -4) return central_charge(GAUGE_SO16, int_num);
    }
    return 0.0;
}

inline double total_ext_cc(const std::vector<ExtInfo>& exts) {
    double total = 0.0;
    for (auto& e : exts) {
        if (e.tag.substr(0, 4) == "nhc_")
            total += nhc_curve_cc(e.tag, e.ext_si, e.int_num, e.target_si);
        else
            total += compute_ext_cc(e.ext_si, e.int_num, e.is_hat1, e.target_si);
    }
    return total;
}

inline double total_ext_cc(const std::vector<ExtInfo>& exts,
                           int new_ext_si, int new_target_si, int new_int_num, bool new_is_hat1) {
    return total_ext_cc(exts) + compute_ext_cc(new_ext_si, new_int_num, new_is_hat1, new_target_si);
}

// Check if non-SUGRA entry can possibly recover; if not, skip saving
inline bool should_save_nonsugra(const AnomalyResult& anom, const std::vector<ExtInfo>& all_exts) {
    double cur_cc = total_ext_cc(all_exts);
    return !recovery_impossible(anom.H_neutral, g_max_ext_cc - cur_cc);
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

        // Peek next line to detect format: Phase1 (SPEC) vs Phase2 (COMBO)
        std::getline(fin, line);
        if (line.substr(0, 4) == "SPEC") {
            // Phase 1 format: SPEC spec_id tag ext_si target_si int_num hat1
            std::istringstream ss(line.substr(5));
            ExtInfo ei; int h;
            ss >> ei.spec_id >> ei.tag >> ei.ext_si >> ei.target_si >> ei.int_num >> h;
            ei.is_hat1 = (h != 0);
            e.combo_key = ei.tag;

            // BASE catalog_id type base_T n_base
            std::getline(fin, line);
            { std::istringstream ss2(line.substr(5)); ss2 >> e.catalog_id >> e.catalog_type >> e.base_T; }

            // ATTACH target_idx ext_curve_idx [multi:...]
            std::getline(fin, line);
            {
                std::istringstream ss2(line.substr(7));
                int target_idx, ext_curve_idx;
                ss2 >> target_idx >> ext_curve_idx;
                ei.curve_idx = ext_curve_idx;
                e.externals.push_back(ei);
            }

            // PHYSICS
            std::getline(fin, line);
            { std::istringstream ss2(line.substr(8));
              ss2 >> e.T >> e.H_charged >> e.V >> e.H_neutral >> e.det >> e.sig_pos >> e.sig_neg >> e.sig_zero; }
        } else if (line.substr(0, 5) == "COMBO") {
            // Phase 2 format
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
            continue; // unknown format, skip
        }

        // IF
        std::getline(fin, line);
        int n = std::stoi(line.substr(3));
        e.final_IF = Eigen::MatrixXi::Zero(n, n);
        for (int i = 0; i < n; i++) {
            std::getline(fin, line);
            std::istringstream ss(line);
            for (int j = 0; j < n; j++) ss >> e.final_IF(i, j);
        }
        // REMAINING
        std::getline(fin, line);
        int nr = std::stoi(line.substr(10));
        for (int i = 0; i < nr; i++) {
            std::getline(fin, line);
            std::istringstream ss(line);
            CatEntry::RemTarget rt;
            ss >> rt.curve_idx >> rt.self_int >> rt.available_cc;
            e.remaining.push_back(rt);
        }
        std::getline(fin, line); // END
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
// PART 3: Dedup + remaining computation
// ============================================================================

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
                // Skip gaugeless (-1) curves, but include (-1) with gauge (e.g. hat1)
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
// PART 4: .cat write
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
// PART 5: Shared attachment helpers
// ============================================================================

inline bool enhance_all_externals(
    NHCResult& nhc, const Eigen::MatrixXi& new_IF,
    const std::vector<ExtInfo>& old_exts, int new_ext_idx,
    int new_ext_si, int new_target_si, bool new_is_hat1,
    double cc_budget)
{
    for (auto& ei : old_exts) {
        if (ei.is_hat1)
            enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
        else if (ei.ext_si == -2 && ei.target_si == -1)
            enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
        if (!nhc.passes) return false;
    }
    if (new_is_hat1)
        enhance_hat1_gauge(nhc, new_IF, new_ext_idx, new_target_si, cc_budget);
    else if (new_ext_si == -2 && new_target_si == -1)
        enhance_external_m2_gauge(nhc, new_IF, new_ext_idx, new_target_si, cc_budget);
    return nhc.passes;
}

inline CatEntry make_result(
    const CatEntry& parent, const std::string& spec_tag, int spec_id,
    int ext_idx, int ext_si, int target_si, int int_num, bool is_hat1,
    const Eigen::MatrixXi& new_IF,
    const NHCResult& nhc, const AnomalyResult& anom,
    const SigInfo& sig, double cc_budget)
{
    CatEntry r;
    r.externals = parent.externals;
    r.externals.push_back({ext_idx, spec_id, spec_tag, ext_si, target_si, int_num, is_hat1});
    // Build sorted combo_key for dedup (order-independent)
    {
        std::vector<std::string> tags;
        for (auto& ei : r.externals) tags.push_back(ei.tag);
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
    r.final_IF = new_IF;
    r.T = anom.T; r.H_charged = anom.H_charged; r.V = anom.V; r.H_neutral = anom.H_neutral;
    r.det = sig.det; r.sig_pos = sig.sig_pos; r.sig_neg = sig.sig_neg; r.sig_zero = sig.sig_zero;
    std::set<int> ei_set;
    for (auto& ei : r.externals) ei_set.insert(ei.curve_idx);
    r.remaining = compute_remaining(new_IF, nhc, ei_set, cc_budget);
    return r;
}

// Helper: compute T_min quickly. Returns > 193 if hopeless.
inline int compute_tmin(const Eigen::MatrixXi& IF, const std::set<int>& hat1_set) {
    int n = IF.rows();
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n+1, n+1);
    M.block(0, 0, n, n) = IF.cast<double>();
    for (int i = 0; i < n; i++) {
        double b = hat1_set.count(i) ? -1.0 : (IF(i,i) + 2.0);
        M(i, n) = b; M(n, i) = b;
    }
    M(n, n) = 0; double B = M.determinant();
    M(n, n) = 1; double A = M.determinant() - B;
    if (std::abs(A) < 1e-10) return 9999;
    double T_crit = 9.0 - (-B / A);
    return (int)std::ceil(T_crit - 1e-9);
}

// Helper: check if entry or new spec involves hat1
inline bool has_hat1(const CatEntry& entry, bool new_is_hat1 = false) {
    if (new_is_hat1) return true;
    for (auto& ei : entry.externals) if (ei.is_hat1) return true;
    return false;
}

// ============================================================================
// PART 6: Mode A — Standard single-target
// ============================================================================

inline void attach_single_target(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    // External indices for this entry
    std::set<int> ext_set;
    for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);

    for (auto& rt : entry.remaining) {
        if (rt.self_int != spec.target_si) continue;
        // hat1→(-2): target must be a standalone (-2), not part of NHC chain
        if (spec.is_hat1 && spec.target_si == -2) {
            bool standalone = true;
            int n = entry.final_IF.rows();
            for (int j = 0; j < n; j++) {
                if (j == rt.curve_idx || entry.final_IF(rt.curve_idx, j) == 0) continue;
                if (ext_set.count(j)) continue;  // ignore external neighbors
                if (entry.final_IF(j, j) != -1) { standalone = false; break; }  // connected to non-(-1) base
            }
            if (!standalone) continue;
        }
        if (spec.target_si == -1) {
            GaugeInfo g = gauge_from_si(spec.ext_si);
            double need = 0;
            if (g.dim > 0) need = central_charge(g, spec.int_num);
            else if (spec.ext_si == -2) need = central_charge(GAUGE_SU2, spec.int_num);
            // hat1 → (-1): su(8) gauge
            if (spec.is_hat1 && spec.target_si == -1) need = central_charge(GAUGE_SU8, spec.int_num);
            if (need > rt.available_cc + 1e-9) continue;
        }
        // Total external cc check (early, before IF construction)
        if (total_ext_cc(entry.externals, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1) > g_max_ext_cc + 1e-9) continue;
        tried++;
        int n = entry.final_IF.rows();
        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
        new_IF.block(0,0,n,n) = entry.final_IF;
        int ext_idx = n;
        new_IF(ext_idx, ext_idx) = spec.ext_si;
        new_IF(ext_idx, rt.curve_idx) = spec.int_num;
        new_IF(rt.curve_idx, ext_idx) = spec.int_num;

        if (!spec.is_hat1 && sig_pos_exceeds_one_fast(new_IF, 193)) continue;
        auto sig = compute_sig(new_IF);
        if (!spec.is_hat1) {
            if (sig.sig_pos != 1) continue;
        }

        // hat1: unimodular cannot reach Hirzebruch → reject
        if (spec.is_hat1 && std::abs(sig.det) == 1) continue;
        // |det|=n² filter (skip if hat1 involved)
        if (g_det_sq_mode && !has_hat1(entry, spec.is_hat1) && !is_perfect_square(std::abs(sig.det))) continue;

        NHCResult nhc = check_nhc(new_IF, cc_budget);
        if (!nhc.passes) continue;
        if (!enhance_all_externals(nhc, new_IF, entry.externals,
                ext_idx, spec.ext_si, spec.target_si, spec.is_hat1, cc_budget))
            continue;

        // Proper c_ext ≤ 16 prune (NHC-aware bridge-level)
        std::set<int> ext_set;
        for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);
        ext_set.insert(ext_idx);
        if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
            continue;

        // True (fiber-weighted) c_ext ≤ 16 — REJECT (no recovery)
        {
            double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
            if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) continue;
        }

        AnomalyResult anom = compute_anomaly(new_IF, nhc);

        // Anomaly + reject_nonuni: fail → save to non-sugra if enabled
        bool is_sugra = true;
        if (anom.H_neutral < 0) is_sugra = false;
        if (is_sugra) {
            std::set<int> h1;
            for (auto& ei : entry.externals) if (ei.is_hat1) h1.insert(ei.curve_idx);
            if (spec.is_hat1) h1.insert(ext_idx);
            if (reject_nonuni(new_IF, sig, anom, h1)) is_sugra = false;
        }
        if (!is_sugra) {
            if (g_save_nonsugra && g_nonsugra_results) {
                double cur_cc = total_ext_cc(entry.externals, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1);
                if (!recovery_impossible(anom.H_neutral, g_max_ext_cc - cur_cc))
                    g_nonsugra_results->push_back(make_result(entry, spec.tag, spec.id,
                        ext_idx, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1,
                        new_IF, nhc, anom, sig, cc_budget));
            }
            continue;
        }

        results.push_back(make_result(entry, spec.tag, spec.id,
            ext_idx, spec.ext_si, spec.target_si, spec.int_num, spec.is_hat1,
            new_IF, nhc, anom, sig, cc_budget));
        passed++;
    }
}

// ============================================================================
// PART 7: Mode B — v6 multi-target (standard ext → multiple -1 curves)
// ============================================================================

inline void attach_v6_multi_target(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    if (spec.is_hat1 || spec.target_si != -1) return;

    GaugeInfo ext_gauge = gauge_from_si(spec.ext_si);
    if (ext_gauge.dim == 0 && spec.ext_si == -2) ext_gauge = GAUGE_SU2;

    std::vector<int> m1_curves;
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
    if ((int)m1_curves.size() < 2) return;

    // Physically motivated limits per gauge algebra (matching Phase 1)
    int max_k, max_tgt;
    if (spec.ext_si == -2) {
        max_k = 4;
        max_tgt = (int)m1_curves.size();
    } else if (spec.ext_si == -3) {
        max_k = 4;
        max_tgt = (int)m1_curves.size();
    } else {
        // per-target: cc budget=8 (each -1 curve limit)
        max_k = std::min(max_level_for_cc(ext_gauge, 8.0), 9);
        // total targets: bounded by total level from cc budget=20
        int max_level = std::min(max_level_for_cc(ext_gauge, 20.0), 9);
        max_tgt = std::min({(int)m1_curves.size(), max_level});
    }

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](const std::vector<int>& targets) {
        enumerate_int_nums((int)targets.size(), max_k, [&](const std::vector<int>& int_nums) {
            tried++;
            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
            new_IF.block(0,0,n,n) = entry.final_IF;
            int ext_idx = n;
            new_IF(ext_idx, ext_idx) = spec.ext_si;
            for (int i = 0; i < (int)targets.size(); i++) {
                new_IF(ext_idx, targets[i]) = int_nums[i];
                new_IF(targets[i], ext_idx) = int_nums[i];
            }

            // LDLT → eigensolve (cheap filters first)
            if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
            auto sig = compute_sig(new_IF);
            if (sig.sig_pos != 1) return;
            if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) return;
            if (!enhance_all_externals(nhc, new_IF, entry.externals,
                    ext_idx, spec.ext_si, -1, false, cc_budget))
                return;

            AnomalyResult anom = compute_anomaly(nhc, sig);
            if (total_ext_cc(entry.externals, spec.ext_si, -1, int_nums[0], false) > g_max_ext_cc + 1e-9) return;

            bool v6_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto v6_result = make_result(entry, spec.tag, spec.id,
                ext_idx, spec.ext_si, -1, int_nums[0], false,
                new_IF, nhc, anom, sig, cc_budget);
            if (v6_sugra) { results.push_back(std::move(v6_result)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, v6_result.externals)) g_nonsugra_results->push_back(std::move(v6_result));
        });
    });
}

// ============================================================================
// PART 7b: Mode D — Mixed multi-target: su2 → (-3) + (-1) simultaneously
// ============================================================================

inline void attach_mixed_multi_target(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    // Only for su2 (ext_si=-2) targeting -1
    if (spec.ext_si != -2 || spec.target_si != -1 || spec.is_hat1) return;

    // Collect available (-1) and (-3) curves from remaining
    std::vector<int> m1_curves, m3_curves;
    for (auto& rt : entry.remaining) {
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -3) m3_curves.push_back(rt.curve_idx);
    }
    if (m1_curves.empty() || m3_curves.empty()) return;

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    int int_max = 4;  // max int_num for (-1) targets

    for (int m3 : m3_curves) {
        int max_m1 = std::min((int)m1_curves.size(), 4);
        enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
            enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
                // Total ext cc check
                if (total_ext_cc(entry.externals, -2, -3, 1, false) > g_max_ext_cc + 1e-9) return;
                tried++;
                int n = entry.final_IF.rows();
                Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                new_IF.block(0,0,n,n) = entry.final_IF;
                int ext_idx = n;
                new_IF(ext_idx, ext_idx) = -2;  // su2 external
                // Connect to (-3) with int=1
                new_IF(ext_idx, m3) = 1; new_IF(m3, ext_idx) = 1;
                // Connect to (-1)s with varying int
                for (int i = 0; i < (int)m1_targets.size(); i++) {
                    new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                    new_IF(m1_targets[i], ext_idx) = int_nums[i];
                }

                if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
                auto sig = compute_sig(new_IF);
                if (sig.sig_pos != 1) return;

                NHCResult nhc = check_nhc(new_IF, cc_budget);
                if (!nhc.passes) return;
                for (auto& ei : entry.externals) {
                    if (ei.is_hat1) {
                        enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                    } else if (ei.ext_si == -2 && ei.target_si == -1) {
                        enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                    }
                }
                if (!nhc.passes) return;
                if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

                AnomalyResult anom = compute_anomaly(nhc, sig);
                bool mx1_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                auto mx1_result = make_result(entry, "su2n3mix", spec.id,
                    ext_idx, -2, -3, 1, false, new_IF, nhc, anom, sig, cc_budget);
                if (mx1_sugra) { results.push_back(std::move(mx1_result)); passed++; }
                else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, mx1_result.externals)) g_nonsugra_results->push_back(std::move(mx1_result));
            });
        });
    }

    // --- Also: su2 → (-4) with int=2 + (-1) simultaneously ---
    std::vector<int> m4_curves;
    for (auto& rt : entry.remaining)
        if (rt.self_int == -4) m4_curves.push_back(rt.curve_idx);

    if (!m1_curves.empty() && !m4_curves.empty()) {
        for (int m4 : m4_curves) {
            int max_m1 = std::min((int)m1_curves.size(), 4);
            enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
                    if (total_ext_cc(entry.externals, -2, -4, 2, false) > g_max_ext_cc + 1e-9) return;
                    tried++;
                    int n = entry.final_IF.rows();
                    Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                    new_IF.block(0,0,n,n) = entry.final_IF;
                    int ext_idx = n;
                    new_IF(ext_idx, ext_idx) = -2;
                    // Connect to (-4) with int=2
                    new_IF(ext_idx, m4) = 2; new_IF(m4, ext_idx) = 2;
                    // Connect to (-1)s
                    for (int i = 0; i < (int)m1_targets.size(); i++) {
                        new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                        new_IF(m1_targets[i], ext_idx) = int_nums[i];
                    }


                    NHCResult nhc = check_nhc(new_IF, cc_budget);
                    if (!nhc.passes) return;
                    for (auto& ei : entry.externals) {
                        if (ei.is_hat1)
                            enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                        else if (ei.ext_si == -2 && ei.target_si == -1)
                            enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                    }
                    if (!nhc.passes) return;

                    if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
                    auto sig = compute_sig(new_IF);
                    if (sig.sig_pos != 1) return;
                    if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

                    AnomalyResult anom = compute_anomaly(nhc, sig);
                    bool mx2_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                    auto mx2_result = make_result(entry, "su8mix", spec.id,
                        ext_idx, -2, -4, 2, false, new_IF, nhc, anom, sig, cc_budget);
                    if (mx2_sugra) { results.push_back(std::move(mx2_result)); passed++; }
                    else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, mx2_result.externals)) g_nonsugra_results->push_back(std::move(mx2_result));
                });
            });
        }
    }
}

// ============================================================================
// PART 7c: Generic mixed multi-target helper
// ext_si curve → special_target (non -1) + multiple (-1) curves
// ============================================================================

inline void attach_generic_mixed(
    const CatEntry& entry, double cc_budget,
    int ext_si, int special_target, int special_int,
    const std::vector<int>& m1_curves,
    const std::string& mix_tag,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    int int_max = 4;
    int max_m1 = std::min((int)m1_curves.size(), 4);
    enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
        enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
            tried++;
            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
            new_IF.block(0,0,n,n) = entry.final_IF;
            int ext_idx = n;
            new_IF(ext_idx, ext_idx) = ext_si;
            new_IF(ext_idx, special_target) = special_int;
            new_IF(special_target, ext_idx) = special_int;
            for (int i = 0; i < (int)m1_targets.size(); i++) {
                new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                new_IF(m1_targets[i], ext_idx) = int_nums[i];
            }

            if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
            auto sig = compute_sig(new_IF);
            if (sig.sig_pos != 1) return;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) return;
            for (auto& ei : entry.externals) {
                if (ei.is_hat1)
                    enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                else if (ei.ext_si == -2 && ei.target_si == -1)
                    enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
            }
            if (!nhc.passes) return;
            if (g_det_sq_mode && !has_hat1(entry) && !is_perfect_square(std::abs(sig.det))) return;

            AnomalyResult anom = compute_anomaly(nhc, sig);
            bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto result = make_result(entry, mix_tag, 0,
                ext_idx, ext_si, entry.final_IF(special_target, special_target), special_int, false,
                new_IF, nhc, anom, sig, cc_budget);
            if (is_sugra) { results.push_back(std::move(result)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, result.externals)) g_nonsugra_results->push_back(std::move(result));
        });
    });
}

// so7: (-3) → two (-2) curves simultaneously (no (-1) connection)
// Forms (-2)-(-3)-(-2) = su2+so7+su2
inline void attach_so7(
    const CatEntry& entry, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    std::vector<int> m1_curves, m2_curves;
    for (auto& rt : entry.remaining) {
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -2) m2_curves.push_back(rt.curve_idx);
    }
    if ((int)m2_curves.size() < 2) return;

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    // Enumerate pairs of (-2) curves
    for (int a = 0; a < (int)m2_curves.size(); a++) {
        for (int b = a + 1; b < (int)m2_curves.size(); b++) {
            tried++;
            int n = entry.final_IF.rows();
            Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
            new_IF.block(0,0,n,n) = entry.final_IF;
            int ext_idx = n;
            new_IF(ext_idx, ext_idx) = -3;  // so7 sits on (-3)
            new_IF(ext_idx, m2_curves[a]) = 1; new_IF(m2_curves[a], ext_idx) = 1;
            new_IF(ext_idx, m2_curves[b]) = 1; new_IF(m2_curves[b], ext_idx) = 1;

            if (sig_pos_exceeds_one_fast(new_IF, 193)) continue;
            auto sig = compute_sig(new_IF);
            if (sig.sig_pos != 1) continue;

            NHCResult nhc = check_nhc(new_IF, cc_budget);
            if (!nhc.passes) continue;
            for (auto& ei : entry.externals) {
                if (ei.is_hat1) enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                else if (ei.ext_si == -2 && ei.target_si == -1)
                    enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
            }
            if (!nhc.passes) continue;

            AnomalyResult anom = compute_anomaly(nhc, sig);
            bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
            auto result = make_result(entry, "so7", 0,
                ext_idx, -3, -2, 1, false, new_IF, nhc, anom, sig, cc_budget);
            if (is_sugra) { results.push_back(std::move(result)); passed++; }
            else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, result.externals)) g_nonsugra_results->push_back(std::move(result));
        }
    }
}

// so7mix: (-3) → two (-2) curves + (-1) curves simultaneously
inline void attach_so7mix(
    const CatEntry& entry, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    std::vector<int> m1_curves, m2_curves;
    for (auto& rt : entry.remaining) {
        if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
        else if (rt.self_int == -2) m2_curves.push_back(rt.curve_idx);
    }
    if ((int)m2_curves.size() < 2 || m1_curves.empty()) return;

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    int int_max = 4;
    for (int a = 0; a < (int)m2_curves.size(); a++) {
        for (int b = a + 1; b < (int)m2_curves.size(); b++) {
            int max_m1 = std::min((int)m1_curves.size(), 4);
            enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                enumerate_int_nums((int)m1_targets.size(), int_max, [&](const std::vector<int>& int_nums) {
                    tried++;
                    int n = entry.final_IF.rows();
                    Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                    new_IF.block(0,0,n,n) = entry.final_IF;
                    int ext_idx = n;
                    new_IF(ext_idx, ext_idx) = -3;
                    new_IF(ext_idx, m2_curves[a]) = 1; new_IF(m2_curves[a], ext_idx) = 1;
                    new_IF(ext_idx, m2_curves[b]) = 1; new_IF(m2_curves[b], ext_idx) = 1;
                    for (int i = 0; i < (int)m1_targets.size(); i++) {
                        new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                        new_IF(m1_targets[i], ext_idx) = int_nums[i];
                    }

                    if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
                    auto sig = compute_sig(new_IF);
                    if (sig.sig_pos != 1) return;

                    NHCResult nhc = check_nhc(new_IF, cc_budget);
                    if (!nhc.passes) return;
                    for (auto& ei : entry.externals) {
                        if (ei.is_hat1) enhance_hat1_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                        else if (ei.ext_si == -2 && ei.target_si == -1)
                            enhance_external_m2_gauge(nhc, new_IF, ei.curve_idx, ei.target_si, cc_budget);
                    }
                    if (!nhc.passes) return;

                    AnomalyResult anom = compute_anomaly(nhc, sig);
                    bool is_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, hat1_set);
                    auto result = make_result(entry, "so7mix", 0,
                        ext_idx, -3, -2, 1, false, new_IF, nhc, anom, sig, cc_budget);
                    if (is_sugra) { results.push_back(std::move(result)); passed++; }
                    else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, result.externals)) g_nonsugra_results->push_back(std::move(result));
                });
            });
        }
    }
}

// ============================================================================
// PART 8: Mode C — hat1 multi-target (hat1 → multiple -1 curves)
// ============================================================================

inline void attach_hat1_multi_target(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
    if (!spec.is_hat1 || spec.target_si != -1) return;

    // hat1→(-1): su(8) cc=7.0. Only (-1) with available_cc >= 7.0
    double hat1_cc = central_charge(GAUGE_SU8, 1);
    std::vector<int> m1_curves;
    for (auto& rt : entry.remaining)
        if (rt.self_int == -1 && rt.available_cc >= hat1_cc - 1e-9)
            m1_curves.push_back(rt.curve_idx);
    if ((int)m1_curves.size() < 2) return;

    int max_tgt = std::min((int)m1_curves.size(), 6);

    std::set<int> hat1_set;
    for (auto& ei : entry.externals)
        if (ei.is_hat1) hat1_set.insert(ei.curve_idx);

    enumerate_subsets_v2(m1_curves, 2, max_tgt, [&](const std::vector<int>& targets) {
        tried++;
        int n = entry.final_IF.rows();
        Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
        new_IF.block(0,0,n,n) = entry.final_IF;
        int ext_idx = n;
        new_IF(ext_idx, ext_idx) = -1;  // hat1: C² = -1
        for (int t : targets) { new_IF(ext_idx, t) = 1; new_IF(t, ext_idx) = 1; }

        if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
        auto sig = compute_sig(new_IF);

        // hat1: unimodular cannot reach Hirzebruch
        if (std::abs(sig.det) == 1) return;
        // hat1: skip |det|=n² filter (not Hirzebruch)

        NHCResult nhc = check_nhc(new_IF, cc_budget);
        if (!nhc.passes) return;
        if (!enhance_all_externals(nhc, new_IF, entry.externals,
                ext_idx, -1, -1, true, cc_budget))
            return;

        AnomalyResult anom = compute_anomaly(nhc, sig);
        if (total_ext_cc(entry.externals, -1, -1, 1, true) > g_max_ext_cc + 1e-9) return;

        std::set<int> h1_all = hat1_set; h1_all.insert(ext_idx);
        bool h1_sugra = (anom.H_neutral >= 0) && !reject_nonuni(new_IF, sig, anom, h1_all);
        auto h1_result = make_result(entry, spec.tag, spec.id,
            ext_idx, -1, -1, 1, true, new_IF, nhc, anom, sig, cc_budget);
        if (h1_sugra) { results.push_back(std::move(h1_result)); passed++; }
        else if (g_save_nonsugra && g_nonsugra_results && should_save_nonsugra(anom, h1_result.externals)) g_nonsugra_results->push_back(std::move(h1_result));
    });
}

// ============================================================================
// PART 9: Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_cat_dir> [output_suffix] [--det-sq]\n";
        return 1;
    }
    std::string input_dir = argv[1];
    std::string suffix = "";
    for (int i = 2; i < argc; i++) {
        if (std::string(argv[i]) == "--det-sq") g_det_sq_mode = true;
        else if (std::string(argv[i]) == "--use-lst-T") g_use_lst_T = true;
        else if (std::string(argv[i]) == "--no-su2") g_no_su2 = true;
        else if (std::string(argv[i]) == "--save-nonsugra") g_save_nonsugra = true;
        else if (suffix.empty()) suffix = argv[i];
    }

    std::cout << "=== gen_sugra_phase2 (v6 + hat1 multi) ===\n"
              << "Input: " << input_dir
              << (g_det_sq_mode ? "  [|det|=n² filter ON]" : "") << "\n\n";

    auto all_specs = build_all_specs();
    double cc_budget = 8.0;

    auto input_entries = load_all_cats(input_dir);
    std::cout << "Loaded: " << input_entries.size() << " entries\n\n";

    int n_ext_in = input_entries.empty() ? 0 : (int)input_entries[0].externals.size();
    std::string out_dir = "cat_" + std::to_string(n_ext_in + 1) + "ext";
    if (!suffix.empty()) out_dir += "_" + suffix;
    std::filesystem::remove_all(out_dir);
    std::filesystem::create_directories(out_dir);

    std::vector<CatEntry> results;
    int tried_a = 0, passed_a = 0;
    int tried_v6 = 0, passed_v6 = 0;
    int tried_h1m = 0, passed_h1m = 0;
    int tried_mix = 0, passed_mix = 0;
    int tried_so7 = 0, passed_so7 = 0;

    // Non-sugra results (anomaly/reject_nonuni failed but sig/NHC passed)
    std::vector<CatEntry> nonsugra;
    if (g_save_nonsugra) g_nonsugra_results = &nonsugra;

    {
        for (auto& entry : input_entries) {
            g_current_base_T = (int)entry.final_IF.rows() - (int)entry.externals.size();

            for (auto& spec : all_specs) {
                if (g_no_su2 && spec.ext_si == -2 && spec.target_si == -1 && !spec.is_hat1) continue;
                attach_single_target(entry, spec, cc_budget, results, tried_a, passed_a);
                attach_v6_multi_target(entry, spec, cc_budget, results, tried_v6, passed_v6);
                attach_hat1_multi_target(entry, spec, cc_budget, results, tried_h1m, passed_h1m);
                if (!g_no_su2)
                    attach_mixed_multi_target(entry, spec, cc_budget, results, tried_mix, passed_mix);
            }

            // su3n2mix: (-3) → (-2) + (-1)s
            {
                std::vector<int> m1_curves, m2_curves;
                for (auto& rt : entry.remaining) {
                    if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
                    else if (rt.self_int == -2) m2_curves.push_back(rt.curve_idx);
                }
                for (int m2 : m2_curves)
                    attach_generic_mixed(entry, cc_budget, -3, m2, 1, m1_curves,
                                         "su3n2mix", results, tried_mix, passed_mix);
            }

            // so16n2mix: (-4) → (-2) + (-1)s
            {
                std::vector<int> m1_curves, m2_curves;
                for (auto& rt : entry.remaining) {
                    if (rt.self_int == -1) m1_curves.push_back(rt.curve_idx);
                    else if (rt.self_int == -2) m2_curves.push_back(rt.curve_idx);
                }
                for (int m2 : m2_curves)
                    attach_generic_mixed(entry, cc_budget, -4, m2, 2, m1_curves,
                                         "so16n2mix", results, tried_mix, passed_mix);
            }

            // so7: (-3) → two (-2) curves
            attach_so7(entry, cc_budget, results, tried_so7, passed_so7);

            // so7mix: (-3) → two (-2) + (-1)s
            attach_so7mix(entry, cc_budget, results, tried_so7, passed_so7);
        }
    }

    std::cout << "Single-target:    tried=" << tried_a << ", passed=" << passed_a << "\n"
              << "v6 multi-target:  tried=" << tried_v6 << ", passed=" << passed_v6 << "\n"
              << "hat1 multi-target: tried=" << tried_h1m << ", passed=" << passed_h1m << "\n"
              << "mixed multi-target: tried=" << tried_mix << ", passed=" << passed_mix << "\n"
              << "so7/so7mix:       tried=" << tried_so7 << ", passed=" << passed_so7 << "\n"
              << "Total raw: " << results.size() << "\n";

    // Dedup
    std::set<std::string> seen;
    std::vector<CatEntry> filtered;
    for (auto& r : results) {
        // Include catalog_id (LST source) so different build paths → different keys
        auto key = dedup_key_v2(r.final_IF, r.combo_key + "|cid:" + std::to_string(r.catalog_id));
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
    std::map<std::string, std::vector<const CatEntry*>> by_combo;
    for (auto& r : filtered) by_combo[r.combo_key].push_back(&r);
    int tw = 0;
    for (auto& [combo, entries] : by_combo) {
        std::ofstream out(out_dir + "/" + combo + ".cat");
        out << "# Phase2 catalog: " << combo << "\n# Externals: "
            << entries[0]->externals.size() << "\n# Total: " << entries.size() << "\n#\n";
        int eid = 0;
        for (auto* rp : entries) { write_cat_entry(out, eid++, *rp); tw++; }
    }

    // INDEX
    { std::ofstream idx(out_dir + "/INDEX.txt");
      idx << "# Phase2 Index | Input: " << input_dir << " | Total: " << tw
          << " bases, " << by_combo.size() << " combos\n#\n";
      for (auto& [c, e] : by_combo) idx << c << " | " << e.size() << "\n"; }

    int n_uni = 0, n_nonuni = 0;
    for (auto& r : filtered) { if (std::abs(r.det) == 1) n_uni++; else n_nonuni++; }
    std::cout << "\nWritten " << tw << " bases (" << by_combo.size() << " combos) to "
              << out_dir << "/\nUnimodular: " << n_uni << ", Non-unimodular: " << n_nonuni << "\n";

    // LaTeX (compact — just summary + first 500 entries)
    { std::string texfile = out_dir + "/results.tex";
      std::ofstream tex(texfile);
      tex << "\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n"
          << "\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n"
          << "\\usepackage[dvipsnames]{xcolor}\n\\setlength{\\parindent}{0pt}\n"
          << "\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n";
      int n_ext = filtered.empty() ? 0 : (int)filtered[0].externals.size();
      tex << "\\section*{" << n_ext << "-External SUGRA Bases (v6+hat1 multi)}\n\n"
          << "Total: " << filtered.size() << " (" << n_uni << " uni, " << n_nonuni << " non-uni). "
          << "Modes: A=" << passed_a << ", B(v6)=" << passed_v6 << ", C(hat1m)=" << passed_h1m << ".\n\n";
      std::map<int,int> Tc; for (auto& r : filtered) Tc[r.T]++;
      tex << "\\begin{center}\\begin{tabular}{c|c}\\toprule $T$ & Count \\\\\\midrule\n";
      for (auto& [T,c] : Tc) tex << T << " & " << c << " \\\\\n";
      tex << "\\bottomrule\\end{tabular}\\end{center}\n\n";
      std::string cur_combo; int cur_T = -1, shown = 0;
      for (auto& r : filtered) {
          if (shown >= 500) break;
          if (r.combo_key != cur_combo) { cur_combo = r.combo_key; cur_T = -1;
              tex << "\n\\subsection*{" << cur_combo << "}\n\n"; }
          if (r.T != cur_T) { cur_T = r.T; tex << "\\paragraph{$T=" << cur_T << "$}\n\n"; }
          int delta = r.H_charged - r.V + 29*r.T - 273;
          int T_max = std::min((273+r.V-r.H_charged)/29, 193);
          tex << "$T=" << r.T << ",\\;\\Delta=" << delta << ",\\;\\det=" << r.det
              << ",\\;(+,-,0)=(" << r.sig_pos << "," << r.sig_neg << "," << r.sig_zero << ")$\n\n";
          shown++;
      }
      if ((int)filtered.size() > 500)
          tex << "\\bigskip\\ldots " << filtered.size()-500 << " more omitted.\n\n";
      tex << "\\end{document}\n";
      std::cout << "Written LaTeX to " << texfile << "\n";
    }

    // Write non-SUGRA catalog if enabled
    if (g_save_nonsugra && !nonsugra.empty()) {
        // Dedup nonsugra
        std::set<std::string> ns_seen;
        std::vector<CatEntry> ns_filtered;
        for (auto& r : nonsugra) {
            // Include catalog_id (LST source) so different build paths → different keys
        auto key = dedup_key_v2(r.final_IF, r.combo_key + "|cid:" + std::to_string(r.catalog_id));
            if (ns_seen.insert(key).second) ns_filtered.push_back(std::move(r));
        }
        std::sort(ns_filtered.begin(), ns_filtered.end(), [](const CatEntry& a, const CatEntry& b) {
            if (a.combo_key != b.combo_key) return a.combo_key < b.combo_key;
            return a.T < b.T;
        });

        std::string ns_dir = out_dir + "_nonsugra";
        std::filesystem::remove_all(ns_dir);
        std::filesystem::create_directories(ns_dir);

        std::map<std::string, std::vector<const CatEntry*>> ns_by_combo;
        for (auto& r : ns_filtered) ns_by_combo[r.combo_key].push_back(&r);
        int ns_tw = 0;
        for (auto& [combo, entries] : ns_by_combo) {
            std::ofstream out2(ns_dir + "/" + combo + ".cat");
            out2 << "# Non-SUGRA catalog: " << combo << "\n#\n";
            int eid = 0;
            for (auto* rp : entries) { write_cat_entry(out2, eid++, *rp); ns_tw++; }
        }
        std::cout << "Non-SUGRA: " << ns_tw << " bases (" << ns_by_combo.size()
                  << " combos) to " << ns_dir << "/\n";
    }

    return 0;
}
