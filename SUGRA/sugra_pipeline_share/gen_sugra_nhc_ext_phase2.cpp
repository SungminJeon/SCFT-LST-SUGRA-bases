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
    struct RemTarget { int curve_idx, self_int; double available_cc; };
    std::vector<RemTarget> remaining;
};

static std::vector<CatEntry>* g_nonsugra_results = nullptr;

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
        e.final_IF = Eigen::MatrixXi::Zero(n, n);
        for (int i = 0; i < n; i++) {
            std::getline(fin, line);
            std::istringstream ss(line);
            for (int j = 0; j < n; j++) ss >> e.final_IF(i, j);
        }
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
// Dedup + remaining
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
    // Build sorted combo_key for dedup (order-independent, matching Phase 2)
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

inline bool has_hat1(const CatEntry& entry, bool new_is_hat1 = false) {
    if (new_is_hat1) return true;
    for (auto& ei : entry.externals) if (ei.is_hat1) return true;
    return false;
}

// ============================================================================
// Mode A: Single-curve external attachment
// ============================================================================

inline void attach_single(
    const CatEntry& entry, const ExternalSpec& spec, double cc_budget,
    std::vector<CatEntry>& results, int& tried, int& passed)
{
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
        if (spec.is_hat1 && std::abs(sig.det) == 1) continue;
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
            double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
            if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) continue;
        }

        AnomalyResult anom = compute_anomaly(new_IF, nhc);
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
        max_k = 4;
        max_tgt = (int)m1_curves.size();
    } else {
        max_k = std::min(max_level_for_cc(ext_gauge, 8.0), 9);
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


            if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
            auto sig = compute_sig(new_IF);
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
                double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
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
            new_IF.block(0, 0, n, n) = entry.final_IF;

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

            // Quick cc check
            if (!quick_cc_check(new_IF, cc_budget)) return;

            // LDLT → eigensolve (cheap filters first)
            if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
            auto sig = compute_sig(new_IF);
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
                double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
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
    // c_ext pre-check
    if (total_ext_cc(entry.externals, ext_si, -1, 1, false) > g_max_ext_cc + 1e-9) return;

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
            // Connect to special target
            new_IF(ext_idx, tgt_curve) = tgt_int;
            new_IF(tgt_curve, ext_idx) = tgt_int;
            // Connect to (-1)s
            for (int i = 0; i < (int)m1_targets.size(); i++) {
                new_IF(ext_idx, m1_targets[i]) = int_nums[i];
                new_IF(m1_targets[i], ext_idx) = int_nums[i];
            }

            if (sig_pos_exceeds_one_fast(new_IF, 193)) return;
            auto sig = compute_sig(new_IF);
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
                double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
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
    if (m1_curves.empty()) return;

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

    // so7: (-3) → two (-2) curves (forms (-2)-(-3)-(-2) = su2+so7+su2)
    if ((int)m2_curves.size() >= 2) {
        std::set<int> hat1_set;
        for (auto& ei : entry.externals)
            if (ei.is_hat1) hat1_set.insert(ei.curve_idx);
        for (int a = 0; a < (int)m2_curves.size(); a++) {
            for (int b = a + 1; b < (int)m2_curves.size(); b++) {
                tried++;
                int n = entry.final_IF.rows();
                Eigen::MatrixXi new_IF = Eigen::MatrixXi::Zero(n+1, n+1);
                new_IF.block(0,0,n,n) = entry.final_IF;
                int ext_idx = n;
                new_IF(ext_idx, ext_idx) = -3;
                new_IF(ext_idx, m2_curves[a]) = 1; new_IF(m2_curves[a], ext_idx) = 1;
                new_IF(ext_idx, m2_curves[b]) = 1; new_IF(m2_curves[b], ext_idx) = 1;

                if (sig_pos_exceeds_one_fast(new_IF, 193)) continue;
                auto sig = compute_sig(new_IF);
                if (sig.sig_pos != 1) continue;

                NHCResult nhc = check_nhc(new_IF, cc_budget);
                if (!nhc.passes) continue;
                if (!enhance_all_externals(nhc, new_IF, entry.externals, cc_budget)) continue;

                // True c_ext ≤ 16 — REJECT
                {
                    std::set<int> ext_set;
                    for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);
                    ext_set.insert(ext_idx);
                    double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
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
    if ((int)m2_curves.size() >= 2 && !m1_curves.empty()) {
        std::set<int> hat1_set;
        for (auto& ei : entry.externals)
            if (ei.is_hat1) hat1_set.insert(ei.curve_idx);
        for (int a = 0; a < (int)m2_curves.size(); a++) {
            for (int b = a + 1; b < (int)m2_curves.size(); b++) {
                int max_m1 = std::min((int)m1_curves.size(), 4);
                enumerate_subsets_v2(m1_curves, 1, max_m1, [&](const std::vector<int>& m1_targets) {
                    enumerate_int_nums((int)m1_targets.size(), 4, [&](const std::vector<int>& int_nums) {
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
                        if (!enhance_all_externals(nhc, new_IF, entry.externals, cc_budget)) return;

                        // True c_ext ≤ 16 — REJECT
                        {
                            std::set<int> ext_set;
                            for (auto& ei : entry.externals) ext_set.insert(ei.curve_idx);
                            ext_set.insert(ext_idx);
                            double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
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
        else if (std::string(argv[i]) == "--use-lst-T") g_use_lst_T = true;
        else if (std::string(argv[i]) == "--save-nonsugra") g_save_nonsugra = true;
        else if (std::string(argv[i]) == "--no-su2") g_no_su2 = true;
        else if (std::string(argv[i]) == "--no-223") g_no_223 = true;
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
    if (g_save_nonsugra) g_nonsugra_results = &nonsugra;

    for (size_t ei = 0; ei < input_entries.size(); ei++) {
        auto& entry = input_entries[ei];
        g_current_base_T = (int)entry.final_IF.rows() - (int)entry.externals.size();

        // Single-curve externals (no ordering constraint — try all specs)
        for (auto& spec : all_specs) {
            if (g_no_su2 && spec.ext_si == -2 && spec.target_si == -1 && !spec.is_hat1) continue;
            attach_single(entry, spec, cc_budget, results, tried_single, passed_single);
            attach_v6_multi(entry, spec, cc_budget, results, tried_v6, passed_v6);
        }

        // Mixed multi-target (all variants: su2n3mix, su8mix, su3n2mix, so16n2mix, so7, so7mix)
        attach_mixed_multi(entry, cc_budget, results, tried_mix, passed_mix);

        // NHC cluster externals
        for (auto& ns : nhc_cluster_specs) {
            if (g_no_223 && ns.tag == "nhc_2_2_3") continue;
            attach_nhc_cluster(entry, ns, cc_budget, results, tried_nhc, passed_nhc);
        }

        if ((ei + 1) % 500 == 0)
            std::cout << "  processed " << (ei + 1) << "/" << input_entries.size()
                      << " entries, results so far: " << results.size() << "\n";
    }

    std::cout << "\nSingle-target:    tried=" << tried_single << ", passed=" << passed_single << "\n"
              << "v6 multi-target:  tried=" << tried_v6 << ", passed=" << passed_v6 << "\n"
              << "NHC cluster:      tried=" << tried_nhc << ", passed=" << passed_nhc << "\n"
              << "Mixed multi:      tried=" << tried_mix << ", passed=" << passed_mix << "\n"
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
        out << "# NHC ext Phase2 catalog: " << combo << "\n# Round: " << round_str
            << "\n# Total: " << entries.size() << "\n#\n";
        int eid = 0;
        for (auto* rp : entries) { write_cat_entry(out, eid++, *rp); tw++; }
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
        std::set<std::string> ns_seen;
        std::vector<CatEntry> ns_filtered;
        for (auto& r : nonsugra) {
            // Include catalog_id (LST source) so different build paths → different keys
        auto key = dedup_key_v2(r.final_IF, r.combo_key + "|cid:" + std::to_string(r.catalog_id));
            if (ns_seen.insert(key).second) ns_filtered.push_back(std::move(r));
        }
        std::string ns_dir = out_dir + "_nonsugra";
        std::filesystem::remove_all(ns_dir);
        std::filesystem::create_directories(ns_dir);
        std::map<std::string, std::vector<const CatEntry*>> ns_by_combo;
        for (auto& r : ns_filtered) ns_by_combo[r.combo_key].push_back(&r);
        int ns_tw = 0;
        for (auto& [combo, entries] : ns_by_combo) {
            std::ofstream out2(ns_dir + "/" + combo + ".cat");
            out2 << "# Non-SUGRA: " << combo << "\n#\n";
            int eid = 0;
            for (auto* rp : entries) { write_cat_entry(out2, eid++, *rp); ns_tw++; }
        }
        std::cout << "Non-SUGRA: " << ns_tw << " bases (" << ns_by_combo.size()
                  << " combos) to " << ns_dir << "/\n";
    }

    return 0;
}
