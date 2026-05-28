// gen_sugra_nhc_ext.cpp — NHC cluster external attachment
//
// Attaches multi-curve NHC clusters (e.g., -2-3, -2-3-2, -2-2-3) as externals
// to LST bases. Each curve in the NHC connects to base curve(s).
//
// Usage: ./gen_sugra_nhc_ext <catalog> [T_max=10] [T_min=0]

#include "sugra_generator.h"
#include <fstream>
#include <cstdlib>
#include <set>
#include <iomanip>
#include <filesystem>

// ============================================================================
// NHC External Spec
// ============================================================================

struct NHCExtSpec {
    int id;
    std::string tag;
    std::string label;
    std::vector<int> self_ints;              // e.g., {-2, -3}
    struct InternalInt { int a, b, k; };     // curve a ↔ curve b with int=k
    std::vector<InternalInt> internal_ints;
    int H, V;   // NHC gauge contributions
    double cc_total;  // total central charge of this cluster
};

inline std::vector<NHCExtSpec> build_nhc_ext_specs() {
    std::vector<NHCExtSpec> specs;
    int id = 0;

    // -2-3 chain: su(2) x G2, H=8, V=17, cc=1.0+2.8=3.8
    specs.push_back({id++, "nhc_2_3", "(-2)-(-3) NHC",
                     {-2, -3}, {{0, 1, 1}}, 8, 17, 3.8});

    // -2-3-2 chain: su(2) x so(7) x su(2), H=16, V=27, cc=1.0+3.5+1.0=5.5
    specs.push_back({id++, "nhc_2_3_2", "(-2)-(-3)-(-2) NHC",
                     {-2, -3, -2}, {{0, 1, 1}, {1, 2, 1}}, 16, 27, 5.5});

    // -2-2-3 chain: . x sp(1) x G2, H=8, V=17, cc=0+1.0+2.8=3.8
    // Curve order: [-3, -2, -2] (so -3 processed first → cc-prune early)
    specs.push_back({id++, "nhc_2_2_3", "(-2)-(-2)-(-3) NHC",
                     {-3, -2, -2}, {{0, 1, 1}, {1, 2, 1}}, 8, 17, 3.8});

    // -2-4 chain: su(8) x so(16), H=128, V=183, cc=12.6+15.0=27.6
    specs.push_back({id++, "nhc_2_4", "(-2)-(-4) NHC",
                     {-2, -4}, {{0, 1, 2}}, 128, 183, 27.6});

    return specs;
}

// ============================================================================
// Attachment: add NHC cluster to base IF
// ============================================================================

// Attach NHC external to base IF.
// target_map[i] = {(base_curve_idx, int_num), ...} for NHC curve i
// Returns new IF with NHC curves appended.
inline Eigen::MatrixXi attach_nhc(
    const Eigen::MatrixXi& base_IF,
    const NHCExtSpec& spec,
    const std::vector<std::vector<std::pair<int,int>>>& target_map)
{
    int nb = base_IF.rows();
    int ne = (int)spec.self_ints.size();
    int n = nb + ne;

    Eigen::MatrixXi IF = Eigen::MatrixXi::Zero(n, n);
    IF.block(0, 0, nb, nb) = base_IF;

    // NHC internal structure
    for (int i = 0; i < ne; i++)
        IF(nb + i, nb + i) = spec.self_ints[i];
    for (auto& ii : spec.internal_ints) {
        IF(nb + ii.a, nb + ii.b) = ii.k;
        IF(nb + ii.b, nb + ii.a) = ii.k;
    }

    // NHC → base connections
    for (int i = 0; i < ne; i++) {
        for (auto& [tgt, inum] : target_map[i]) {
            IF(nb + i, tgt) += inum;
            IF(tgt, nb + i) += inum;
        }
    }

    return IF;
}

// ============================================================================
// Result struct
// ============================================================================

struct NHCExtResult {
    int nhc_spec_id;
    int catalog_id;
    std::string catalog_type;
    int base_T;

    std::vector<int> ext_curve_indices;   // indices of NHC curves in final_IF
    std::vector<std::vector<std::pair<int,int>>> target_map;  // per NHC curve: (base_idx, int_num)

    Eigen::MatrixXi base_IF;
    Eigen::MatrixXi final_IF;

    NHCResult nhc;
    AnomalyResult anomaly;
    SigInfo sig;
    std::vector<double> eigenvalues;
};

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <catalog> [T_max=10] [T_min=0] [--det-sq]\n";
        return 1;
    }

    std::string catalog_file = argv[1];
    int T_max = (argc > 2) ? std::atoi(argv[2]) : 10;
    int T_min = (argc > 3) ? std::atoi(argv[3]) : 0;
    bool det_sq_mode = false;
    bool save_nonsugra = false;
    bool skip_223 = false;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--det-sq") det_sq_mode = true;
        if (std::string(argv[i]) == "--use-lst-T") g_use_lst_T = true;
        if (std::string(argv[i]) == "--save-nonsugra") save_nonsugra = true;
        if (std::string(argv[i]) == "--no-223") skip_223 = true;
    }

    std::cout << "=== gen_sugra_nhc_ext ===\n";
    std::cout << "Catalog: " << catalog_file << "\n";
    std::cout << "T range: [" << T_min << ", " << T_max << "]\n\n";

    auto catalog = load_catalog(catalog_file);
    std::cout << "Catalog loaded: " << catalog.size() << " entries\n";

    auto nhc_specs = build_nhc_ext_specs();
    std::cout << "NHC external specs: " << nhc_specs.size() << " types\n\n";

    SUGRAConfig config;
    config.T_max = 193;
    config.check_nhc = true;
    config.check_anomaly = true;
    config.check_determinant = false;
    config.check_det_perfect_square = det_sq_mode;

    // Add DM (dummy) bases to catalog
    {
        // DM:A(n): linear chain of n curves, all -1 except ends, T = n-1
        int a_min = std::max(2, T_min + 1);
        int a_max = T_max + 1;
        for (int nc = a_min; nc <= a_max; nc++) {
            CatalogEntry dm;
            dm.id = -nc;
            dm.type = "DM:A(" + std::to_string(nc) + ")";
            dm.T = nc - 1;
            dm.topo = "";
            catalog.push_back(dm);
        }
        // DM:D(n): D-type
        int d_min = std::max(3, T_min + 1);
        int d_max = T_max + 1;
        for (int nc = d_min; nc <= d_max; nc++) {
            CatalogEntry dm;
            dm.id = -(1000 + nc);
            dm.type = "DM:D(" + std::to_string(nc) + ")";
            dm.T = nc - 1;
            dm.topo = "";
            catalog.push_back(dm);
        }
        std::cout << "Added DM bases: A(" << a_min << "~" << a_max << "), D(" << d_min << "~" << d_max << ")\n";
    }

    // Precompute IFs
    struct CachedIF { Eigen::MatrixXi IF; SigInfo sig; bool valid = false; };
    std::vector<CachedIF> cached_ifs(catalog.size());
    int n_cached = 0;
    for (size_t i = 0; i < catalog.size(); i++) {
        auto& entry = catalog[i];
        cached_ifs[i].valid = false;
        if (entry.T < T_min || entry.T > T_max) continue;
        // DM bases: build IF directly
        if (entry.id < 0 && entry.id > -1000) {
            // DM:A(n)
            int nc = -entry.id;
            cached_ifs[i].IF = build_dummy_A(nc);
        } else if (entry.id <= -1000) {
            // DM:D(n)
            int nc = -(entry.id + 1000);
            cached_ifs[i].IF = build_dummy_D(nc);
        } else {
            cached_ifs[i].IF = reconstruct_IF(entry);
        }
        if (cached_ifs[i].IF.rows() == 0) continue;
        cached_ifs[i].sig = compute_sig(cached_ifs[i].IF);
        cached_ifs[i].valid = true;
        n_cached++;
    }
    std::cout << "Precomputed IFs: " << n_cached << " entries (incl. DM)\n\n";

    std::string dir = "cat_nhc_ext";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    int total_bases = 0;
    int total_nonsugra = 0;

    for (auto& nhc_spec : nhc_specs) {
        if (skip_223 && nhc_spec.tag == "nhc_2_2_3") {
            std::cout << "[" << nhc_spec.id << "] " << nhc_spec.label << " (" << nhc_spec.tag << ") SKIPPED (--no-223)\n";
            continue;
        }
        std::cout << "[" << nhc_spec.id << "] " << nhc_spec.label << " (" << nhc_spec.tag << ")";
        std::cout.flush();

        int ne = (int)nhc_spec.self_ints.size();
        std::set<std::string> seen;          // SUGRA dedup
        std::set<std::string> nonsugra_seen; // non-SUGRA dedup (separate so cospectral pairs aren't merged)
        std::vector<NHCExtResult> results;
        std::vector<NHCExtResult> nonsugra_results;

        for (size_t ci = 0; ci < catalog.size(); ci++) {
            if (!cached_ifs[ci].valid) continue;
            auto& entry = catalog[ci];
            auto& base_IF = cached_ifs[ci].IF;
            int nb = base_IF.rows();
            g_current_base_T = nb;  // LST curve count for --use-lst-T

            // Find (-1) curves with cc budget room + (-2) curves as targets
            std::vector<int> m1_curves;  // (-1) with cc room
            std::vector<int> m2_curves;  // (-2) curves
            for (int i = 0; i < nb; i++) {
                if (base_IF(i, i) == -1) {
                    double used = used_cc_on_minus1(base_IF, i);
                    if (used + 1.0 <= config.cc_budget + 1e-9)
                        m1_curves.push_back(i);
                } else if (base_IF(i, i) == -2) {
                    m2_curves.push_back(i);
                }
            }

            if (m1_curves.empty()) continue;

            // Combined targets: (-1) and (-2) curves (only for nhc_2_3)
            std::vector<int> all_targets = m1_curves;
            if (nhc_spec.tag == "nhc_2_3")
                all_targets.insert(all_targets.end(), m2_curves.begin(), m2_curves.end());

            // For simplicity: enumerate assignments of NHC curves to base (-1) curves
            // Each NHC curve connects to one base (-1) curve with int=1
            // Allow same base curve for multiple NHC curves

            // For 2-curve NHC: enumerate all pairs (t0, t1) from m1_curves (allow repeats)
            // For 3-curve NHC: enumerate all triples (t0, t1, t2)
            // Cap: if too many combos, sample

            // Pre-compute available cc per (-1) curve
            std::map<int, double> avail_cc;
            for (int i : m1_curves)
                avail_cc[i] = config.cc_budget - used_cc_on_minus1(base_IF, i);

            // assignment[i] = list of (target, int_num) for NHC curve i
            // cc_used tracks incremental cc usage per (-1) curve from assignments so far
            std::map<int, double> cc_used;  // curve_idx -> accumulated cc from NHC assignments

            std::function<void(int, std::vector<std::vector<std::pair<int,int>>>&)> enumerate;
            enumerate = [&](int nhc_idx, std::vector<std::vector<std::pair<int,int>>>& assignment) {
                if (nhc_idx == ne) {
                    auto& target_map = assignment;
                    Eigen::MatrixXi new_IF = attach_nhc(base_IF, nhc_spec, target_map);

                    // Quick cc check (full, catches non-(-1) interactions too)
                    if (!quick_cc_check(new_IF, config.cc_budget)) return;

                    // Fast signature prefilter (hybrid LDLT + eigensolve fallback): rejects
                    // the vast majority of leaves cheaply before the full eigensolve below.
                    if (sig_pos_exceeds_one_fast(new_IF, config.T_max)) return;

                    // Full eigensolve
                    std::vector<double> ev;
                    auto sig = compute_sig(new_IF, ev);
                    if (sig.sig_pos != 1) return;

                    // |det|=n² filter
                    if (config.check_det_perfect_square && !is_perfect_square(std::abs(sig.det))) return;

                    // NHC check
                    NHCResult nhc = check_nhc(new_IF, config.cc_budget);
                    if (!nhc.passes) return;

                    // Proper c_ext ≤ 16 prune (NHC-aware bridge-level)
                    std::set<int> ext_set;
                    for (int i = 0; i < ne; i++) ext_set.insert(nb + i);
                    if (compute_proper_c_ext(new_IF, nhc, ext_set) > g_proper_c_ext_max + 1e-9)
                        return;

                    // True (fiber-weighted) c_ext ≤ 16 — REJECT (no recovery: c_ext only grows
                    // with further attachments, so non-SUGRA chain can't rescue this)
                    {
                        double tc = compute_true_c_ext_fiber(new_IF, nhc, ext_set);
                        if (tc >= 0 && tc > g_true_c_ext_max + 1e-9) return;
                    }

                    // Anomaly + reject_nonuni
                    AnomalyResult anom = compute_anomaly(nhc, sig);
                    bool is_sugra = true;
                    if (anom.H_neutral < 0) is_sugra = false;
                    if (is_sugra) {
                        std::set<int> h1;
                        if (reject_nonuni(new_IF, sig, anom, h1)) is_sugra = false;
                    }

                    // Dedup: eigenvalues + diagonal — separate per sugra/non-sugra
                    // (cospectral non-isomorphic IFs can fall in same key; the SUGRA
                    //  representative must be kept regardless of enumeration order.)
                    std::string dedup_key;
                    {
                        std::ostringstream ks;
                        ks << std::fixed << std::setprecision(6);
                        auto sev = ev; std::sort(sev.begin(), sev.end());
                        for (auto e : sev) ks << (std::abs(e) < 1e-10 ? 0.0 : e) << ",";
                        ks << "|";
                        int nn = new_IF.rows();
                        for (int ii = 0; ii < nn; ii++) ks << new_IF(ii,ii) << ",";
                        ks << "|" << nhc_spec.tag;
                        // Include catalog_id (LST source) so different build paths → different keys
                        ks << "|cid:" << entry.id;
                        dedup_key = ks.str();
                    }
                    if (is_sugra) {
                        if (!seen.insert(dedup_key).second) return;
                    } else {
                        if (!nonsugra_seen.insert(dedup_key).second) return;
                    }

                    NHCExtResult r;
                    r.nhc_spec_id = nhc_spec.id;
                    r.catalog_id = entry.id;
                    r.catalog_type = entry.type;
                    r.base_T = entry.T;
                    r.ext_curve_indices.resize(ne);
                    for (int i = 0; i < ne; i++) r.ext_curve_indices[i] = nb + i;
                    r.target_map = target_map;
                    r.base_IF = base_IF;
                    r.final_IF = new_IF;
                    r.nhc = nhc;
                    r.anomaly = anom;
                    r.sig = sig;
                    r.eigenvalues = ev;
                    if (is_sugra) {
                        results.push_back(std::move(r));
                    } else if (save_nonsugra) {
                        if (!recovery_impossible(anom.H_neutral, 16.0 - nhc_spec.cc_total))
                            nonsugra_results.push_back(std::move(r));
                    }
                    return;
                }

                // Generate all possible (target, k) assignments for this NHC curve
                // Single-target: one (target, k) pair
                // Multi-target: 2+ different (-1) curves with k=1 or 2 each

                // Helper: check if assignment fits cc budget on all (-1) targets
                // Returns false if any target (-1) would exceed cc_budget
                int nhc_si = nhc_spec.self_ints[nhc_idx];
                // NHC-aware gauge for cc accounting. A (-2) inside an NHC chain carries
                // SU2 (NOT none), so it DOES consume cc on its (-1) target. gauge_from_si(-2)
                // returns NONE -> add=0 -> incremental prune never fires -> multi-target
                // enumeration explodes (esp. nhc_2_3_2). Charging the true NHC gauge here makes
                // the prune fire and matches check_nhc's per-(-1) accounting (output-preserving).
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
                    // Check cc on each (-1) target
                    for (auto& [t, k] : tgts) {
                        if (base_IF(t, t) != -1) continue;
                        double add = (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
                        if (cc_used[t] + add > avail_cc[t] + 1e-9) return;
                    }
                    // Commit cc
                    for (auto& [t, k] : tgts) {
                        if (base_IF(t, t) != -1) continue;
                        cc_used[t] += (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
                    }
                    assignment.push_back(tgts);
                    enumerate(nhc_idx + 1, assignment);
                    assignment.pop_back();
                    // Undo cc
                    for (auto& [t, k] : tgts) {
                        if (base_IF(t, t) != -1) continue;
                        cc_used[t] -= (nhc_g.dim > 0) ? central_charge(nhc_g, k) : 0;
                    }
                };

                // Compute this NHC curve's minimum cc (= NHC-aware gauge at k=1).
                // Use this to pre-filter targets that can't possibly accept this curve.
                double curve_min_cc = (nhc_g.dim > 0) ? central_charge(nhc_g, 1) : 0.0;
                if (nhc_spec.tag == "nhc_2_3_2") {
                    static const double tbl[] = {1.0, 3.5, 1.0};  // SU2, SO7, SU2 at k=1
                    if (nhc_idx < 3) curve_min_cc = tbl[nhc_idx];
                } else if (nhc_spec.tag == "nhc_2_2_3") {
                    static const double tbl[] = {2.8, 1.0, 0.0};  // G2, SU2, NONE (spec order)
                    if (nhc_idx < 3) curve_min_cc = tbl[nhc_idx];
                } else if (nhc_spec.tag == "nhc_2_4") {
                    static const double tbl[] = {6.3, 15.0};  // SU8, SO16
                    if (nhc_idx < 2) curve_min_cc = tbl[nhc_idx];
                } else if (nhc_spec.tag == "nhc_2_3") {
                    int si = nhc_spec.self_ints[nhc_idx];
                    if (si == -2) curve_min_cc = 1.0;        // SU2
                    else if (si == -3) curve_min_cc = 2.8;   // G2
                }

                // Per-curve target list: only (-1) curves with avail_cc ≥ curve_min_cc
                std::vector<int> curve_targets;
                std::vector<int> curve_m1_targets;
                for (int t : all_targets) {
                    if (base_IF(t, t) == -1) {
                        if (avail_cc[t] >= curve_min_cc - 1e-9) {
                            curve_targets.push_back(t);
                            curve_m1_targets.push_back(t);
                        }
                    } else {
                        curve_targets.push_back(t);  // (-2) target for nhc_2_3
                    }
                }

                // nhc_2_2_3: restrict to single-target with k=1 only (combinatorics control)
                bool restrict_223 = (nhc_spec.tag == "nhc_2_2_3");

                // Single-target options
                for (int t : curve_targets) {
                    int kmax = restrict_223 ? 1 : (base_IF(t,t) == -1 ? 2 : 1);
                    for (int k = 1; k <= kmax; k++) {
                        try_assignment({{t, k}});
                    }
                }

                // Multi-target: from per-curve filtered list (skipped for nhc_2_2_3)
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

        // Sort
        std::sort(results.begin(), results.end(),
                  [](const NHCExtResult& a, const NHCExtResult& b) {
            if (a.anomaly.T != b.anomaly.T) return a.anomaly.T < b.anomaly.T;
            return std::abs(a.sig.det) < std::abs(b.sig.det);
        });

        // Write .cat file
        std::string fname = dir + "/" + nhc_spec.tag + ".cat";
        std::ofstream out(fname);
        out << "# NHC External catalog: " << nhc_spec.tag << "\n";
        out << "# NHC: " << nhc_spec.label << "\n";
        out << "# Total: " << results.size() << " bases\n#\n";

        for (int ri = 0; ri < (int)results.size(); ri++) {
            auto& r = results[ri];
            int n = r.final_IF.rows();
            int ne_r = (int)r.ext_curve_indices.size();

            out << "ENTRY " << ri << "\n";
            out << "COMBO " << nhc_spec.tag << "\n";
            out << "EXTERNALS " << ne_r << "\n";
            for (int i = 0; i < ne_r; i++) {
                // curve_idx spec_id tag ext_si target_si int_num is_hat1
                int tgt = r.target_map[i][0].first;
                int inum = r.target_map[i][0].second;
                int tgt_si = r.final_IF(tgt, tgt);
                out << r.ext_curve_indices[i] << " " << nhc_spec.id << " " << nhc_spec.tag
                    << " " << nhc_spec.self_ints[i] << " " << tgt_si << " " << inum << " 0\n";
            }
            out << "BASE " << r.catalog_id << " " << r.catalog_type
                << " " << r.base_T << " " << r.base_IF.rows() << "\n";
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
            // REMAINING: compute available cc for each non-external curve
            {
                std::set<int> ext_set(r.ext_curve_indices.begin(), r.ext_curve_indices.end());
                std::vector<std::tuple<int,int,double>> remaining;
                for (int ci = 0; ci < n; ci++) {
                    if (ext_set.count(ci)) continue;
                    int si = r.final_IF(ci, ci);
                    double avail = config.cc_budget;
                    if (si == -1) {
                        double used = 0;
                        for (int j = 0; j < n; j++) {
                            if (j == ci || r.final_IF(ci, j) == 0 || r.final_IF(j, j) == -1) continue;
                            GaugeInfo g = (j < (int)r.nhc.curve_gauges.size() && r.nhc.curve_gauges[j].dim > 0)
                                ? r.nhc.curve_gauges[j] : gauge_from_si(r.final_IF(j, j));
                            if (g.dim > 0) used += central_charge(g, std::abs(r.final_IF(ci, j)));
                        }
                        avail = config.cc_budget - used;
                    }
                    if (si == -1 && avail > 0.1)
                        remaining.push_back({ci, si, avail});
                    else if (si <= -2)
                        remaining.push_back({ci, si, avail});
                }
                out << "REMAINING " << remaining.size() << "\n";
                for (auto& [idx, si, cc] : remaining) {
                    out << idx << " " << si << " "
                        << std::fixed << std::setprecision(2) << cc << "\n";
                }
            }
            out << "END\n\n";
        }
        out.close();

        std::cout << " → " << results.size() << " bases";
        if (save_nonsugra && !nonsugra_results.empty())
            std::cout << " (nonsugra: " << nonsugra_results.size() << ")";
        std::cout << "\n";
        total_bases += (int)results.size();

        // Write non-SUGRA
        if (save_nonsugra && !nonsugra_results.empty()) {
            std::string ns_dir = dir + "_nonsugra";
            std::filesystem::create_directories(ns_dir);
            std::string ns_fname = ns_dir + "/" + nhc_spec.tag + ".cat";
            std::ofstream ns_out(ns_fname);
            ns_out << "# NHC External catalog (non-SUGRA): " << nhc_spec.tag << "\n";
            ns_out << "# Total: " << nonsugra_results.size() << " bases\n#\n";
            for (int ri = 0; ri < (int)nonsugra_results.size(); ri++) {
                auto& r = nonsugra_results[ri];
                int n = r.final_IF.rows();
                int ne_r = (int)r.ext_curve_indices.size();
                ns_out << "ENTRY " << ri << "\n";
                ns_out << "COMBO " << nhc_spec.tag << "\n";
                ns_out << "EXTERNALS " << ne_r << "\n";
                for (int i = 0; i < ne_r; i++) {
                    int tgt = r.target_map[i][0].first;
                    int inum = r.target_map[i][0].second;
                    int tgt_si = r.final_IF(tgt, tgt);
                    ns_out << r.ext_curve_indices[i] << " " << nhc_spec.id << " " << nhc_spec.tag
                        << " " << nhc_spec.self_ints[i] << " " << tgt_si << " " << inum << " 0\n";
                }
                ns_out << "BASE " << r.catalog_id << " " << r.catalog_type
                    << " " << r.base_T << " " << r.base_IF.rows() << "\n";
                ns_out << "PHYSICS " << r.anomaly.T << " " << r.anomaly.H_charged
                    << " " << r.anomaly.V << " " << r.anomaly.H_neutral
                    << " " << r.sig.det << " " << r.sig.sig_pos
                    << " " << r.sig.sig_neg << " " << r.sig.sig_zero << "\n";
                ns_out << "IF " << n << "\n";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) { if (j) ns_out << " "; ns_out << r.final_IF(i,j); }
                    ns_out << "\n";
                }
                {
                    std::set<int> ext_set(r.ext_curve_indices.begin(), r.ext_curve_indices.end());
                    std::vector<std::tuple<int,int,double>> remaining;
                    for (int ci = 0; ci < n; ci++) {
                        if (ext_set.count(ci)) continue;
                        int si = r.final_IF(ci, ci);
                        double avail = config.cc_budget;
                        if (si == -1) {
                            double used = 0;
                            for (int j = 0; j < n; j++) {
                                if (j == ci || r.final_IF(ci, j) == 0 || r.final_IF(j, j) == -1) continue;
                                GaugeInfo g = (j < (int)r.nhc.curve_gauges.size() && r.nhc.curve_gauges[j].dim > 0)
                                    ? r.nhc.curve_gauges[j] : gauge_from_si(r.final_IF(j, j));
                                if (g.dim > 0) used += central_charge(g, std::abs(r.final_IF(ci, j)));
                            }
                            avail = config.cc_budget - used;
                        }
                        if (si == -1 && avail > 0.1)
                            remaining.push_back({ci, si, avail});
                        else if (si <= -2)
                            remaining.push_back({ci, si, avail});
                    }
                    ns_out << "REMAINING " << remaining.size() << "\n";
                    for (auto& [idx, si, cc] : remaining)
                        ns_out << idx << " " << si << " " << std::fixed << std::setprecision(2) << cc << "\n";
                }
                ns_out << "END\n\n";
            }
            ns_out.close();
            total_nonsugra += (int)nonsugra_results.size();
        }
    }

    std::cout << "\n=== Summary ===\n";
    std::cout << "Total: " << total_bases << " NHC external bases\n";
    if (save_nonsugra) std::cout << "Non-SUGRA: " << total_nonsugra << " bases to " << dir << "_nonsugra/\n";
    std::cout << "Output: " << dir << "/\n";

    // LaTeX output
    {
        std::string texfile = dir + "/results.tex";
        std::ofstream tex(texfile);
        tex << "\\documentclass[10pt]{article}\n\\usepackage[a4paper,margin=0.6in]{geometry}\n"
            << "\\usepackage{longtable,booktabs,amsmath,amssymb,array}\n\\usepackage[dvipsnames]{xcolor}\n"
            << "\\setlength{\\parindent}{0pt}\n\\setlength{\\parskip}{2pt}\n\\begin{document}\n\n"
            << "\\section*{NHC External SUGRA Bases}\n\n"
            << "Total: " << total_bases << " bases.\n\n";

        for (auto& nhc_spec : nhc_specs) {
            std::string catf = dir + "/" + nhc_spec.tag + ".cat";
            std::ifstream fin(catf);
            if (!fin) continue;

            // Count entries
            std::string text((std::istreambuf_iterator<char>(fin)),
                             std::istreambuf_iterator<char>());
            fin.close();

            int n_entries = 0;
            { size_t pos = 0; while ((pos = text.find("ENTRY ", pos)) != std::string::npos) { n_entries++; pos++; } }

            tex << "\\subsection*{" << nhc_spec.label << " (" << nhc_spec.tag << ") --- " << n_entries << " bases}\n\n";

            // Parse and show entries
            int shown = 0;
            size_t pos = 0;
            while ((pos = text.find("ENTRY ", pos)) != std::string::npos && shown < 500) {
                size_t end = text.find("END", pos);
                if (end == std::string::npos) break;
                std::string block = text.substr(pos, end - pos);

                // Parse PHYSICS line
                size_t ppos = block.find("PHYSICS ");
                if (ppos != std::string::npos) {
                    std::istringstream ss(block.substr(ppos + 8));
                    int T, H, V, Hn, det, sp, sn, sz;
                    ss >> T >> H >> V >> Hn >> det >> sp >> sn >> sz;
                    int Tm = std::min((273 + V - H) / 29, 193);

                    // Parse IF
                    size_t ipos = block.find("IF ");
                    int n = 0;
                    if (ipos != std::string::npos) {
                        std::istringstream iss(block.substr(ipos + 3));
                        iss >> n;
                    }

                    // Get diagonal
                    std::vector<int> diag;
                    if (ipos != std::string::npos) {
                        std::istringstream iss(block.substr(ipos));
                        std::string line;
                        std::getline(iss, line); // skip "IF n"
                        for (int row = 0; row < n; row++) {
                            std::getline(iss, line);
                            std::istringstream rss(line);
                            for (int col = 0; col <= row; col++) {
                                int v; rss >> v;
                                if (col == row) diag.push_back(v);
                            }
                        }
                    }

                    tex << "$";
                    for (int d = 0; d < (int)diag.size(); d++) {
                        if (d) tex << "\\;";
                        tex << "(" << -diag[d] << ")";
                    }
                    tex << "$ \\quad $T=" << T << ",\\;\\det=" << det
                        << ",\\;H_n=" << Hn << ",\\;T_{\\max}=" << Tm << "$\n\n";

                    // Show IF as smallmatrix
                    if (ipos != std::string::npos && n > 0) {
                        tex << "\\hspace{2em}$\\left(\\begin{smallmatrix}";
                        std::istringstream iss2(block.substr(ipos));
                        std::string line2;
                        std::getline(iss2, line2); // skip "IF n"
                        for (int row = 0; row < n; row++) {
                            if (row) tex << " \\\\ ";
                            std::getline(iss2, line2);
                            std::istringstream rss(line2);
                            for (int col = 0; col < n; col++) {
                                if (col) tex << " & ";
                                int v; rss >> v;
                                tex << v;
                            }
                        }
                        tex << "\\end{smallmatrix}\\right)$\n\n";
                    }

                    shown++;
                }
                pos = end + 3;
            }
            if (n_entries > 500) {
                tex << "\\bigskip\\ldots " << n_entries - 500 << " more omitted.\n\n";
            }
        }

        tex << "\\end{document}\n";
        tex.close();
        std::cout << "LaTeX: " << texfile << "\n";
    }

    return 0;
}
