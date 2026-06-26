// gen_sugra_special_m1.cpp
// ----------------------------------------------------------------------------
// Special single-external phase1: attach a GAUGE-LESS (-1) external curve
// (intersection number 1) to every LST base curve i whose positive-integer
// kernel component X_i == 1, over the ENTIRE LST catalog. No combo / stacking.
//
// For each LST entry:
//   base_IF = reconstruct_IF(entry)
//   X       = find_positive_integer_kernel(base_IF)   // fiber/null class
//   for each curve i with X_i == 1:
//       new_IF = attach_curve(base_IF, i, -1, 1)       // gauge-less (-1), int 1
//       keep iff  sig_pos == 1  &&  NHC passes  &&  H_neutral >= 0   (our filter)
//   dedup by (graph-iso IF + cid)   // same key family as phase1
//
// Anomaly uses the canonical --use-lst-T convention: T = LST curve count
// (g_current_base_T = base_IF.rows()), external does NOT add to the anomaly T.
//
// Usage: ./gen_sugra_special_m1 <catalog> [T_max=193] [T_min=0]
//                               [--out FILE] [--cc-budget 8.0]
// ----------------------------------------------------------------------------
#include <Eigen/Dense>
#include "theory_sugra.h"
#include "sugra_generator.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>
#ifdef _OPENMP
#include <omp.h>
#endif

// --- graph + eigenvalue + cid dedup key (self-contained copy; helpers live in
//     sugra_generator.h). Identical fingerprint to phase1's dedup_key_v2. ---
static DedupKey128 dedup_key_cached(const Eigen::MatrixXi& IF,
                                    const std::vector<double>& eigenvalues,
                                    const std::string& extra) {
    int n = IF.rows();
    std::vector<uint8_t> buf;
    buf.reserve(64 + n * 32 + extra.size());
    for (int i = 0; i < n; i++) {
        double v = eigenvalues[i];
        if (std::abs(v) < 1e-10) v = 0.0;
        buf_append_i64(buf, (int64_t)std::llround(v * 1e8));
    }
    buf.push_back(0xFF);
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
            if (j != i && IF(i, j) != 0) { vds[i].deg++; vds[i].nbr_si.push_back(IF(j, j)); }
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
    std::vector<std::tuple<int, int, int>> edges;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (IF(i, j) != 0) {
                int a = IF(i, i), b = IF(j, j);
                if (a > b) std::swap(a, b);
                edges.push_back({a, b, IF(i, j)});
            }
    std::sort(edges.begin(), edges.end());
    for (auto& [a, b, w] : edges) { buf_append_i32(buf, a); buf_append_i32(buf, b); buf_append_i32(buf, w); }
    buf.push_back(0xFD);
    buf_append_bytes(buf, extra.data(), extra.size());
    return hash128(buf);
}

struct Res {
    int catalog_id;
    std::string catalog_type;
    int base_T;        // entry.T (sig_neg of LST base)
    int n_base;        // LST curve count (= anomaly T)
    int target_idx;    // curve in base where (-1) attached (X_i==1)
    int target_si;     // self-int of that target curve
    Eigen::MatrixXi final_IF;
    NHCResult nhc;
    AnomalyResult anom;
    SigInfo sig;
    std::vector<double> evs;
};

// Attach a gauge-less (-1) (int_num 1) to every curve i with kernel X_i==1,
// apply our filter (sig_pos==1 -> NHC -> H_neutral>=0), append survivors to `out`.
// ctr accumulates [tried, sig_pass, nhc_pass, anom_pass]. Shared by catalog + dummy bases.
static void process_base(const Eigen::MatrixXi& base_IF, int cid, const std::string& ctype,
                         int base_T, double cc_budget, std::vector<Res>& out, long long ctr[4]) {
    int n_base = base_IF.rows();
    Eigen::VectorXi X = find_positive_integer_kernel(base_IF);
    if (X.size() == 0) return;                 // no clean positive integer kernel
    g_current_base_T = n_base;                 // --use-lst-T: anomaly T = LST curve count
    for (int i = 0; i < n_base; i++) {
        if (X(i) != 1) continue;               // attach only on kernel-component-1 curves
        ctr[0]++;
        Eigen::MatrixXi new_IF = attach_curve(base_IF, i, -1, 1);
        std::vector<double> evs;
        SigInfo sig = compute_sig(new_IF, evs);
        if (sig.sig_pos != 1) continue;
        ctr[1]++;
        NHCResult nhc = check_nhc(new_IF, cc_budget);
        if (!nhc.passes) continue;
        ctr[2]++;
        AnomalyResult anom = compute_anomaly(new_IF, nhc);
        if (anom.H_neutral < 0) continue;
        ctr[3]++;
        Res r;
        r.catalog_id = cid; r.catalog_type = ctype;
        r.base_T = base_T; r.n_base = n_base;
        r.target_idx = i; r.target_si = base_IF(i, i);
        r.final_IF = std::move(new_IF);
        r.nhc = nhc; r.anom = anom; r.sig = sig; r.evs = std::move(evs);
        out.push_back(std::move(r));
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <catalog> [T_max=193] [T_min=0] [--out FILE] [--cc-budget 8.0]\n";
        return 1;
    }
    std::string catalog_file = argv[1];
    int T_max = (argc > 2 && argv[2][0] != '-') ? std::atoi(argv[2]) : 193;
    int T_min = (argc > 3 && argv[3][0] != '-') ? std::atoi(argv[3]) : 0;
    std::string out_file = "special_m1.cat";
    double cc_budget = 8.0;
    for (int i = 2; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--out" && i + 1 < argc) out_file = argv[++i];
        else if (a == "--cc-budget" && i + 1 < argc) cc_budget = std::atof(argv[++i]);
    }

    auto catalog = load_catalog(catalog_file);
    std::cerr << "Loaded " << catalog.size() << " LST entries; T in [" << T_min
              << "," << T_max << "], cc_budget=" << cc_budget << "\n";

    std::vector<Res> all;
    long long ctr[4] = {0, 0, 0, 0};   // [tried, sig_pass, nhc_pass, anom_pass]

    // Dummy LST work-list — A-type (nc=2..T_max+1) and D-type (nc=3..T_max+1),
    // matching gen_sugra_phase1's dummy enumeration. These are the small LST bases
    // (incl. base_T 1,2) absent from unified.cat (whose min T is 3).
    struct DW { bool isA; int nc; };
    std::vector<DW> dummies;
    for (int nc = std::max(2, T_min + 1); nc <= T_max + 1; nc++) dummies.push_back({true, nc});
    for (int nc = std::max(3, T_min + 1); nc <= T_max + 1; nc++) dummies.push_back({false, nc});

    #pragma omp parallel
    {
        std::vector<Res> local;
        long long lc[4] = {0, 0, 0, 0};
        // --- unified.cat entries ---
        #pragma omp for schedule(dynamic) nowait
        for (int ei = 0; ei < (int)catalog.size(); ei++) {
            const auto& entry = catalog[ei];
            if (entry.T < T_min || entry.T > T_max) continue;
            Eigen::MatrixXi base_IF = reconstruct_IF(entry);
            process_base(base_IF, entry.id, entry.type, entry.T, cc_budget, local, lc);
        }
        // --- dummy LSTs (A/D) ---
        #pragma omp for schedule(dynamic)
        for (int di = 0; di < (int)dummies.size(); di++) {
            Eigen::MatrixXi base_IF = dummies[di].isA ? build_dummy_A(dummies[di].nc)
                                                      : build_dummy_D(dummies[di].nc);
            SigInfo bs = compute_sig_fast(base_IF);
            if (bs.sig_neg + 1 > T_max) continue;   // matches phase1 try_dummy T gate
            int cid = dummies[di].isA ? -dummies[di].nc : -(1000 + dummies[di].nc);
            std::string ctype = (dummies[di].isA ? "DM:A(" : "DM:D(")
                                + std::to_string(dummies[di].nc) + ")";
            process_base(base_IF, cid, ctype, bs.sig_neg, cc_budget, local, lc);
        }
        #pragma omp critical
        {
            all.insert(all.end(), std::make_move_iterator(local.begin()),
                       std::make_move_iterator(local.end()));
            for (int k = 0; k < 4; k++) ctr[k] += lc[k];
        }
    }
    long long n_attach_tried = ctr[0], n_pass_sig = ctr[1], n_pass_nhc = ctr[2], n_pass_anom = ctr[3];

    // Deterministic order before dedup (parallel collection order varies run-to-run).
    std::sort(all.begin(), all.end(), [](const Res& a, const Res& b){
        if (a.catalog_id != b.catalog_id) return a.catalog_id < b.catalog_id;
        if (a.base_T != b.base_T) return a.base_T < b.base_T;
        return a.target_idx < b.target_idx;
    });

    // --- dedup by (graph-iso IF + cid) ---
    std::unordered_set<DedupKey128, DedupKey128Hash> seen;
    std::vector<Res> ded;
    ded.reserve(all.size());
    for (auto& r : all) {
        std::string extra = "specialm1|cid:" + std::to_string(r.catalog_id);
        DedupKey128 key = dedup_key_cached(r.final_IF, r.evs, extra);
        if (seen.insert(key).second) ded.push_back(std::move(r));
    }

    std::cerr << "attach tried=" << n_attach_tried
              << "  sig_pos==1=" << n_pass_sig
              << "  NHC=" << n_pass_nhc
              << "  H_neutral>=0=" << n_pass_anom
              << "  after dedup=" << ded.size() << "\n";

    // --- write .cat (phase1 format; REMAINING 0 since no combo) ---
    std::ofstream out(out_file);
    out << "# special_m1 SUGRA catalog: gauge-less (-1) external (int_num 1)\n";
    out << "# attached to LST base curves with kernel component X_i == 1\n";
    out << "# Total: " << ded.size() << " bases\n#\n";
    for (int ri = 0; ri < (int)ded.size(); ri++) {
        auto& r = ded[ri];
        int n = r.final_IF.rows();
        out << "ENTRY " << ri << "\n";
        // SPEC spec_id tag ext_si target_si int_num is_hat1
        // tag "1" = gauge-less (-1) external (mirrors "2" = gauge-less (-2))
        out << "SPEC 0 1 -1 " << r.target_si << " 1 0\n";
        out << "BASE " << r.catalog_id << " " << r.catalog_type
            << " " << r.base_T << " " << r.n_base << "\n";
        out << "ATTACH " << r.target_idx << " " << (n - 1) << "\n";
        out << "PHYSICS " << r.anom.T << " " << r.anom.H_charged
            << " " << r.anom.V << " " << r.anom.H_neutral
            << " " << r.sig.det << " " << r.sig.sig_pos
            << " " << r.sig.sig_neg << " " << r.sig.sig_zero << "\n";
        out << "IF " << n << "\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) { if (j) out << " "; out << r.final_IF(i, j); }
            out << "\n";
        }
        out << "REMAINING 0\n";
        out << "END\n\n";
    }
    out.close();
    std::cerr << "Wrote " << ded.size() << " bases to " << out_file << "\n";
    return 0;
}
