// Validate schur_sig (incremental Schur signature) against compute_sig (full eigensolve).
#include "sugra_generator.h"
#include <random>
#include <iostream>

int main() {
    std::mt19937 rng(20260617);
    std::uniform_int_distribution<int> offd(-4, 4);
    std::uniform_int_distribution<int> dg(-12, -1);
    int tests = 0, pos_fail = 0, full_fail = 0;
    for (int trial = 0; trial < 5000; trial++) {
        int n = 4 + (int)(rng() % 24);   // base size
        int k = 1 + (int)(rng() % 3);    // appended curves
        int N = n + k;
        Eigen::MatrixXi M(N, N);
        for (int i = 0; i < N; i++) {
            M(i, i) = dg(rng);
            for (int j = i + 1; j < N; j++) {
                int v = (rng() % 4 == 0) ? offd(rng) : 0;  // sparse-ish off-diagonal
                M(i, j) = v; M(j, i) = v;
            }
        }
        Eigen::MatrixXi B = M.block(0, 0, n, n);
        BaseSchur bs = precompute_base_schur(B);
        SigInfo ss = schur_sig(bs, M);
        std::vector<double> ev;
        SigInfo cs = compute_sig(M, ev);
        tests++;
        if (ss.sig_pos != cs.sig_pos) {
            pos_fail++;
            if (pos_fail <= 8)
                std::cout << "POS MISMATCH n=" << n << " k=" << k
                          << " schur_pos=" << ss.sig_pos << " full_pos=" << cs.sig_pos
                          << " (schur z=" << ss.sig_zero << ", full z=" << cs.sig_zero << ")\n";
        }
        if (ss.sig_pos != cs.sig_pos || ss.sig_neg != cs.sig_neg || ss.sig_zero != cs.sig_zero)
            full_fail++;
    }
    std::cout << "tests=" << tests
              << "  sig_pos mismatches=" << pos_fail
              << "  full-inertia mismatches=" << full_fail << "\n";
    return pos_fail ? 1 : 0;
}
