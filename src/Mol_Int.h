// Molecular Integrals over Cartesian Gaussian Functions
#include <vector>
#include <cmath>
#include "umath.h"
#include <Eigen/Dense>

namespace Elevia {

    /*** =========================== ***/
    /*** Cartesian Gaussian Function ***/
    /*** =========================== ***/

    class PGF {
    public:
        PGF() = default;
        PGF(const Eigen::Vector3d& a, const double& b,
            const Eigen::Vector3d& c, const Eigen::Vector3d& d) :
            r(a), zeta(b), n(c), R(d) {}

        Eigen::Vector3d r;
        double zeta;
        Eigen::Vector3d n;
        Eigen::Vector3d R;

        double lambda = n(0) + n(1) + n(2);
        double N = pow((2 / uovie::pi), 0.75) * pow(2.0, lambda)
            * pow(zeta, (2 * lambda + 3) / 4) * pow(uovie::factorial2(2 * std::round(n(0)) - 1)
                * uovie::factorial2(2 * std::round(n(1)) - 1) * uovie::factorial2(2 * std::round(n(2)) - 1), -0.5);
        double uvalue = pow(r(0) - R(0), n(0)) * pow(r(1) - R(1), n(1)) * pow(r(2) - R(2), n(2))
            * exp(-zeta * ((r - R).squaredNorm()));
        double nvalue = N * uvalue;
    };

    /*** Adjust n in the PGF ***/
    PGF n_adj(const PGF pgf, const int site, const double num) {
        Eigen::Vector3d c;
        for (int i = 0; i < 3; i++) {
            if (i == site)
                c(i) = num;
            else
                c(i) = 0;
        }
        Eigen::Vector3d n_new = pgf.n + c;
        PGF pgf_new(pgf.r, pgf.zeta, n_new, pgf.R);
        return pgf_new;
    }

    class Mol_Int {
    public:
        Mol_Int() = default;
        double TCO(std::vector<PGF> pgf);

    private:
        ;
    };

    /*** ============================= ***/
    /*** Three-Center Overlap Integral ***/
    /*** ============================= ***/

    double Mol_Int::TCO(std::vector<PGF> pgf) {
        double zeta_abc = pgf[0].zeta + pgf[1].zeta + pgf[2].zeta;
        Eigen::Vector3d R_abc = (pgf[0].zeta * pgf[0].R + pgf[1].zeta * pgf[1].R + pgf[2].zeta * pgf[2].R) / zeta_abc;
        double N_abc = exp(zeta_abc * (R_abc.squaredNorm()) - (pgf[0].zeta * (pgf[0].R.squaredNorm())
            + pgf[1].zeta * (pgf[1].R.squaredNorm()) + pgf[2].zeta * (pgf[2].R.squaredNorm())));
        Eigen::Matrix3d n_set;
        n_set << pgf[0].n(0), pgf[0].n(1), pgf[0].n(2),
            pgf[1].n(0), pgf[1].n(1), pgf[1].n(2),
            pgf[2].n(0), pgf[2].n(1), pgf[2].n(2);

        double int_val = pow(uovie::pi / zeta_abc, 1.5) * N_abc;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (std::round(n_set(i, j)) > 0) {
                    // (a-1^j|b|c), (a|b-1^j|c), (a|b|c-1^j)
                    std::vector<PGF> pgf_adj(3);
                    for (int k = 0; k < 3; k++) {
                        if (k == j)
                            pgf_adj[k] = n_adj(pgf[k], k, -1);
                        else
                            pgf_adj[k] = pgf[k];
                    }
                    // (a-2^j|b|c)
                    std::vector<PGF> pgf_a_minus2(3);
                    pgf_a_minus2[0] = n_adj(pgf[0], j, -2);
                    pgf_a_minus2[1] = pgf[1];
                    pgf_a_minus2[2] = pgf[2];

                    // (a|b-2^j|c)
                    std::vector<PGF> pgf_b_minus2(3);
                    pgf_b_minus2[0] = pgf[1];
                    pgf_b_minus2[1] = n_adj(pgf[1], j, -2);
                    pgf_b_minus2[2] = pgf[2];

                    // (a|b|c-2^j)
                    std::vector<PGF> pgf_c_minus2(3);
                    pgf_c_minus2[0] = pgf[1];
                    pgf_c_minus2[1] = n_adj(pgf[1], j, -2);
                    pgf_c_minus2[2] = pgf[2];

                    int_val = (R_abc(j) - pgf[i].R(j)) * TCO(pgf_adj)
                        + (pgf[0].n(j) * TCO(pgf_a_minus2) + pgf[1].n(j) * TCO(pgf_b_minus2)
                            + pgf[2].n(j) * TCO(pgf_c_minus2)) / (2 * zeta_abc);
                }
                else if (std::round(n_set(i, j)) == 0) {
                    continue;
                }
                else if (std::round(n_set(i, j)) < 0) {
                    return 0;
                }
            }
        }
        return int_val;
    }
}
