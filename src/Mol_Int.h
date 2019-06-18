/*
 *  Molecular Integrals over Cartesian Gaussian Functions
 */

// standard C++ headers
#include <vector>
#include <cmath>
#include <chrono>

// Elevia headers
#include "Global.h"
#include "umath.h"

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Libint Gaussian integrals library
#include <libint2.hpp>

using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

extern Elevia::fio vie;

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

    /*** ========================================== ***/
    /*** Three-Center Overlap Integral (unfinished) ***/
    /*** ========================================== ***/

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

    /*** =================================== ***/
    /*** Many-body Integrals (invoke libint) ***/
    /*** =================================== ***/

    // function declarations
    size_t nbasis(const std::vector<libint2::Shell>& shells);
    size_t max_nprim(const std::vector<libint2::Shell>& shells);
    int max_l(const std::vector<libint2::Shell>& shells);
    std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
    Matrix one_body_ints(const std::vector<libint2::Shell>& shells, libint2::Operator obtype,
        const std::vector<Elevia::atom>& atoms);
    Matrix two_body_fock(const std::vector<libint2::Shell>& shells, const Matrix& D);

    // number of basis functions
    size_t nbasis(const std::vector<libint2::Shell>& shells) {
        size_t n = 0;
        for (const auto& shell : shells)
            n += shell.size();
        return n;
    }

    // max number of primitives in shells
    size_t max_nprim(const std::vector<libint2::Shell>& shells) {
        size_t n = 0;
        for (auto shell : shells)
            n = std::max(shell.nprim(), n);
        return n;
    }

    // max angular momentum of shells
    int max_l(const std::vector<libint2::Shell>& shells) {
        int l = 0;
        for (auto shell : shells)
            for (auto c : shell.contr)
                l = std::max(c.l, l);
        return l;
    }

    // map shell to basis function
    std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
        std::vector<size_t> result;
        result.reserve(shells.size());

        size_t n = 0;
        for (auto shell : shells) {
            result.push_back(n);
            n += shell.size();
        }

        return result;
    }

    // computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
    // in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
    Matrix compute_soad(const std::vector<Elevia::atom>& atoms) {

        // number of atomic orbitals
        size_t nao = 0;
        for (const auto& atom : atoms) {
            const auto Z = atom.atom_num;
            if (Z == 1 || Z == 2) // H, He
                nao += 1;
            else if (Z <= 10) // Li - Ne
                nao += 5;
            else
                throw "SOAD with Z > 10 is not yet supported";
        }

        // compute the minimal basis density
        Matrix D(nao, nao);
        size_t ao_offset = 0; // first AO of this atom
        for (const auto& atom : atoms) {
            const auto Z = atom.atom_num;
            if (Z == 1 || Z == 2) { // H, He
                D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
                ao_offset += 1;
            }
            else if (Z <= 10) {
                D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
                D(ao_offset + 1, ao_offset + 1) = (Z == 3) ? 1 : 2; // Li? only 1 electron in 2s, else 2 electrons
                // smear the remaining electrons in 2p orbitals
                const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4) / 3 : 0;
                for (auto xyz = 0; xyz != 3; ++xyz)
                    D(ao_offset + 2 + xyz, ao_offset + 2 + xyz) = num_electrons_per_2p;
                ao_offset += 5;
            }
        }

        return D * 0.5; // we use densities normalized to # of electrons/2
    }

    // one-body Integrals
    Matrix one_body_ints(const std::vector<libint2::Shell> & shells,
        libint2::Operator obtype,
        const std::vector<Elevia::atom> & atoms)
    {
        using libint2::Shell;
        using libint2::Engine;
        using libint2::Operator;

        const auto n = nbasis(shells);
        Matrix result(n, n);

        // construct the overlap integrals engine
        Engine engine(obtype, max_nprim(shells), max_l(shells), 0);

        if (obtype == Operator::nuclear) {
            std::vector<std::pair<real_t, std::array<real_t, 3>>> q;
            for (const auto& atom : atoms) {
                q.push_back({ static_cast<real_t>(atom.atom_num), {{atom.R[0], atom.R[1], atom.R[2]}} });
            }
            engine.set_params(q);
        }

        auto shell2bf = map_shell_to_basis_function(shells);

        // buf[0] points to the target shell set after every call  to engine.compute()
        const auto& buf = engine.results();

        // loop over unique shell pairs, {s1,s2} such that s1 >= s2
        // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
        for (auto s1 = 0; s1 != shells.size(); ++s1) {

            auto bf1 = shell2bf[s1]; // first basis function in this shell
            auto n1 = shells[s1].size();

            for (auto s2 = 0; s2 <= s1; ++s2) {

                auto bf2 = shell2bf[s2];
                auto n2 = shells[s2].size();

                // compute shell pair
                engine.compute(shells[s1], shells[s2]);

                // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
                Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
                result.block(bf1, bf2, n1, n2) = buf_mat;
                if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
                    result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

            }
        }

        return result;
    }


    Matrix two_body_fock(const std::vector<libint2::Shell>& shells,const Matrix& D)
    {
        using libint2::Shell;
        using libint2::Engine;
        using libint2::Operator;

        auto time_elapsed = std::chrono::duration<double>::zero();

        const auto n = nbasis(shells);
        Matrix G(n, n);

        // construct the 2-electron repulsion integrals engine
        Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

        auto shell2bf = map_shell_to_basis_function(shells);

        const auto& buf = engine.results();

        /*** take permutational symmetries into account ***/
        // loop over permutationally-unique set of shells
        for (auto s1 = 0; s1 != shells.size(); ++s1) {

            auto bf1_first = shell2bf[s1]; // first basis function in this shell
            auto n1 = shells[s1].size();   // number of basis functions in this shell

            for (auto s2 = 0; s2 <= s1; ++s2) {

                auto bf2_first = shell2bf[s2];
                auto n2 = shells[s2].size();

                for (auto s3 = 0; s3 <= s1; ++s3) {

                    auto bf3_first = shell2bf[s3];
                    auto n3 = shells[s3].size();

                    const auto s4_max = (s1 == s3) ? s2 : s3;
                    for (auto s4 = 0; s4 <= s4_max; ++s4) {

                        auto bf4_first = shell2bf[s4];
                        auto n4 = shells[s4].size();

                        // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
                        auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
                        auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
                        auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
                        auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

                        const auto tstart = std::chrono::high_resolution_clock::now();

                        engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                        const auto* buf_1234 = buf[0];
                        if (buf_1234 == nullptr)
                            continue; // if all integrals screened out, skip to next quartet

                        const auto tstop = std::chrono::high_resolution_clock::now();
                        time_elapsed += tstop - tstart;

                        // ANSWER
                        // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
                        //    F(a,b) += (ab|cd) * D(c,d)
                        //    F(c,d) += (ab|cd) * D(a,b)
                        //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
                        //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
                        //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
                        //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
                        // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
                        //    i.e. the number of the integrals/sets equivalent to it
                        // 3) the end result must be symmetrized
                        for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                            const auto bf1 = f1 + bf1_first;
                            for (auto f2 = 0; f2 != n2; ++f2) {
                                const auto bf2 = f2 + bf2_first;
                                for (auto f3 = 0; f3 != n3; ++f3) {
                                    const auto bf3 = f3 + bf3_first;
                                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                                        const auto bf4 = f4 + bf4_first;

                                        const auto value = buf_1234[f1234];

                                        const auto value_scal_by_deg = value * s1234_deg;

                                        G(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                                        G(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                                        G(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                                        G(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                                        G(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                                        G(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                                    }
                                }
                            }
                        }

                    }
                }
            }
        }

        // symmetrize the result and return
        Matrix Gt = G.transpose();
        return 0.5 * (G + Gt);
    }
}
