// Windows SDK pollution guard - must come before any Windows headers
#if defined(_WIN32) || defined(_WIN64)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  ifndef NOGDI
#    define NOGDI
#  endif
#  ifndef NOUSER
#    define NOUSER
#  endif
#endif

/*!
 * \file  Decoder_FRS2.cpp
 * \brief List decoder for Folded RS — Algorithm 17.3.1 (sliding window) implementation.
 *
 * Differences from Decoder_FRS1:
 *   1. Q is (s+1)-variate instead of (m+1)-variate.
 *   2. Number of interpolation constraints: N_frs * (m_fold - s_win + 1)
 *      (each folded position contributes m-s+1 sliding-window constraints).
 *   3. D = floor( (N_frs*(m-s+1) - K + 1) / (s+1) )
 *   4. Root finding: same upper-triangular structure but diagonal is B(gamma^j)
 *      of degree s-1 => at most s-1 free vars => list size <= q^(s-1).
 */
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <limits>

#include "Tools/Exception/exception.hpp"
#include "Tools/Perf/common/hard_decide.h"
#include "Tools/Algo/Bit_packer.hpp"
#include "Tools/Math/utils.h"
#include "Decoder_FRS2.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

// ─────────────────────────────────────────────────────────────────────────────
// GF helpers
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
inline int Decoder_FRS2<B,R>::gf_mul(int a, int b) const
{
    if (a == 0 || b == 0) return 0;
    const int ia = index_of[a], ib = index_of[b];
    if (ia == -1 || ib == -1) return 0;
    return alpha_to[(ia + ib) % N_p2_1];
}

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
Decoder_FRS2<B,R>::Decoder_FRS2(const int K, const int N,
                                  const tools::RS_polynomial_generator& GF,
                                  const int m_fold,
                                  const int s_win,
                                  const int list_limit,
                                  SelectorFn selector)
: Decoder_SIHO<B,R>(K * GF.get_m(), N * GF.get_m()),
  K_rs        (K                                                ),
  N_rs        (N                                                ),
  m_gf        (GF.get_m()                                       ),
  n_rdncy     (GF.get_n_rdncy()                                 ),
  n_rdncy_bits(GF.get_n_rdncy() * GF.get_m()                    ),
  N_p2_1      (tools::next_power_of_2(N) - 1                    ),
  alpha_to    (GF.get_alpha_to()                                ),
  index_of    (GF.get_index_of()                                ),
  m_fold      (m_fold                                           ),
  s_win       (s_win                                            ),
  N_frs       (N / m_fold                                       ),
  D           ( ((N/m_fold)*(m_fold - s_win + 1) - K + 1)
                / (s_win + 1)                                   ),
  list_limit  (list_limit                                       ),
  selector_fn (selector                                         ),
  YH_Nb       (N * GF.get_m()                                   ),
  YH_N        (N                                                ),
  y_folded    (N/m_fold, std::vector<S>(m_fold, 0)              ),
  A_coeffs    (s_win + 1                                        )
{
    this->set_name("Decoder_FRS2");

    if (m_fold < 1 || N % m_fold != 0)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "m_fold must be >= 1 and divide N_rs.");
    if (s_win < 1 || s_win > m_fold)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "s_win must satisfy 1 <= s_win <= m_fold.");
    if (N_rs - K_rs != n_rdncy)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "N_rs - K_rs != n_rdncy.");

    A_coeffs[0].assign(D + K_rs, 0);
    for (int j = 1; j <= s_win; ++j)
        A_coeffs[j].assign(D + 1, 0);

    // Default selector: minimum Hamming distance to received (symbol level)
    if (!selector_fn)
    {
        selector_fn = [](const std::vector<std::vector<S>>& cands,
                         const std::vector<std::vector<S>>& y_fld,
                         const R* /*llr*/, int N_f, int /*K*/) -> int
        {
            if (cands.empty()) return -1;
            if (cands.size() == 1) return 0;
            int best_idx = 0, best_dist = std::numeric_limits<int>::max();
            for (int c = 0; c < (int)cands.size(); ++c)
            {
                int dist = 0;
                for (int j = 0; j < N_f; ++j)
                    for (int i = 0; i < (int)y_fld[j].size(); ++i)
                        if (y_fld[j][i] != cands[c][j*(int)y_fld[j].size()+i])
                            ++dist;
                if (dist < best_dist) { best_dist = dist; best_idx = c; }
            }
            return best_idx;
        };
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// clone
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
Decoder_FRS2<B,R>* Decoder_FRS2<B,R>::clone() const
{
    auto d = new Decoder_FRS2(*this);
    d->deep_copy(*this);
    return d;
}

// ─────────────────────────────────────────────────────────────────────────────
// _decode_siho / _decode_hiho
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int Decoder_FRS2<B,R>::_decode_siho(const R* Y_N, int8_t* CWD,
                                     B* V_K, const size_t frame_id)
{
    tools::hard_decide(Y_N, YH_Nb.data(), this->N);
    tools::Bit_packer::pack(YH_Nb.data(), YH_N.data(), this->N, 1, false, m_gf);
    for (int j = 0; j < N_frs; ++j)
        for (int i = 0; i < m_fold; ++i)
            y_folded[j][i] = YH_N[j * m_fold + i];

    std::vector<S> V_K_sym(K_rs, 0);
    int status = _decode_gf(y_folded, Y_N, V_K_sym.data(), frame_id);
    CWD[0] = (status == 0) ? 1 : 0;
    tools::Bit_packer::unpack(V_K_sym.data(), V_K, this->K, 1, false, m_gf);
    return status;
}

template <typename B, typename R>
int Decoder_FRS2<B,R>::_decode_hiho(const B* Y_N, int8_t* CWD,
                                     B* V_K, const size_t frame_id)
{
    tools::Bit_packer::pack(Y_N, YH_N.data(), this->N, 1, false, m_gf);
    for (int j = 0; j < N_frs; ++j)
        for (int i = 0; i < m_fold; ++i)
            y_folded[j][i] = YH_N[j * m_fold + i];

    std::vector<S> V_K_sym(K_rs, 0);
    int status = _decode_gf(y_folded, nullptr, V_K_sym.data(), frame_id);
    CWD[0] = (status == 0) ? 1 : 0;
    tools::Bit_packer::unpack(V_K_sym.data(), V_K, this->K, 1, false, m_gf);
    return status;
}

// ─────────────────────────────────────────────────────────────────────────────
// _decode_gf
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int Decoder_FRS2<B,R>::_decode_gf(const std::vector<std::vector<S>>& y,
                                    const R* llr_soft,
                                    S* V_K_sym,
                                    const size_t /*frame_id*/)
{
    candidate_list.clear();

    if (!interpolate(y))
    {
        std::fill(V_K_sym, V_K_sym + K_rs, 0);
        return 1;
    }

    find_roots();

    if (candidate_list.empty())
    {
        std::fill(V_K_sym, V_K_sym + K_rs, 0);
        return 1;
    }

    int chosen = selector_fn(candidate_list, y, llr_soft, N_frs, K_rs);
    if (chosen < 0 || chosen >= (int)candidate_list.size())
    {
        std::fill(V_K_sym, V_K_sym + K_rs, 0);
        return 1;
    }

    const auto& best = candidate_list[chosen];
    for (int i = 0; i < K_rs; ++i)
        V_K_sym[i] = best[i];
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// interpolate  — Step 1 of Algorithm 17.3.1
//
// Q(X, Y1,..,Ys) = A0(X) + A1(X)*Y1 + .. + As(X)*Ys
// Constraints (eq 17.5):
//   Q(gamma^(im+j), y_{im+j}, .., y_{im+j+s-1}) = 0
//   for i in [0, N_frs), j in [0, m_fold - s_win]
//
// Number of constraints = N_frs * (m_fold - s_win + 1)
// Number of variables   = (D+K) + s*(D+1) = (s+1)(D+1) + K - 1
// D chosen so variables > constraints.
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
bool Decoder_FRS2<B,R>::interpolate(const std::vector<std::vector<S>>& y)
{
    const int n_constraints_per_pos = m_fold - s_win + 1;
    const int n_eqs  = N_frs * n_constraints_per_pos;
    const int n_vars = (D + K_rs) + s_win * (D + 1);

    std::vector<std::vector<int>> M(n_eqs, std::vector<int>(n_vars, 0));

    int eq = 0;
    for (int i = 0; i < N_frs; ++i)
    {
        for (int jj = 0; jj <= m_fold - s_win; ++jj, ++eq)
        {
            // Evaluation point: X = gamma^(i*m_fold + jj)
            const int x_idx = (static_cast<long long>(i) * m_fold + jj) % N_p2_1;

            // A0(X) at X = gamma^(im+jj)
            for (int d = 0; d < D + K_rs; ++d)
            {
                int val;
                if (d == 0) val = 1;
                else        val = alpha_to[(static_cast<long long>(x_idx) * d) % N_p2_1];
                M[eq][d] = val;
            }

            // Aj(X) * y_{im+jj+j-1} for j = 1..s_win
            for (int j = 1; j <= s_win; ++j)
            {
                const int col_base = (D + K_rs) + (j-1) * (D + 1);
                // y index: flat RS position = i*m_fold + jj + (j-1)
                const int flat_idx = i * m_fold + jj + (j - 1);
                const int yval = (flat_idx < N_rs) ? y[flat_idx / m_fold][flat_idx % m_fold] : 0;

                for (int d = 0; d <= D; ++d)
                {
                    int xd;
                    if (d == 0) xd = 1;
                    else        xd = alpha_to[(static_cast<long long>(x_idx) * d) % N_p2_1];
                    M[eq][col_base + d] = gf_mul(xd, yval);
                }
            }
        }
    }

    // Gaussian elimination (same as in Decoder_FRS1)
    std::vector<int> pivot_col(n_eqs, -1);
    int current_row = 0;

    for (int col = 0; col < n_vars && current_row < n_eqs; ++col)
    {
        int pivot_row = -1;
        for (int row = current_row; row < n_eqs; ++row)
            if (M[row][col] != 0) { pivot_row = row; break; }
        if (pivot_row == -1) continue;

        std::swap(M[current_row], M[pivot_row]);
        pivot_col[current_row] = col;

        const int inv_p = alpha_to[(N_p2_1 - index_of[M[current_row][col]]) % N_p2_1];
        for (int c = col; c < n_vars; ++c)
            M[current_row][c] = gf_mul(M[current_row][c], inv_p);

        for (int row = 0; row < n_eqs; ++row)
        {
            if (row == current_row || M[row][col] == 0) continue;
            const int factor = M[row][col];
            for (int c = col; c < n_vars; ++c)
                M[row][c] = gf_add(M[row][c], gf_mul(factor, M[current_row][c]));
        }
        ++current_row;
    }

    const int rank = current_row;
    std::vector<bool> is_pivot(n_vars, false);
    for (int r = 0; r < rank; ++r)
        if (pivot_col[r] >= 0) is_pivot[pivot_col[r]] = true;

    int free_var = -1;
    for (int col = n_vars - 1; col >= 0; --col)
        if (!is_pivot[col]) { free_var = col; break; }
    if (free_var == -1) return false;

    std::vector<int> sol(n_vars, 0);
    sol[free_var] = 1;
    for (int r = rank - 1; r >= 0; --r)
    {
        int pc = pivot_col[r];
        if (pc < 0) continue;
        int val = 0;
        for (int c = pc + 1; c < n_vars; ++c)
            val = gf_add(val, gf_mul(M[r][c], sol[c]));
        sol[pc] = val;
    }

    for (int d = 0; d < D + K_rs; ++d)
        A_coeffs[0][d] = sol[d];
    for (int j = 1; j <= s_win; ++j)
    {
        const int base = (D + K_rs) + (j-1) * (D + 1);
        for (int d = 0; d <= D; ++d)
            A_coeffs[j][d] = sol[base + d];
    }

    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// eval_B  — B(x) = a_{1,0} + a_{2,0}*x + .. + a_{s,0}*x^{s-1}
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int Decoder_FRS2<B,R>::eval_B(int x_poly) const
{
    int result = 0;
    int x_pow  = 1;
    for (int l = 1; l <= s_win; ++l)
    {
        result = gf_add(result, gf_mul(A_coeffs[l][0], x_pow));
        x_pow  = gf_mul(x_pow, x_poly);
    }
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// find_roots  — enumerate affine subspace
//   Identical logic to Decoder_FRS1::find_roots but uses s_win instead of m_fold
//   => list size <= q^(s_win - 1)
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
void Decoder_FRS2<B,R>::find_roots()
{
    candidate_list.clear();

    // Factor out common X power
    int ell = std::numeric_limits<int>::max();
    for (int j = 0; j <= s_win; ++j)
    {
        const auto& coef = A_coeffs[j];
        for (int d = 0; d < (int)coef.size(); ++d)
            if (coef[d] != 0) { ell = std::min(ell, d); break; }
    }
    if (ell == std::numeric_limits<int>::max()) return;
    if (ell > 0)
    {
        for (int j = 0; j <= s_win; ++j)
        {
            auto& coef = A_coeffs[j];
            if (ell < (int)coef.size())
                coef.erase(coef.begin(), coef.begin() + ell);
            else
                coef.assign(1, 0);
        }
    }

    // Determine free variables: B(gamma^j) == 0 => f_j is free
    std::vector<int>  B_at_gammaj(K_rs);
    std::vector<bool> is_free(K_rs, false);
    for (int j = 0; j < K_rs; ++j)
    {
        int gj_poly = (j == 0) ? 1 : alpha_to[j % N_p2_1];
        B_at_gammaj[j] = eval_B(gj_poly);
        if (B_at_gammaj[j] == 0) is_free[j] = true;
    }

    std::vector<int> free_vars;
    for (int j = 0; j < K_rs; ++j)
        if (is_free[j]) free_vars.push_back(j);

    const int n_free = (int)free_vars.size(); // <= s_win - 1
    const int q = 1 << m_gf;

    long long total = 1;
    for (int i = 0; i < n_free; ++i)
    {
        total *= q;
        if (list_limit > 0 && total > list_limit) { total = list_limit; break; }
    }
    if (total == 0) total = 1; // at least try the zero assignment

    std::vector<int> free_vals(n_free, 0);

    auto enumerate_one = [&]() -> std::vector<S>
    {
        std::vector<S> f(K_rs, 0);
        for (int fi = 0; fi < n_free; ++fi)
            f[free_vars[fi]] = free_vals[fi];

        for (int j = 0; j < K_rs; ++j)
        {
            if (is_free[j]) continue;

            int a0j = (j < (int)A_coeffs[0].size()) ? A_coeffs[0][j] : 0;
            int rhs = a0j;

            for (int l = 0; l < j; ++l)
            {
                // b^(j)_l = sum_{iota=1}^{s} A_coeffs[iota][j-l] * gamma^(l*(iota-1))
                int bjl = 0;
                for (int iota = 1; iota <= s_win; ++iota)
                {
                    const int coef_idx = j - l;
                    int aval = (coef_idx < (int)A_coeffs[iota].size())
                               ? A_coeffs[iota][coef_idx] : 0;
                    if (aval == 0) continue;
                    int exp  = (static_cast<long long>(l) * (iota - 1)) % N_p2_1;
                    int gpow = (exp == 0) ? 1 : alpha_to[exp];
                    bjl = gf_add(bjl, gf_mul(aval, gpow));
                }
                rhs = gf_add(rhs, gf_mul(f[l], bjl));
            }

            const int Bj_idx = index_of[B_at_gammaj[j]];
            const int inv_Bj = alpha_to[(N_p2_1 - Bj_idx) % N_p2_1];
            f[j] = gf_mul(rhs, inv_Bj);
        }
        return f;
    };

    for (long long cnt = 0; cnt < total; ++cnt)
    {
        candidate_list.push_back(enumerate_one());

        for (int fi = n_free - 1; fi >= 0; --fi)
        {
            free_vals[fi]++;
            if (free_vals[fi] < q) break;
            free_vals[fi] = 0;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Explicit instantiation
// ─────────────────────────────────────────────────────────────────────────────
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Decoder_FRS2<B_8,  Q_8 >;
template class aff3ct::module::Decoder_FRS2<B_16, Q_16>;
template class aff3ct::module::Decoder_FRS2<B_32, Q_32>;
template class aff3ct::module::Decoder_FRS2<B_64, Q_64>;
#else
template class aff3ct::module::Decoder_FRS2<B, Q>;
#endif
