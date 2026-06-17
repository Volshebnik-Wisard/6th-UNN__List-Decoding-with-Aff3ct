// Windows SDK pollution guard
#if defined(_WIN32) || defined(_WIN64)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#endif

/*!
 * \file  Decoder_FRS1.cpp
 * \brief List decoder for Folded RS — Algorithm 17.2.1 (clean rewrite)
 *
 * RS-КОД: evaluation code.
 *   f(X) = f_0 + f_1*X + .. + f_{K-1}*X^{K-1}  — инф. полином степени < K
 *   Кодовое слово: cw[j] = f(α^j), j=0..N-1,  α=alpha_to[1]
 *   Свёрнутое слово: y[i][r] = cw[i*m + r]
 *
 * СОГЛАСИЕ: f согласуется с y в позиции i если
 *   f(α^{im+r}) = y[i][r]  для всех r=0..m-1
 *
 * ШАГ 1 — ИНТЕРПОЛЯЦИЯ:
 *   Ищем Q(X,Y1,..,Ym) = A0(X) + A1(X)*Y1 + .. + Am(X)*Ym
 *   deg(A0) ≤ D+K-1,  deg(Aj) ≤ D
 *   Условие: Q(α^{im}, y[i][0], .., y[i][m-1]) = 0  для i=0..N_frs-1
 *   D = floor((N_frs - K + 1) / (m + 1))
 *   Число переменных (D+K) + m*(D+1) > N_frs → ненулевое Q существует.
 *   Решаем методом Гаусса над GF.
 *
 * ШАГ 2 — НАХОЖДЕНИЕ КОРНЕЙ:
 *   Подставляем f(X), f(αX), .., f(α^{m-1}X) в Q(X,...):
 *     R(X) = A0(X) + A1(X)*f(X) + A2(X)*f(αX) + .. + Am(X)*f(α^{m-1}X)
 *   R ≡ 0 тогда и только тогда когда f — корень Q.
 *   Раскрываем R(X) по степеням X и приравниваем коэффициенты 0..K-1 к нулю.
 *   Это даёт верхнетреугольную систему на f_0..f_{K-1} с диагональю B(α^j).
 *
 * ВЫХОДНЫЕ ДАННЫЕ:
 *   V_K_sym[i] = f(α^{n_rdncy + i}), i=0..K-1
 *   (информационные символы = значения f в точках α^{n_rdncy}..α^{N-1})
 */
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>

#include "Tools/Exception/exception.hpp"
#include "Tools/Perf/common/hard_decide.h"
#include "Tools/Algo/Bit_packer.hpp"
#include "Tools/Math/utils.h"
#include "Decoder_FRS1.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

// ─── GF arithmetic ───────────────────────────────────────────────────────────

template <typename B, typename R>
inline int Decoder_FRS1<B,R>::gf_mul(int a, int b) const
{
    if (a == 0 || b == 0) return 0;
    const int ia = index_of[a], ib = index_of[b];
    if (ia == -1 || ib == -1) return 0;
    return alpha_to[(ia + ib) % N_p2_1];
}

template <typename B, typename R>
inline int Decoder_FRS1<B,R>::gf_pow(int base_idx, int exp) const
{
    // base_idx — показатель (index form), возвращает α^(base_idx*exp) в poly form
    if (exp == 0) return 1;
    if (base_idx == -1) return 0; // 0^exp = 0
    return alpha_to[(static_cast<long long>(base_idx) * exp) % N_p2_1];
}

// Схема Горнера: p[0] + p[1]*x + p[2]*x^2 + ...
template <typename B, typename R>
int Decoder_FRS1<B,R>::eval_poly(const std::vector<int>& p, int x_poly) const
{
    if (p.empty()) return 0;
    // Горнер: result = p[n] * x^{n} + .. + p[0]
    // = (...((p[n-1]*x + p[n-2])*x + p[n-3])*x + ...) * x + p[0]  — НЕВЕРНО
    // Правильно: p(x) = p[0] + p[1]*x + ... = p[0] + x*(p[1] + x*(p[2]+...))
    // Горнер снизу вверх (от старшего к младшему, потом умножаем на x):
    int result = p.back(); // p[deg]
    for (int d = (int)p.size() - 2; d >= 0; --d)
        result = gf_add(gf_mul(result, x_poly), p[d]);
    return result;
}

// ─── Constructor ─────────────────────────────────────────────────────────────

template <typename B, typename R>
Decoder_FRS1<B,R>::Decoder_FRS1(const int K, const int N,
                                  const tools::RS_polynomial_generator& GF,
                                  const int m_fold,
                                  const int list_limit,
                                  SelectorFn selector)
: Decoder_SIHO<B,R>(K * GF.get_m(), N * GF.get_m()),
  K_rs        (K),
  N_rs        (N),
  m_gf        (GF.get_m()),
  n_rdncy     (GF.get_n_rdncy()),
  n_rdncy_bits(GF.get_n_rdncy() * GF.get_m()),
  N_p2_1      (tools::next_power_of_2(N) - 1),
  alpha_to    (GF.get_alpha_to()),
  index_of    (GF.get_index_of()),
  g_poly      (GF.get_g()),
  m_fold      (m_fold),
  N_frs       (N / m_fold),
  D           ((N/m_fold - K + 1) / (m_fold + 1)),
  list_limit  (list_limit),
  selector_fn (selector),
  YH_Nb       (N * GF.get_m()),
  YH_N        (N),
  y_folded    (N/m_fold, std::vector<S>(m_fold, 0)),
  A_coeffs    (m_fold + 1)
{
    this->set_name("Decoder_FRS1");

    if (m_fold < 1 || N % m_fold != 0)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "m_fold must be >= 1 and divide N_rs.");
    if (N_rs - K_rs != n_rdncy)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "N_rs - K_rs != n_rdncy.");
    if (K_rs >= N_frs)
    {
        std::stringstream ss;
        ss << "K_rs(" << K_rs << ") >= N_frs(" << N_frs
           << "): need K_sym < N_sym/m_fold.";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, ss.str());
    }
    if (D <= 0)
    {
        std::stringstream ss;
        ss << "D=(N_frs-K+1)/(m+1)=(" << N_frs << "-" << K_rs
           << "+1)/(" << m_fold << "+1)=" << D << "<=0.";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, ss.str());
    }

    A_coeffs[0].assign(D + K_rs, 0);   // A0: deg <= D+K-1
    for (int j = 1; j <= m_fold; ++j)
        A_coeffs[j].assign(D + 1, 0);  // Aj: deg <= D

    if (!selector_fn)
    {
        // cands[c] = коэффициенты f[0..K-1].
        // Для выбора лучшего кандидата считаем Хэмминг-расстояние
        // между f(α^{j*m+i}) и y[j][i] для каждой позиции.
        const int mf = m_fold, Nf = N_frs, Np = N_p2_1;
        const std::vector<int>* at_ptr = &alpha_to;
        selector_fn = [mf, Nf, Np, at_ptr, this](
                                const std::vector<std::vector<S>>& cands,
                                const std::vector<std::vector<S>>& y_fld,
                                const R*, int, int) -> int
        {
            if (cands.empty()) return -1;
            int best = 0, best_d = std::numeric_limits<int>::max();
            for (int c = 0; c < (int)cands.size(); ++c)
            {
                int dist = 0;
                for (int j = 0; j < Nf; ++j)
                    for (int i = 0; i < mf; ++i)
                    {
                        const int pos = j*mf + i;
                        const int xp  = (pos == 0) ? 1 : (*at_ptr)[pos % Np];
                        if (y_fld[j][i] != eval_poly(cands[c], xp)) ++dist;
                    }
                if (dist < best_d) { best_d = dist; best = c; }
            }
            return best;
        };
    }
}

template <typename B, typename R>
Decoder_FRS1<B,R>* Decoder_FRS1<B,R>::clone() const
{
    auto d = new Decoder_FRS1(*this);
    d->deep_copy(*this);
    return d;
}

// ─── decode_siho / decode_hiho ───────────────────────────────────────────────

template <typename B, typename R>
int Decoder_FRS1<B,R>::_decode_siho(const R* Y_N, int8_t* CWD,
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
int Decoder_FRS1<B,R>::_decode_hiho(const B* Y_N, int8_t* CWD,
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

// ─── _decode_gf ──────────────────────────────────────────────────────────────

template <typename B, typename R>
int Decoder_FRS1<B,R>::_decode_gf(const std::vector<std::vector<S>>& y,
                                    const R* llr_soft,
                                    S* V_K_sym,
                                    const size_t frame_id)
{
    static int dbg = 0;
    const bool do_dbg = (dbg < 2);

    candidate_list.clear();

    if (!interpolate(y))
    {
        if (do_dbg) { std::cout << "[DBG] interpolate FAILED\n"; ++dbg; }
        std::fill(V_K_sym, V_K_sym + K_rs, 0);
        return 1;
    }

    if (do_dbg)
    {
        std::cout << "[DBG frame=" << dbg << "] K=" << K_rs
                  << " N=" << N_rs << " N_frs=" << N_frs
                  << " m=" << m_fold << " D=" << D
                  << " n_rdncy=" << n_rdncy << "\n";
        std::cout << "  A0[0..4]=";
        for (int d=0;d<std::min(5,D+K_rs);++d) std::cout<<A_coeffs[0][d]<<" ";
        std::cout << "\n  A1[0..4]=";
        for (int d=0;d<std::min(5,D+1);++d) std::cout<<A_coeffs[1][d]<<" ";
        std::cout << "\n";

        // Q-constraint check
        int qok=0, qbad=0;
        for (int i=0;i<N_frs;++i)
        {
            const int xi_idx=(int)((static_cast<long long>(m_fold)*i)%N_p2_1);
            const int xi = (xi_idx==0)?1:alpha_to[xi_idx];
            int Qv=eval_poly(A_coeffs[0],xi);
            for (int j=1;j<=m_fold;++j)
                Qv=gf_add(Qv,gf_mul(eval_poly(A_coeffs[j],xi),y[i][j-1]));
            if (Qv==0) ++qok; else ++qbad;
        }
        std::cout << "  Q-constraint: " << qok << "/" << N_frs
                  << " OK, " << qbad << " violations\n";
    }

    find_roots(y);

    if (do_dbg)
    {
        std::cout << "  candidates=" << candidate_list.size() << "\n";
        if (!candidate_list.empty())
        {
            const auto& f0 = candidate_list[0];
            std::cout << "  f[0..5]=";
            for(int i=0;i<std::min(6,K_rs);++i) std::cout<<f0[i]<<" ";
            std::cout<<"\n";

            // Проверка: f(α^0), f(α^1), f(α^2) vs y[0][0], y[1][0], y[2][0]
            for(int pos=0;pos<3;++pos){
                const int xp=(pos==0)?1:alpha_to[pos%N_p2_1];
                int fval=eval_poly(f0,xp);
                int yval=y[pos/m_fold][pos%m_fold];
                std::cout<<"  f(a^"<<pos<<")="<<fval<<" y="<<yval
                         <<(fval==yval?" MATCH":" DIFF")<<"\n";
            }

            // Проверка R(X)=A0(X)+A1(X)*f(X) — первые 5 коэфф должны быть 0
            // Коэфф X^n: A0[n] + sum_{r=0}^{n} A1[n-r]*f[r]
            std::cout<<"  R[0..4]=";
            for(int n=0;n<std::min(5,D+K_rs);++n){
                int rn=(n<(int)A_coeffs[0].size())?A_coeffs[0][n]:0;
                for(int r=0;r<=n&&r<K_rs;++r){
                    int nr=n-r;
                    int a1=(nr<(int)A_coeffs[1].size())?A_coeffs[1][nr]:0;
                    int fr=(r<(int)f0.size())?f0[r]:0;
                    rn=gf_add(rn,gf_mul(a1,fr));
                }
                std::cout<<rn<<" ";
            }
            std::cout<<"\n";

            // agreements через eval_poly
            int agr=0;
            for (int j=0;j<N_frs;++j)
                for (int i=0;i<m_fold;++i)
                {
                    const int pos=j*m_fold+i;
                    const int xp=(pos==0)?1:alpha_to[pos%N_p2_1];
                    if (y[j][i]==eval_poly(f0,xp)) ++agr;
                }
            std::cout << "  cand[0] agreements=" << agr << "/" << N_rs << "\n";
        }
        ++dbg;
    }

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

    // Информационные символы — коэффициенты f[0..K-1].
    // Кодер: packed_U_K[i] = коэффициент f_i, они пишутся в X_N напрямую.
    // Декодер должен вернуть те же коэффициенты.
    const auto& best = candidate_list[chosen];
    for (int i = 0; i < K_rs; ++i)
        V_K_sym[i] = best[i];

    return 0;
}

// ─── interpolate ─────────────────────────────────────────────────────────────
// Строим линейную систему над GF и находим ненулевой вектор из ядра матрицы.
// Строка i: Q(α^{m*i}, y[i][0], .., y[i][m-1]) = 0
//   => sum_{d=0}^{D+K-1} a0_d * (α^{mi})^d
//    + sum_{j=1}^{m} y[i][j-1] * sum_{d=0}^{D} a{j}_d * (α^{mi})^d = 0
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
bool Decoder_FRS1<B,R>::interpolate(const std::vector<std::vector<S>>& y)
{
    // Сбросить коэффициенты
    for (int j = 0; j <= m_fold; ++j)
        std::fill(A_coeffs[j].begin(), A_coeffs[j].end(), 0);

    const int n_vars = (D + K_rs) + m_fold * (D + 1);
    const int n_eqs  = N_frs;

    // Матрица [n_eqs × n_vars]
    std::vector<std::vector<int>> M(n_eqs, std::vector<int>(n_vars, 0));

    for (int i = 0; i < n_eqs; ++i)
    {
        // xi = α^{m*i}  (в poly form)
        const long long exp_mi = (static_cast<long long>(m_fold) * i) % N_p2_1;
        const int xi = (exp_mi == 0) ? 1 : alpha_to[(int)exp_mi];

        // Переменные A0: столбцы 0..D+K-1
        // Вклад a0_d: (α^{mi})^d = xi^d
        int xi_pow = 1; // xi^0
        for (int d = 0; d < D + K_rs; ++d)
        {
            M[i][d] = xi_pow;
            xi_pow = gf_mul(xi_pow, xi);
        }

        // Переменные Aj (j=1..m): столбцы (D+K) + (j-1)*(D+1) .. + D
        // Вклад a{j}_d: y[i][j-1] * (α^{mi})^d
        for (int j = 1; j <= m_fold; ++j)
        {
            const int col_base = (D + K_rs) + (j-1)*(D+1);
            const int yij = y[i][j-1];
            int xi_pow2 = 1;
            for (int d = 0; d <= D; ++d)
            {
                M[i][col_base + d] = gf_mul(yij, xi_pow2);
                xi_pow2 = gf_mul(xi_pow2, xi);
            }
        }
    }

    // Гаусс: приведение к row echelon form, поиск вектора из ядра
    std::vector<int> pivot_col(n_eqs, -1);
    int cur_row = 0;

    for (int col = 0; col < n_vars && cur_row < n_eqs; ++col)
    {
        int prow = -1;
        for (int row = cur_row; row < n_eqs; ++row)
            if (M[row][col] != 0) { prow = row; break; }
        if (prow < 0) continue;

        std::swap(M[cur_row], M[prow]);
        pivot_col[cur_row] = col;

        const int pv  = M[cur_row][col];
        const int inv = alpha_to[(N_p2_1 - index_of[pv]) % N_p2_1];
        for (int c = col; c < n_vars; ++c)
            M[cur_row][c] = gf_mul(M[cur_row][c], inv);

        for (int row = 0; row < n_eqs; ++row)
        {
            if (row == cur_row || M[row][col] == 0) continue;
            const int f = M[row][col];
            for (int c = col; c < n_vars; ++c)
                M[row][c] = gf_add(M[row][c], gf_mul(f, M[cur_row][c]));
        }
        ++cur_row;
    }

    const int rank = cur_row;

    // Найти свободную переменную
    std::vector<bool> is_pivot(n_vars, false);
    for (int r = 0; r < rank; ++r)
        if (pivot_col[r] >= 0) is_pivot[pivot_col[r]] = true;

    // ВАЖНО: выбираем ПЕРВУЮ свободную переменную (наименьший индекс).
    // Если выбрать последнюю — получим A со старшим свободным коэффициентом,
    // тогда a_{j,0}=0 для всех j, B(x)=0, и find_roots не работает.
    // Первая свободная переменная гарантирует ненулевой младший коэффициент
    // хотя бы у одного A_j после факторизации X^ell.
    int free_var = -1;
    for (int col = 0; col < n_vars; ++col)
        if (!is_pivot[col]) { free_var = col; break; }
    if (free_var < 0) return false;

    // Обратная подстановка: free_var = 1
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

    // Распаковка в A_coeffs
    for (int d = 0; d < D + K_rs; ++d)
        A_coeffs[0][d] = sol[d];
    for (int j = 1; j <= m_fold; ++j)
    {
        const int base = (D + K_rs) + (j-1)*(D+1);
        for (int d = 0; d <= D; ++d)
            A_coeffs[j][d] = sol[base + d];
    }
    return true;
}

// ─── eval_B ──────────────────────────────────────────────────────────────────
// B(x) = a_{1,0} + a_{2,0}*x + a_{3,0}*x^2 + .. + a_{m,0}*x^{m-1}
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int Decoder_FRS1<B,R>::eval_B(int x_poly) const
{
    int result = 0, xpow = 1;
    for (int l = 1; l <= m_fold; ++l)
    {
        result = gf_add(result, gf_mul(A_coeffs[l][0], xpow));
        xpow = gf_mul(xpow, x_poly);
    }
    return result;
}

// ─── find_roots ──────────────────────────────────────────────────────────────
// Раскладываем R(X) = A0(X) + sum_{j=1}^{m} Aj(X)*f(α^{j-1}*X) по степеням.
// Коэффициент при X^n в Aj(X)*f(α^{j-1}*X):
//   = sum_{d+r=n} a_{j,d} * f_r * (α^{j-1})^r
// Условие c_n = 0 для n=0..K-1 даёт верхнетреугольную систему на f_0..f_{K-1}.
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
void Decoder_FRS1<B,R>::find_roots(const std::vector<std::vector<S>>& y_ref)
{
    candidate_list.clear();

    // Убрать общий множитель X^ell
    int ell = std::numeric_limits<int>::max();
    for (int j = 0; j <= m_fold; ++j)
        for (int d = 0; d < (int)A_coeffs[j].size(); ++d)
            if (A_coeffs[j][d] != 0) { ell = std::min(ell, d); break; }
    if (ell == std::numeric_limits<int>::max()) return;
    if (ell > 0)
        for (int j = 0; j <= m_fold; ++j)
        {
            auto& c = A_coeffs[j];
            if (ell < (int)c.size())
                c.erase(c.begin(), c.begin() + ell);
            else
                c.assign(1, 0); // всё нулевое — заменяем на [0]
        }

    // Проверяем что хотя бы у одного A_j (j>=1) константный член ненулевой
    // (иначе B(x)=0 и find_roots вернёт пустой список)
    bool b_nonzero = false;
    for (int j = 1; j <= m_fold; ++j)
        if (!A_coeffs[j].empty() && A_coeffs[j][0] != 0) { b_nonzero = true; break; }
    if (!b_nonzero)
    {
        // Редкий случай: Q = X^ell * Q', B ≡ 0.
        // Пробуем взять другое решение из пространства ядра.
        // Для простоты: возвращаем пустой список (нет корней).
        return;
    }

    // B(α^j) = a_{1,0} + a_{2,0}*α^j + .. + a_{m,0}*α^{j*(m-1)}
    // Если B(α^j)=0 → f_j свободна
    std::vector<int>  Bval(K_rs);
    std::vector<bool> is_free(K_rs, false);
    for (int j = 0; j < K_rs; ++j)
    {
        // α^j: index form = j % N_p2_1
        const int gj_poly = (j == 0) ? 1 : alpha_to[j % N_p2_1];
        Bval[j] = eval_B(gj_poly);
        if (Bval[j] == 0) is_free[j] = true;
    }

    std::vector<int> free_idx;
    for (int j = 0; j < K_rs; ++j)
        if (is_free[j]) free_idx.push_back(j);

    const int n_free = (int)free_idx.size();
    const int q      = 1 << m_gf;

    long long total = 1;
    for (int i = 0; i < n_free; ++i)
    {
        if (list_limit > 0 && total > (long long)list_limit) break;
        total *= q;
    }
    if (list_limit > 0 && total > (long long)list_limit) total = list_limit;
    if (n_free == 0) total = 1;

    std::vector<int> free_vals(n_free, 0);

    for (long long cnt = 0; cnt < total; ++cnt)
    {
        std::vector<S> f(K_rs, 0);
        for (int fi = 0; fi < n_free; ++fi)
            f[free_idx[fi]] = free_vals[fi];

        bool valid = true;
        for (int n = 0; n < K_rs && valid; ++n)
        {
            if (is_free[n]) continue;

            // c_n = 0:
            // a_{0,n} + f_n * B(α^n) + sum_{r=0}^{n-1} f_r * b_r^(n) = 0
            // b_r^(n) = sum_{j=1}^{m} a_{j, n-r} * (α^{j-1})^r
            //         = sum_{j=1}^{m} A_coeffs[j][n-r] * α^{r*(j-1)}

            int a0n = (n < (int)A_coeffs[0].size()) ? A_coeffs[0][n] : 0;
            int rhs = a0n; // GF(2^m): -x = x

            for (int r = 0; r < n; ++r)
            {
                // b_r^(n)
                int brn = 0;
                for (int j = 1; j <= m_fold; ++j)
                {
                    const int nr = n - r; // index into Aj
                    if (nr >= (int)A_coeffs[j].size()) continue;
                    const int ajnr = A_coeffs[j][nr];
                    if (ajnr == 0) continue;
                    // α^{r*(j-1)} in poly form
                    const long long exp_val = (static_cast<long long>(r)*(j-1)) % N_p2_1;
                    const int gpow = (exp_val == 0) ? 1 : alpha_to[(int)exp_val];
                    brn = gf_add(brn, gf_mul(ajnr, gpow));
                }
                rhs = gf_add(rhs, gf_mul(f[r], brn));
            }

            if (index_of[Bval[n]] == -1) { valid = false; break; }
            const int inv_B = alpha_to[(N_p2_1 - index_of[Bval[n]]) % N_p2_1];
            f[n] = gf_mul(rhs, inv_B);
        }

        if (valid)
        {
            // Храним КОЭФФИЦИЕНТЫ f[0..K-1] как кандидата.
            // Кодовое слово при необходимости вычисляется в selector.
            candidate_list.push_back(f);
        }

        // Счётчик base-q
        for (int fi = n_free-1; fi >= 0; --fi)
        {
            if (++free_vals[fi] < q) break;
            free_vals[fi] = 0;
        }
    }
}

// ─── Explicit instantiation ───────────────────────────────────────────────────
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Decoder_FRS1<B_8,  Q_8 >;
template class aff3ct::module::Decoder_FRS1<B_16, Q_16>;
template class aff3ct::module::Decoder_FRS1<B_32, Q_32>;
template class aff3ct::module::Decoder_FRS1<B_64, Q_64>;
#else
template class aff3ct::module::Decoder_FRS1<B, Q>;
#endif
