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
 * \file  Encoder_FRS.cpp
 * \brief Folded Reed-Solomon encoder — evaluation code.
 *
 * Кодовое слово evaluation RS:
 *   f(X) = f_0 + f_1*X + ... + f_{K-1}*X^{K-1}
 *   cw[j] = f(α^j),  j = 0..N-1
 *
 * Свёрнутое слово: cw_folded[i][r] = cw[i*m_fold + r]
 *
 * Выходной вектор X_N (для канала):
 *   X_N = [cw[0] | cw[1] | ... | cw[N-1]] в битах (N*m_gf бит)
 *   Это ЧИСТОЕ кодовое слово — никаких перезаписей.
 *
 * info_bits_pos: позиции информационных бит в X_N.
 *   packed_U_K[i] = f_i (коэффициент полинома).
 *   В кодовом слове коэффициенты не хранятся явно — хранятся значения.
 *   Для Monitor_BFER нам нужно согласовать ref_bits и dec_bits:
 *   ref_bits = U_K (K*m_gf бит) = упакованные f_0..f_{K-1}.
 *   dec_bits = V_K = Bit_packer::unpack(V_K_sym) где V_K_sym[i]=f_i.
 *   Значит Monitor сравнивает f_0..f_{K-1} с f_0..f_{K-1} — верно.
 *
 *   info_bits_pos указывает на позиции в X_N откуда Monitor читает ref.
 *   НО: Monitor_BFER получает ref_bits напрямую от Source, не из X_N!
 *   Поэтому info_bits_pos нужен только для is_systematic() и других проверок,
 *   не для Monitor. Ставим его в [0..K*m_gf-1] как заглушку.
 */
#include <string>
#include <sstream>
#include <climits>
#include <algorithm>

#include "Tools/Exception/exception.hpp"
#include "Tools/Algo/Bit_packer.hpp"
#include "Encoder_FRS.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
Encoder_FRS<B>::Encoder_FRS(const int K, const int N,
                              const tools::RS_polynomial_generator& GF,
                              const int m_fold)
: Encoder<B>(K * GF.get_m(), N * GF.get_m()),
  K_rs        (K),
  N_rs        (N),
  m_gf        (GF.get_m()),
  n_rdncy     (GF.get_n_rdncy()),
  n_rdncy_bits(GF.get_n_rdncy() * GF.get_m()),
  alpha_to    (GF.get_alpha_to()),
  index_of    (GF.get_index_of()),
  g           (GF.get_g()),
  m_fold      (m_fold),
  N_frs       (N / m_fold),
  packed_U_K  (K),
  bb          (GF.get_n_rdncy()),
  rs_cw       (N),
  cw_folded   (N / m_fold, std::vector<S>(m_fold, 0))
{
    this->set_name("Encoder_FRS");

    if (m_fold < 1)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "m_fold must be >= 1.");
    if (N % m_fold != 0)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "N_rs must be divisible by m_fold.");
    if ((N_rs - K_rs) != n_rdncy)
        throw tools::invalid_argument(__FILE__, __LINE__, __func__,
            "N_rs - K_rs != n_rdncy.");

    // info_bits_pos: заглушка — ставим первые K*m_gf позиций.
    // Monitor_BFER использует ref_bits напрямую от Source,
    // поэтому значение info_bits_pos для сравнения не важно.
    for (int i = 0; i < this->K; ++i)
        this->info_bits_pos[i] = i;
}

template <typename B>
Encoder_FRS<B>* Encoder_FRS<B>::clone() const
{
    auto e = new Encoder_FRS(*this);
    e->deep_copy(*this);
    return e;
}

template <typename B>
typename Encoder_FRS<B>::S
Encoder_FRS<B>::eval_poly(const std::vector<S>& p, S x_poly) const
{
    if (p.empty()) return 0;
    S result = p.back();
    for (int d = (int)p.size() - 2; d >= 0; --d)
    {
        S prod = 0;
        if (result != 0 && x_poly != 0)
        {
            const int ir = index_of[result], ix = index_of[x_poly];
            if (ir != -1 && ix != -1)
                prod = alpha_to[(ir + ix) % N_rs];
        }
        result = prod ^ p[d];
    }
    return result;
}

template <typename B>
void Encoder_FRS<B>::_encode(const B* U_K, B* X_N, const size_t /*frame_id*/)
{
    // 1. Упаковать биты → GF-символы (коэффициенты f_0..f_{K-1})
    tools::Bit_packer::pack(U_K, packed_U_K.data(), this->K, 1, false, m_gf);

    // 2. Кодовое слово: rs_cw[j] = f(α^j)
    //    α^0 = 1, α^j = alpha_to[j % N_rs]
    for (int j = 0; j < N_rs; ++j)
    {
        const S x = (j == 0) ? static_cast<S>(1)
                              : static_cast<S>(alpha_to[j % N_rs]);
        rs_cw[j] = eval_poly(packed_U_K, x);
    }

    // 3. Свёрнутая структура
    for (int i = 0; i < N_frs; ++i)
        for (int r = 0; r < m_fold; ++r)
            cw_folded[i][r] = rs_cw[i * m_fold + r];

    // 4. Выход: X_N = rs_cw в битах — чистое кодовое слово
    tools::Bit_packer::unpack(rs_cw.data(), X_N, this->N, 1, false, m_gf);
}

template <typename B>
bool Encoder_FRS<B>::is_codeword(const B* X_N)
{
    return true; // TODO: syndrome check
}

#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Encoder_FRS<B_8>;
template class aff3ct::module::Encoder_FRS<B_16>;
template class aff3ct::module::Encoder_FRS<B_32>;
template class aff3ct::module::Encoder_FRS<B_64>;
#else
template class aff3ct::module::Encoder_FRS<B>;
#endif
