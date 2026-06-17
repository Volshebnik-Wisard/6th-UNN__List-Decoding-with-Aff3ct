/*!
 * \file  Decoder_FRS2.hpp
 * \brief List decoder for Folded RS — Algorithm 17.3.1 (sliding window, s <= m).
 *
 * KEY DIFFERENCE from Decoder_FRS1:
 *   Introduces parameter s  (1 <= s <= m_fold).
 *   Instead of one constraint per folded position, uses a "sliding window" of
 *   size s over each folded block, giving (m-s+1) constraints per position.
 *
 *   Step 1: Interpolate Q(X, Y1,..,Ys)  — (s+1)-variate polynomial
 *     Constraints: Q(gamma^(im+j), y_{im+j}, .., y_{im+j+s-1}) = 0
 *                  for i in [0,N_frs), j in [0, m-s]
 *     D = floor( (N_frs*(m-s+1) - K + 1) / (s+1) )
 *
 *   Step 2: Root finding identical to FRS1 but with s instead of m.
 *     List size <= q^(s-1).
 *
 *   Error radius: s/(s+1) * (1 - m*R/(m-s+1))
 *
 *   To approach capacity 1-R-eps: choose s=Theta(1/eps), m=Theta(1/eps^2).
 */
#ifndef DECODER_FRS2_HPP
#define DECODER_FRS2_HPP

#include <vector>
#include <functional>

#include "Tools/Code/RS/RS_polynomial_generator.hpp"
#include "Module/Decoder/Decoder_SIHO.hpp"

namespace aff3ct {
namespace module {

template <typename B = int, typename R = float>
class Decoder_FRS2 : public Decoder_SIHO<B, R>
{
public:
    using S = int;

    using SelectorFn = std::function<int(
        const std::vector<std::vector<S>>& /*candidates*/,
        const std::vector<std::vector<S>>& /*y_folded*/,
        const R* /*llr_soft*/,
        int /*N_frs*/,
        int /*K_rs*/
    )>;

protected:
    // RS / GF
    const int K_rs, N_rs, m_gf;
    const int n_rdncy, n_rdncy_bits;
    const int N_p2_1;
    const std::vector<int>& alpha_to;
    const std::vector<int>& index_of;

    // FRS
    const int m_fold;   // book's m
    const int s_win;    // book's s  (sliding window size, s <= m_fold)
    const int N_frs;    // N_rs / m_fold
    const int D;        // interpolation degree
    const int list_limit;

    SelectorFn selector_fn;

    // Buffers
    std::vector<B>              YH_Nb;
    std::vector<S>              YH_N;
    std::vector<std::vector<S>> y_folded; // [N_frs][m_fold]

    // A_coeffs[j] for j=0..s_win  (Q has s_win+1 terms)
    std::vector<std::vector<int>> A_coeffs;

    // Output list: each entry = K_rs info symbols
    std::vector<std::vector<S>> candidate_list;

public:
    /*!
     * \param K          info symbols
     * \param N          codeword length (RS symbols)
     * \param GF         poly generator
     * \param m_fold     folding parameter
     * \param s_win      sliding window size (1 <= s_win <= m_fold)
     * \param list_limit 0 = unlimited (bounded by q^(s-1))
     * \param selector   list selector; nullptr => min-Hamming fallback
     */
    Decoder_FRS2(const int K, const int N,
                 const tools::RS_polynomial_generator& GF,
                 const int m_fold,
                 const int s_win,
                 const int list_limit = 0,
                 SelectorFn selector = nullptr);

    virtual ~Decoder_FRS2() = default;
    virtual Decoder_FRS2<B,R>* clone() const;

    const std::vector<std::vector<S>>& get_last_list() const { return candidate_list; }
    int get_m_fold() const { return m_fold; }
    int get_s_win()  const { return s_win;  }
    int get_N_frs()  const { return N_frs;  }

protected:
    virtual int _decode_siho(const R* Y_N, int8_t* CWD, B* V_K, const size_t frame_id);
    virtual int _decode_hiho(const B* Y_N, int8_t* CWD, B* V_K, const size_t frame_id);

private:
    int _decode_gf(const std::vector<std::vector<S>>& y,
                   const R* llr_soft,
                   S* V_K_sym,
                   const size_t frame_id);

    bool interpolate(const std::vector<std::vector<S>>& y);
    void find_roots();
    int  eval_B(int x_poly) const;

    inline int gf_mul(int a, int b) const;
    inline int gf_add(int a, int b) const { return a ^ b; }
};

} // namespace module
} // namespace aff3ct

#endif // DECODER_FRS2_HPP
