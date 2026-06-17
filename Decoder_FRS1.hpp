/*!
 * \file  Decoder_FRS1.hpp
 * \brief List decoder for Folded RS — Algorithm 17.2.1
 */
#ifndef DECODER_FRS1_HPP
#define DECODER_FRS1_HPP

#include <vector>
#include <functional>
#include <memory>

#include "Tools/Code/RS/RS_polynomial_generator.hpp"
#include "Module/Decoder/Decoder_SIHO.hpp"

namespace aff3ct {
namespace module {

template <typename B = int, typename R = float>
class Decoder_FRS1 : public Decoder_SIHO<B, R>
{
public:
    using S = int; // символ GF

    // Сигнатура селектора:
    //   candidates  — список кандидатов; каждый = N_rs значений f(γ^0)..f(γ^{N-1})
    //   y_folded    — принятое свёрнутое слово [N_frs][m_fold]
    //   llr_soft    — мягкие LLR (N_bits), может быть nullptr
    //   N_frs, K_rs — размеры
    // Возвращает индекс выбранного кандидата или -1
    using SelectorFn = std::function<int(
        const std::vector<std::vector<S>>& /*candidates (N_rs symbols each)*/,
        const std::vector<std::vector<S>>& /*y_folded [N_frs][m_fold]*/,
        const R*                           /*llr_soft*/,
        int                                /*N_frs*/,
        int                                /*K_rs*/
    )>;

protected:
    // RS / GF параметры
    const int K_rs, N_rs, m_gf;
    const int n_rdncy, n_rdncy_bits;
    const int N_p2_1;
    const std::vector<int>& alpha_to;
    const std::vector<int>& index_of;
    const std::vector<int>& g_poly;   // generator polynomial coefficients (index form)

    // FRS параметры
    const int m_fold;
    const int N_frs;
    const int D;           // параметр интерполяции = (N_frs-K+1)/(m+1)
    const int list_limit;  // 0 = без ограничений

    SelectorFn selector_fn;

    // Буферы
    std::vector<B>              YH_Nb;
    std::vector<S>              YH_N;
    std::vector<std::vector<S>> y_folded;   // [N_frs][m_fold]

    // Коэффициенты A_j(X): A_coeffs[j][d] = коэффициент степени d
    std::vector<std::vector<int>> A_coeffs; // [0..m_fold][0..D+K или D]

    // Список кандидатов: каждый — N_rs символов f(γ^0)..f(γ^{N-1})
    std::vector<std::vector<S>> candidate_list;

public:
    Decoder_FRS1(const int K, const int N,
                 const tools::RS_polynomial_generator& GF,
                 const int m_fold,
                 const int list_limit = 0,
                 SelectorFn selector  = nullptr);

    virtual ~Decoder_FRS1() = default;
    virtual Decoder_FRS1<B,R>* clone() const;

    const std::vector<std::vector<S>>& get_last_list() const { return candidate_list; }
    int get_m_fold() const { return m_fold; }
    int get_N_frs()  const { return N_frs;  }
    int get_D()      const { return D;       }

protected:
    virtual int _decode_siho(const R* Y_N, int8_t* CWD, B* V_K, const size_t frame_id);
    virtual int _decode_hiho(const B* Y_N, int8_t* CWD, B* V_K, const size_t frame_id);

private:
    int  _decode_gf(const std::vector<std::vector<S>>& y,
                    const R* llr_soft, S* V_K_sym, const size_t frame_id);
    bool interpolate(const std::vector<std::vector<S>>& y);
    void find_roots(const std::vector<std::vector<S>>& y_ref);
    int  eval_B(int x_poly) const;
    int  eval_poly(const std::vector<int>& p, int x_poly) const;

    inline int gf_mul(int a, int b) const;
    inline int gf_add(int a, int b) const { return a ^ b; }
    inline int gf_pow(int base_idx, int exp) const;
};

} // namespace module
} // namespace aff3ct

#endif // DECODER_FRS1_HPP
