/*!
 * \file   Encoder_FRS.hpp
 * \brief  Folded Reed-Solomon encoder (evaluation code).
 */
#ifndef ENCODER_FRS_HPP
#define ENCODER_FRS_HPP

#include <vector>
#include "Tools/Code/RS/RS_polynomial_generator.hpp"
#include "Module/Encoder/Encoder.hpp"

namespace aff3ct {
namespace module {

template <typename B = int>
class Encoder_FRS : public Encoder<B>
{
public:
    using S = int;

protected:
    const int K_rs, N_rs, m_gf;
    const int n_rdncy, n_rdncy_bits;

    const std::vector<int>& alpha_to;
    const std::vector<int>& index_of;
    const std::vector<int>& g;

    const int m_fold;
    const int N_frs;

    std::vector<S> packed_U_K;
    std::vector<S> bb;
    std::vector<S> rs_cw;

public:
    std::vector<std::vector<S>> cw_folded;

    Encoder_FRS(const int K, const int N,
                const tools::RS_polynomial_generator& GF,
                const int m_fold);

    virtual ~Encoder_FRS() = default;
    virtual Encoder_FRS<B>* clone() const;

    const std::vector<std::vector<S>>& get_cw_folded() const { return cw_folded; }
    int get_m_fold() const { return m_fold; }
    int get_N_frs()  const { return N_frs;  }
    int get_N_rs()   const { return N_rs;   }
    int get_K_rs()   const { return K_rs;   }
    int get_m_gf()   const { return m_gf;   }

    bool is_codeword(const B* X_N);

protected:
    virtual void _encode(const B* U_K, B* X_N, const size_t frame_id);

    // Evaluate polynomial p at point x_poly using Horner scheme
    // p[0] = constant term, p[1] = coeff of x, etc.
    S eval_poly(const std::vector<S>& p, S x_poly) const;
};

} // namespace module
} // namespace aff3ct

#endif // ENCODER_FRS_HPP
