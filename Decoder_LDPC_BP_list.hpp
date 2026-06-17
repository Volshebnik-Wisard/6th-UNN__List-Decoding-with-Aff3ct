/*!
 * \file Decoder_LDPC_BP_list.hpp
 *
 * Post-hoc List-BP LDPC decoder.
 *
 * Inherits from Decoder_LDPC_BP_flooding_SPA<B,R> and overrides only the
 * _decode_siho / _decode_siho_cw paths to add list expansion after the BP run.
 *
 * Algorithm (matches ldpc_code.go :: ListDecode):
 *   1. Run standard BP-flooding/SPA for n_ite iterations (parent _decode()).
 *   2. Collect posterior LLRs (this->post[]) over all N bits.
 *   3. Among the K information bits, pick the L = floor(log2(list_size))
 *      bits whose |post| is smallest  ("least certain" decisions).
 *   4. Enumerate all 2^L bit-flip patterns over those L positions.
 *   5. For each candidate codeword check the LDPC syndrome (H·c = 0 mod 2).
 *      Return the cheapest (minimum flipping cost) valid candidate.
 *      If none is valid, return the cheapest candidate overall.
 *
 * Encoder pairing (from LDPC_old.cpp):
 *   Encoder_LDPC_from_H  (G_method="IDENTITY", default), systematic.
 *   info_bits_pos = {0, 1, …, K-1}.
 *   K=144, N=271, H loaded from .alist file.
 *
 * Usage (see LDPC_list_main.cpp for a full example):
 *   module::Decoder_LDPC_BP_list<int,float> dec(
 *       K, N, n_ite, H, info_bits_pos,
 *       list_size,        // e.g. 8  (actual list = 2^floor(log2(list_size)))
 *       enable_syndrome,  // true
 *       syndrome_depth);  // 1
 *
 *   dec.decode_siho(LLR_vector, dec_bits_vector);
 */

#ifndef DECODER_LDPC_BP_LIST_HPP_
#define DECODER_LDPC_BP_LIST_HPP_

#include <cstdint>
#include <vector>

#include "Tools/Algo/Matrix/Sparse_matrix/Sparse_matrix.hpp"
#include "Module/Decoder/LDPC/BP/Flooding/SPA/Decoder_LDPC_BP_flooding_SPA.hpp"

namespace aff3ct
{
namespace module
{

template <typename B = int, typename R = float>
class Decoder_LDPC_BP_list : public Decoder_LDPC_BP_flooding_SPA<B,R>
{
private:
    int            list_size;    ///< as given by the user
    int            n_ambiguous;  ///< floor(log2(list_size)), clamped to K

    std::vector<B> hard_ref;     ///< hard-decision of post[] after BP  (len N)
    std::vector<B> candidate;    ///< working candidate codeword         (len N)

public:
    /*!
     * \param K               Number of information bits.
     * \param N               Codeword length.
     * \param n_ite           Max BP iterations.
     * \param H               Parity-check matrix (sparse).
     * \param info_bits_pos   Positions of the K info bits inside the codeword.
     * \param list_size       Desired list size (actual = 2^floor(log2(list_size))).
     * \param enable_syndrome Early-stop on valid syndrome check.
     * \param syndrome_depth  Consecutive valid checks before stopping.
     */
    Decoder_LDPC_BP_list(const int                    K,
                         const int                    N,
                         const int                    n_ite,
                         const tools::Sparse_matrix  &H,
                         const std::vector<uint32_t> &info_bits_pos,
                         const int                    list_size       = 8,
                         const bool                   enable_syndrome = true,
                         const int                    syndrome_depth  = 1);

    virtual ~Decoder_LDPC_BP_list() = default;

    virtual Decoder_LDPC_BP_list<B,R>* clone() const override;

protected:
    int _decode_siho   (const R *Y_N, int8_t *CWD, B *V_K,
                         const size_t frame_id) override;

    int _decode_siho_cw(const R *Y_N, int8_t *CWD, B *V_N,
                         const size_t frame_id) override;

private:
    /*!
     * Core list-search: reads this->post[], enumerates candidates, returns
     * the best valid codeword in best_cw. Returns true if syndrome was
     * satisfied by at least one candidate.
     */
    bool _list_search(std::vector<B> &best_cw);
};

} // namespace module
} // namespace aff3ct

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Decoder_LDPC_BP_list.hxx"
#endif

#endif /* DECODER_LDPC_BP_LIST_HPP_ */
