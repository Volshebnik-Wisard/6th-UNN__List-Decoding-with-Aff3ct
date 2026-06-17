/*!
 * \file Decoder_LDPC_BP_list.hxx
 * Implementation of Decoder_LDPC_BP_list<B,R>.
 * Included by Decoder_LDPC_BP_list.hpp — do NOT compile separately.
 */

#ifndef DECODER_LDPC_BP_LIST_HXX_
#define DECODER_LDPC_BP_LIST_HXX_

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "Tools/Code/LDPC/Syndrome/LDPC_syndrome.hpp"
#include "Tools/Perf/common/hard_decide.h"

#include "Decoder_LDPC_BP_list.hpp"

namespace aff3ct
{
namespace module
{

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
Decoder_LDPC_BP_list<B,R>
::Decoder_LDPC_BP_list(const int                    K,
                       const int                    N,
                       const int                    n_ite,
                       const tools::Sparse_matrix  &H,
                       const std::vector<uint32_t> &info_bits_pos,
                       const int                    list_size,
                       const bool                   enable_syndrome,
                       const int                    syndrome_depth)
: Decoder_LDPC_BP_flooding_SPA<B,R>(K, N, n_ite, H, info_bits_pos,
                                    enable_syndrome, syndrome_depth),
  list_size  (list_size),
  n_ambiguous(list_size >= 2
              ? std::min(static_cast<int>(
                    std::floor(std::log2(static_cast<double>(list_size)))),
                    K)
              : 0),
  hard_ref (static_cast<size_t>(N), static_cast<B>(0)),
  candidate(static_cast<size_t>(N), static_cast<B>(0))
{
    this->set_name("Decoder_LDPC_BP_list");
}

// ─────────────────────────────────────────────────────────────────────────────
// clone
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
Decoder_LDPC_BP_list<B,R>*
Decoder_LDPC_BP_list<B,R>::clone() const
{
    auto m = new Decoder_LDPC_BP_list(*this);
    m->deep_copy(*this);
    return m;
}

// ─────────────────────────────────────────────────────────────────────────────
// _list_search
//   Reads this->post[] (filled by parent _decode()), enumerates 2^n_ambiguous
//   candidates over the least-reliable information bits, picks the best
//   syndrome-valid one. Writes best codeword into best_cw[0..N-1].
//   Returns true if at least one candidate satisfies the syndrome.
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
bool
Decoder_LDPC_BP_list<B,R>::_list_search(std::vector<B> &best_cw)
{
    const int N = this->N;
    const int K = this->K;

    // 1. Hard-decide the full posterior vector
    //    post[i] >= 0  →  bit = 0  (positive LLR = more likely 0)
    //    post[i] <  0  →  bit = 1
    for (int i = 0; i < N; i++)
        hard_ref[i] = (this->post[i] < static_cast<R>(0))
                      ? static_cast<B>(1) : static_cast<B>(0);

    // 2. Sort info-bit positions by ascending |post| (least reliable first)
    //    Only need to find the n_ambiguous smallest — use partial_sort.
    struct Entry { int cw_pos; R abs_llr; };
    std::vector<Entry> order;
    order.reserve(static_cast<size_t>(K));
    for (int k = 0; k < K; k++)
    {
        const int pos = static_cast<int>(this->info_bits_pos[k]);
        order.push_back({ pos, std::abs(this->post[pos]) });
    }

    if (n_ambiguous > 0)
    {
        const int n = (n_ambiguous < K) ? n_ambiguous : K;
        std::partial_sort(order.begin(), order.begin() + n, order.end(),
                          [](const Entry &a, const Entry &b)
                          { return a.abs_llr < b.abs_llr; });
    }

    // Collect the codeword positions that will be enumerated
    std::vector<int> flip_pos(static_cast<size_t>(n_ambiguous));
    for (int i = 0; i < n_ambiguous; i++)
        flip_pos[i] = order[i].cw_pos;

    // 3. Enumerate all 2^n_ambiguous candidates
    const int n_cand = 1 << n_ambiguous;  // = 1 when n_ambiguous == 0

    bool found_valid = false;
    R    best_cost   = std::numeric_limits<R>::infinity();

    for (int pattern = 0; pattern < n_cand; pattern++)
    {
        // Start from the BP hard-decision codeword
        std::copy(hard_ref.begin(), hard_ref.end(), candidate.begin());

        // Flip bits selected by this pattern; accumulate "cost"
        R cost = static_cast<R>(0);
        for (int b = 0; b < n_ambiguous; b++)
        {
            if ((pattern >> b) & 1)
            {
                candidate[flip_pos[b]] ^= static_cast<B>(1);
                cost += std::abs(this->post[flip_pos[b]]);
            }
        }

        // Syndrome check
        const bool valid =
            tools::LDPC_syndrome::check_hard(candidate.data(), this->H);

        if (valid)
        {
            // Among valid candidates keep the cheapest (minimum flipping cost)
            if (!found_valid || cost < best_cost)
            {
                best_cost   = cost;
                found_valid = true;
                std::copy(candidate.begin(), candidate.end(), best_cw.begin());
            }
        }
        else if (!found_valid && cost < best_cost)
        {
            // No valid candidate yet — keep the cheapest invalid one as fallback
            best_cost = cost;
            std::copy(candidate.begin(), candidate.end(), best_cw.begin());
        }
    }

    return found_valid;
}

// ─────────────────────────────────────────────────────────────────────────────
// _decode_siho  — returns K info bits in V_K
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int
Decoder_LDPC_BP_list<B,R>::_decode_siho(const R     *Y_N,
                                         int8_t      *CWD,
                                         B           *V_K,
                                         const size_t frame_id)
{
    // Run full BP (populates this->post[])
    const int bp_status = this->_decode(Y_N, frame_id);

    // List search → best codeword
    std::vector<B> best_cw(static_cast<size_t>(this->N));
    const bool valid = _list_search(best_cw);

    // Extract K information bits
    for (int k = 0; k < this->K; k++)
        V_K[k] = best_cw[static_cast<size_t>(this->info_bits_pos[k])];

    // CWD[0] == 0  ↔  decoder believes codeword is valid
    CWD[0] = valid ? static_cast<int8_t>(0) : static_cast<int8_t>(1);

    return valid ? 0 : bp_status;
}

// ─────────────────────────────────────────────────────────────────────────────
// _decode_siho_cw  — returns full codeword (N bits) in V_N
// ─────────────────────────────────────────────────────────────────────────────
template <typename B, typename R>
int
Decoder_LDPC_BP_list<B,R>::_decode_siho_cw(const R     *Y_N,
                                             int8_t      *CWD,
                                             B           *V_N,
                                             const size_t frame_id)
{
    const int bp_status = this->_decode(Y_N, frame_id);

    std::vector<B> best_cw(static_cast<size_t>(this->N));
    const bool valid = _list_search(best_cw);

    std::copy(best_cw.begin(), best_cw.end(), V_N);

    CWD[0] = valid ? static_cast<int8_t>(0) : static_cast<int8_t>(1);
    return valid ? 0 : bp_status;
}

} // namespace module
} // namespace aff3ct

#endif /* DECODER_LDPC_BP_LIST_HXX_ */
