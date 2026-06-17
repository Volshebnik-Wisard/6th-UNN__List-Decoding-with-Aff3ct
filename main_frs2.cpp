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
 * \file  main_frs2.cpp
 * \brief BER/FER simulation: Folded RS + FRS Decoder 2 (Algorithm 17.3.1).
 *
 * PARAMETERS:
 *   K_sym   — information symbols (RS)
 *   N_sym   — codeword length in symbols (RS)
 *   m_fold  — folding parameter (must divide N_sym)
 *   s_win   — sliding window size (1 <= s_win <= m_fold)
 *   list_sz — list size limit (0 = unbounded, bounded by q^(s_win-1))
 *
 * DECODING RADIUS (Theorem 17.3.3):
 *   delta = s/(s+1) * (1 - m*R/(m-s+1))
 *
 * CAPACITY-APPROACHING:
 *   Set s = Theta(1/eps), m = Theta(1/eps^2)  =>  delta >= 1 - R - eps
 */
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include <aff3ct.hpp>
using namespace aff3ct;

#include "Encoder_FRS.hpp"
#include "Decoder_FRS2.hpp"
#include "FRS_List_Selector.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════
// ТАБЛИЦА ПАРАМЕТРОВ  (N_sym=255, GF(2^8), m_gf=8)
//
// Условие корректности Decoder 2:  D = (N_frs*(m-s+1)-K+1)/(s+1) > 0
// Формула границы (Теорема 17.3.3): delta = s/(s+1) * (1 - m*R/(m-s+1))
// Приближение к ёмкости: s=Theta(1/eps), m=Theta(1/eps^2) => delta >= 1-R-eps
//
// ── Пара 1: K_sym=63, N_sym=255, R=0.247 ────────────────────────────────
//   RS baseline:  K=63, N=255, t=96
//
//   FRS Decoder 2:
//     m_fold=15, s_win=5 → N_frs=17, D=(17*11-63+1)/6=17  ✓
//       delta = 5/6*(1 - 15*0.247/11) = 0.833*0.663 = 0.552
//     m_fold=15, s_win=3 → N_frs=17, D=(17*13-63+1)/4=40  ✓
//       delta = 3/4*(1 - 15*0.247/13) = 0.75*0.715 = 0.536
//   Итого для Пары 1: K_sym=63, N_sym=255, m_fold=15, s_win=5
//
// ── Пара 2: K_sym=31, N_sym=255, R=0.122 ────────────────────────────────
//   RS baseline:  K=31, N=255, t=112
//
//   FRS Decoder 2:
//     m_fold=15, s_win=6 → N_frs=17, D=(17*10-31+1)/7=20  ✓
//       delta = 6/7*(1 - 15*0.122/10) = 0.857*0.817 = 0.700
//     m_fold=15, s_win=5 → N_frs=17, D=(17*11-31+1)/6=28  ✓
//       delta = 5/6*(1 - 15*0.122/11) = 0.833*0.833 = 0.694
//   Итого для Пары 2: K_sym=31, N_sym=255, m_fold=15, s_win=6
//
// ── Запуск из командной строки ───────────────────────────────────────────
//   cd E:\1-ННГУ\6-ВКР\2-Практика\1_RS_old
//   cmake --build build --config Release --target frs2
//   .\build\Release\frs2.exe 63 255 15 5 0 llr   <- Пара 1
//   .\build\Release\frs2.exe 31 255 15 6 0 llr   <- Пара 2
// ═══════════════════════════════════════════════════════════════════════════

struct params
{
    // ── RS параметры ──────────────────────────────────────────────────────
    // Пара 1 (умеренная скорость): K_sym=63,  N_sym=255, m_fold=15, s_win=5
    // Пара 2 (низкая скорость):   K_sym=31,  N_sym=255, m_fold=15, s_win=6
    int K_sym   = 63;   // информационных символов RS
    int N_sym   = 255;  // длина кодового слова (= 2^8 - 1 для GF(2^8))

    // ── FRS параметры ─────────────────────────────────────────────────────
    int m_fold  = 15;   // параметр складывания; 255=3*5*17, делители: 3,5,15,17,51,85
    int s_win   = 5;    // размер скользящего окна (1 <= s_win <= m_fold)
                        // чем меньше s/m, тем ближе к границе ёмкости
    int list_sz = 0;    // максимальный размер списка (0 = без ограничений,
                        // теоретический максимум q^(s_win-1) = 256^4 — большой!)

    // ── Параметры симуляции ───────────────────────────────────────────────
    int   fe        = 100;
    int   seed      = 0;
    float ebn0_min  = 0.0f;
    float ebn0_max  = 10.0f;
    float ebn0_step = 0.5f;

    // ── Стратегия выбора из списка ────────────────────────────────────────
    // "llr"     — максимум мягкой метрики LLR (лучше при AWGN)
    // "hamming" — минимум расстояния Хэмминга (запасной вариант)
    std::string selector = "llr";

    // ── Производные (заполняются в main) ─────────────────────────────────
    float R;
    int   t;
    int   m_gf;
    int   K_bits, N_bits;
};

// ─────────────────────────────────────────────────────────────────────────────
void print_params(const params& p)
{
    std::cout << "#----------------------------------------------------------\n";
    std::cout << "# FRS simulation — Decoder 2 (Algorithm 17.3.1)\n";
    std::cout << "#----------------------------------------------------------\n";
    std::cout << "#  RS symbols    K = " << p.K_sym << "\n";
    std::cout << "#  RS symbols    N = " << p.N_sym << "\n";
    std::cout << "#  Correction    t = " << p.t     << "\n";
    std::cout << "#  GF extension  m_gf = " << p.m_gf  << "\n";
    std::cout << "#  Folding   m_fold = " << p.m_fold << "\n";
    std::cout << "#  Window    s_win  = " << p.s_win  << "\n";
    std::cout << "#  FRS block N_frs  = " << p.N_sym/p.m_fold << "\n";
    std::cout << "#  List limit    L  = " << (p.list_sz ? std::to_string(p.list_sz) : "unlimited") << "\n";
    std::cout << "#  Selector       = " << p.selector << "\n";
    std::cout << "#  Code rate     R  = " << p.R << "\n";

    // Theoretical decoding radius (Theorem 17.3.3)
    const double m = p.m_fold, s = p.s_win, R = p.R;
    const double delta = (s/(s+1.0)) * (1.0 - m * R / (m - s + 1.0));
    std::cout << "#  Theory delta (Thm 17.3.3) = " << delta << "\n";

    // Capacity bound
    const double eps_equiv = 1.0 - R - delta;
    std::cout << "#  Capacity gap  eps ~ " << eps_equiv << "\n";
    std::cout << "#\n";
}

// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char** argv)
{
    params p;

    // ./main_frs2 K N m_fold s_win [list_sz] [selector]
    if (argc >= 5)
    {
        p.K_sym  = std::stoi(argv[1]);
        p.N_sym  = std::stoi(argv[2]);
        p.m_fold = std::stoi(argv[3]);
        p.s_win  = std::stoi(argv[4]);
    }
    if (argc >= 6) p.list_sz  = std::stoi(argv[5]);
    if (argc >= 7) p.selector = argv[6];

    // Validate
    if (p.N_sym % p.m_fold != 0)
    {
        std::cerr << "ERROR: N_sym (" << p.N_sym << ") not divisible by m_fold ("
                  << p.m_fold << ").\n";
        return 1;
    }
    if (p.s_win < 1 || p.s_win > p.m_fold)
    {
        std::cerr << "ERROR: s_win must be in [1, m_fold].\n";
        return 1;
    }
    if ((p.N_sym - p.K_sym) % 2 != 0)
    {
        std::cerr << "ERROR: N_sym - K_sym must be even.\n";
        return 1;
    }
    p.t = (p.N_sym - p.K_sym) / 2;

    const int N_p2 = tools::next_power_of_2(p.N_sym) - 1;
    tools::RS_polynomial_generator GF(N_p2, p.t);

    p.m_gf   = GF.get_m();
    p.K_bits = p.K_sym * p.m_gf;
    p.N_bits = p.N_sym * p.m_gf;
    p.R      = static_cast<float>(p.K_sym) / static_cast<float>(p.N_sym);

    print_params(p);

    // ── Selector ─────────────────────────────────────────────────────────
    using Sel = FRS_List_Selector<int, float>;
    Sel::SelectorFn selector_fn;
    if (p.selector == "llr")
    {
        const int n_rdncy_bits = GF.get_n_rdncy() * p.m_gf;
        selector_fn = Sel::make_max_llr(p.m_gf, n_rdncy_bits, p.K_bits);
    }
    else
    {
        selector_fn = Sel::make_min_hamming(p.m_fold);
    }

    // ── Modules ──────────────────────────────────────────────────────────
    auto source  = std::make_unique<module::Source_random<>>(p.K_bits);
    auto encoder = std::make_unique<module::Encoder_FRS<>>(p.K_sym, p.N_sym, GF, p.m_fold);
    auto decoder = std::make_unique<module::Decoder_FRS2<>>(
        p.K_sym, p.N_sym, GF, p.m_fold, p.s_win, p.list_sz, selector_fn);
    auto modem   = std::make_unique<module::Modem_BPSK<>>(p.N_bits);
    auto channel = std::make_unique<module::Channel_AWGN_LLR<>>(p.N_bits);
    auto monitor = std::make_unique<module::Monitor_BFER<>>(p.K_bits, p.fe);

    channel->set_seed(p.seed);

    // ── Reporters ────────────────────────────────────────────────────────
    auto noise = std::make_unique<tools::Sigma<>>();
    std::vector<std::unique_ptr<tools::Reporter>> reporters;
    reporters.push_back(std::make_unique<tools::Reporter_noise<>>(*noise));
    reporters.push_back(std::make_unique<tools::Reporter_BFER<>>(*monitor));
    reporters.push_back(std::make_unique<tools::Reporter_throughput<>>(*monitor));
    auto terminal = std::make_unique<tools::Terminal_std>(reporters);

    // ── Buffers ──────────────────────────────────────────────────────────
    std::vector<int>   ref_bits (p.K_bits);
    std::vector<int>   enc_bits (p.N_bits);
    std::vector<float> symbols  (p.N_bits);
    std::vector<float> noisy    (p.N_bits);
    std::vector<float> llrs     (p.N_bits);
    std::vector<int>   dec_bits (p.K_bits);
    std::vector<float> sigma_vec(1);

    terminal->legend();

    // ── SNR loop ─────────────────────────────────────────────────────────
    for (float ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
    {
        const float esn0  = tools::ebn0_to_esn0(ebn0, p.R);
        const float sigma = tools::esn0_to_sigma(esn0);
        sigma_vec[0] = sigma;
        noise->set_values(sigma, ebn0, esn0);

        terminal->start_temp_report();

        while (!monitor->fe_limit_achieved() && !terminal->is_interrupt())
        {
            source ->generate    (ref_bits                     );
            encoder->encode      (ref_bits,  enc_bits          );
            modem  ->modulate     (enc_bits,  symbols           );
            channel->add_noise    (sigma_vec, symbols,  noisy   );
            modem  ->demodulate   (sigma_vec, noisy,    llrs    );
            decoder->decode_siho  (llrs,      dec_bits          );
            monitor->check_errors (dec_bits,  ref_bits          );
        }

        terminal->final_report();
        monitor->reset();
        terminal->reset();
        if (terminal->is_over()) break;
    }

    std::cout << "# End of simulation\n";
    return 0;
}
