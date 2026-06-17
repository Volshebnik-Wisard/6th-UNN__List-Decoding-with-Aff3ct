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
 * \file  main_frs1.cpp
 * \brief BER/FER simulation: Folded RS + FRS Decoder 1 (Algorithm 17.2.1).
 *
 * PARAMETERS:
 *   K_sym   — information symbols (RS)
 *   N_sym   — codeword length in symbols (RS); must be 2^m_gf - 1 or divisor
 *   t       — RS correction power  (N_sym - K_sym = 2*t)
 *   m_fold  — folding parameter (must divide N_sym)
 *   list_sz — list size limit (0 = unbounded, bounded by q^(m_fold-1))
 *
 * DECODING RADIUS (Theorem 17.2.3):
 *   delta = m/(m+1) * (1 - m*R)   where R = K/N (symbol rate)
 *
 * NOTE ON bit vs symbol domain:
 *   aff3ct works in bits.  RS/FRS encode K_sym*m_gf bits -> N_sym*m_gf bits.
 *   Monitor counts BIT errors over K_sym*m_gf info bits.
 */
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>

#include <aff3ct.hpp>
using namespace aff3ct;

#include "Encoder_FRS.hpp"
#include "Decoder_FRS1.hpp"
#include "FRS_List_Selector.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// Simulation parameters
// ─────────────────────────────────────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════
// ТАБЛИЦА ПАРАМЕТРОВ  (N_sym=255, GF(2^8), m_gf=8)
//
// Условие корректности Decoder 1:  R < 1/m_fold,  K_sym < N_frs = N_sym/m_fold
// Формула границы (Теорема 17.2.3): delta = m/(m+1) * (1 - m*R)
//
// ── Пара 1: K_sym=63, N_sym=255, R=0.247 ────────────────────────────────
//   RS baseline:  K=63, N=255, t=96
//
//   FRS Decoder 1:
//     m_fold=3 → N_frs=85, D=(85-63+1)/4=5,  delta=3/4*(1-3*0.247)=0.194  ✓
//     m_fold=5 → N_frs=51, D=(51-63+1)/6<0                                 ✗
//   Итого для Пары 1: K_sym=63, N_sym=255, m_fold=3
//
// ── Пара 2: K_sym=31, N_sym=255, R=0.122 ────────────────────────────────
//   RS baseline:  K=31, N=255, t=112
//
//   FRS Decoder 1:
//     m_fold=3  → N_frs=85, D=(85-31+1)/4=13,  delta=3/4*(1-3*0.122)=0.475 ✓
//     m_fold=5  → N_frs=51, D=(51-31+1)/6=3,   delta=5/6*(1-5*0.122)=0.322 ✓
//     m_fold=17 → N_frs=15, D=(15-31+1)/18<0                                ✗
//   Итого для Пары 2: K_sym=31, N_sym=255, m_fold=3 (или 5)
//
// ── Запуск из командной строки ───────────────────────────────────────────
//   cd E:\1-ННГУ\6-ВКР\2-Практика\1_RS_old
//   cmake --build build --config Release --target frs1
//   .\build\Release\frs1.exe 63 255 3 0 llr    <- Пара 1
//   .\build\Release\frs1.exe 31 255 3 0 llr    <- Пара 2, m=3
//   .\build\Release\frs1.exe 31 255 5 0 llr    <- Пара 2, m=5
// ═══════════════════════════════════════════════════════════════════════════

struct params
{
    // ── RS параметры ──────────────────────────────────────────────────────
    // Пара 1 (умеренная скорость): K_sym=63, N_sym=255, m_fold=3
    // Пара 2 (низкая скорость):   K_sym=31, N_sym=255, m_fold=3  (или 5)
    int K_sym   = 63;   // информационных символов RS
    int N_sym   = 255;  // длина кодового слова (= 2^8 - 1 для GF(2^8))
    // t выводится автоматически: t = (N_sym - K_sym) / 2

    // ── FRS параметры ─────────────────────────────────────────────────────
    int m_fold  = 3;    // параметр складывания; должен делить N_sym
                        // 255 = 3*5*17 → допустимые: 3, 5, 15, 17, 51, 85
    int list_sz = 2;    // максимальный размер списка (0 = без ограничений,
                        // теоретический максимум q^(m_fold-1) = 256^2 = 65536)

    // ── Параметры симуляции ───────────────────────────────────────────────
    int   fe        = 100;   // целевое число ошибочных кадров
    int   seed      = 0;
    // SNR диапазон: декодер работает только при доле ошибок < delta.
    // Для m_fold=3, K=63, N=255: delta=0.194 => порог примерно при Eb/N0 ~ 3-4 dB.
    // Для m_fold=1, K=63, N=255: delta=0.376 => порог при Eb/N0 ~ 1-2 dB.
    // Рекомендуется начинать с 0 dB и смотреть где BER начинает падать.
    float ebn0_min  = 6.0f;
    float ebn0_max  = 12.0f;
    float ebn0_step = 0.5f;

    // ── Стратегия выбора из списка ────────────────────────────────────────
    // "llr"     — максимум мягкой метрики LLR (лучше при AWGN)
    // "hamming" — минимум расстояния Хэмминга (запасной вариант)
    std::string selector = "llr";

    // ── Производные (заполняются в main) ─────────────────────────────────
    float R;
    int   t;
    int   m_gf;
    int   K_bits;
    int   N_bits;
};

// ─────────────────────────────────────────────────────────────────────────────
void print_params(const params& p)
{
    std::cout << "#----------------------------------------------------------\n";
    std::cout << "# FRS simulation — Decoder 1 (Algorithm 17.2.1)\n";
    std::cout << "#----------------------------------------------------------\n";
    std::cout << "#  RS symbols    K = " << p.K_sym << "\n";
    std::cout << "#  RS symbols    N = " << p.N_sym << "\n";
    std::cout << "#  Correction    t = " << p.t     << "\n";
    std::cout << "#  GF extension  m_gf = " << p.m_gf << "\n";
    std::cout << "#  Folding   m_fold = " << p.m_fold << "\n";
    std::cout << "#  FRS block N_frs = " << p.N_sym/p.m_fold << "\n";
    std::cout << "#  List limit    L = " << (p.list_sz ? std::to_string(p.list_sz) : "unlimited") << "\n";
    std::cout << "#  Selector      = " << p.selector << "\n";
    std::cout << "#  Code rate     R = " << p.R << "\n";

    // Theoretical decoding radius
    const double m = p.m_fold;
    const double R = p.R;
    const double delta = (m/(m+1.0)) * (1.0 - m*R);
    std::cout << "#  Theory delta (Thm 17.2.3) = " << delta << "\n";
    std::cout << "#\n";
}

// ─────────────────────────────────────────────────────────────────────────────
int main(int argc, char** argv)
{
    params p;

    // Simple command-line override: ./main_frs1 K N m_fold [list_sz] [selector]
    if (argc >= 4)
    {
        p.K_sym  = std::stoi(argv[1]);
        p.N_sym  = std::stoi(argv[2]);
        p.m_fold = std::stoi(argv[3]);
    }
    if (argc >= 5) p.list_sz  = std::stoi(argv[4]);
    if (argc >= 6) p.selector = argv[5];

    // Validate
    if (p.N_sym % p.m_fold != 0)
    {
        std::cerr << "ERROR: N_sym (" << p.N_sym << ") must be divisible by m_fold ("
                  << p.m_fold << ").\n";
        return 1;
    }
    if ((p.N_sym - p.K_sym) % 2 != 0)
    {
        std::cerr << "ERROR: N_sym - K_sym must be even (need integer t).\n";
        return 1;
    }
    p.t = (p.N_sym - p.K_sym) / 2;

    // Build GF — RS_polynomial_generator takes N_p2 = next_power_of_2(N_sym)-1
    const int N_p2 = tools::next_power_of_2(p.N_sym) - 1;
    tools::RS_polynomial_generator GF(N_p2, p.t);

    p.m_gf  = GF.get_m();
    p.K_bits = p.K_sym * p.m_gf;
    p.N_bits = p.N_sym * p.m_gf;
    p.R      = static_cast<float>(p.K_sym) / static_cast<float>(p.N_sym);

    print_params(p);

    // ── Build selector ────────────────────────────────────────────────────
    using Sel = FRS_List_Selector<int, float>;
    Sel::SelectorFn selector_fn;
    // Кандидаты хранятся как N_rs значений f(γ^0)..f(γ^{N-1}).
    // Информационные символы — позиции [n_rdncy .. N_rs-1].
    // Селектор LLR сравнивает информационную часть кандидата с мягкими LLR.
    // Селектор Hamming считает символьные ошибки в свёрнутом слове.
    if (p.selector == "llr")
    {
        const int n_rdncy_sym  = GF.get_n_rdncy();  // число символов чётности
        const int n_rdncy_bits = n_rdncy_sym * p.m_gf;
        selector_fn = Sel::make_max_llr(p.m_gf, n_rdncy_bits, p.K_bits);
    }
    else
    {
        selector_fn = Sel::make_min_hamming(p.m_fold);
    }

    // ── Modules ──────────────────────────────────────────────────────────
    auto source  = std::make_unique<module::Source_random<>>(p.K_bits);
    auto encoder = std::make_unique<module::Encoder_FRS<>>(p.K_sym, p.N_sym, GF, p.m_fold);
    auto decoder = std::make_unique<module::Decoder_FRS1<>>(
        p.K_sym, p.N_sym, GF, p.m_fold, p.list_sz, selector_fn);
    auto modem   = std::make_unique<module::Modem_BPSK<>>(p.N_bits);
    auto channel = std::make_unique<module::Channel_AWGN_LLR<>>(p.N_bits);
    auto monitor = std::make_unique<module::Monitor_BFER<>>(p.K_bits, p.fe);

    channel->set_seed(p.seed);

    // ── Noise / reporters ────────────────────────────────────────────────
    auto noise = std::make_unique<tools::Sigma<>>();
    std::vector<std::unique_ptr<tools::Reporter>> reporters;
    reporters.push_back(std::make_unique<tools::Reporter_noise<>>(*noise));
    reporters.push_back(std::make_unique<tools::Reporter_BFER<>>(*monitor));
    reporters.push_back(std::make_unique<tools::Reporter_throughput<>>(*monitor));
    auto terminal = std::make_unique<tools::Terminal_std>(reporters);

    // ── Buffers ──────────────────────────────────────────────────────────
    std::vector<int>   ref_bits    (p.K_bits);
    std::vector<int>   enc_bits    (p.N_bits);
    std::vector<float> symbols     (p.N_bits);
    std::vector<float> noisy       (p.N_bits);
    std::vector<float> llrs        (p.N_bits);
    std::vector<int>   dec_bits    (p.K_bits);
    std::vector<float> sigma_vec   (1);

    terminal->legend();

    // ── SNR loop ─────────────────────────────────────────────────────────
    for (float ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
    {
        const float esn0  = tools::ebn0_to_esn0(ebn0, p.R);
        const float sigma = tools::esn0_to_sigma(esn0);
        sigma_vec[0] = sigma;
        noise->set_values(sigma, ebn0, esn0);

        terminal->start_temp_report();

        static int main_dbg = 0;
        while (!monitor->fe_limit_achieved() && !terminal->is_interrupt())
        {
            source ->generate    (ref_bits              );
            encoder->encode      (ref_bits,  enc_bits   );

            if (main_dbg < 2)
            {
                std::vector<int> ref_sym(3), enc_sym(3);
                aff3ct::tools::Bit_packer::pack(
                    ref_bits.data(), ref_sym.data(), p.m_gf*3, 1, false, p.m_gf);
                aff3ct::tools::Bit_packer::pack(
                    enc_bits.data(), enc_sym.data(), p.m_gf*3, 1, false, p.m_gf);
                std::cout << "[MAIN] ref=" << ref_sym[0] << "," << ref_sym[1]
                          << "," << ref_sym[2] << "  enc=" << enc_sym[0]
                          << "," << enc_sym[1] << "," << enc_sym[2] << "\n";
            }

            modem  ->modulate     (enc_bits,  symbols          );
            channel->add_noise    (sigma_vec, symbols,  noisy  );
            modem  ->demodulate   (sigma_vec, noisy,    llrs   );
            decoder->decode_siho  (llrs,      dec_bits         );

            if (main_dbg < 2)
            {
                std::vector<int> dec_sym(3);
                aff3ct::tools::Bit_packer::pack(
                    dec_bits.data(), dec_sym.data(), p.m_gf*3, 1, false, p.m_gf);
                std::cout << "[MAIN] dec=" << dec_sym[0] << "," << dec_sym[1]
                          << "," << dec_sym[2] << "\n";
                ++main_dbg;
            }

            monitor->check_errors (dec_bits,  ref_bits         );
        }

        terminal->final_report();
        monitor->reset();
        terminal->reset();
        if (terminal->is_over()) break;
    }

    std::cout << "# End of simulation\n";
    return 0;
}
