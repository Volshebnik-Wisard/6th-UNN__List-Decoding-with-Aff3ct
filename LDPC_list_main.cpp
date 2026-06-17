/*!
 * \file LDPC_list_main.cpp
 *
 * BER/FER симуляция для Decoder_LDPC_BP_list<int,float>.
 * Написан по образцу LDPC_old.cpp — использует стандартные AFF3CT-модули.
 *
 * ╔══════════════════════════════════════════════════════╗
 * ║  Все настройки — в структуре params (см. ниже).     ║
 * ║  Аргументы командной строки не нужны.               ║
 * ╚══════════════════════════════════════════════════════╝
 */

#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <aff3ct.hpp>
using namespace aff3ct;

#include "Decoder_LDPC_BP_list.hpp"

// ─────────────────────────────────────────────────────────────────────────────
// Параметры симуляции — меняйте здесь
// ─────────────────────────────────────────────────────────────────────────────
struct params
{
    // ── Путь к матрице H и размер листа ─────────────────────────────────
    // Путь к .alist файлу. Используем std::wstring — это единственный
    // надёжный способ открыть файл с кириллицей в пути на Windows.
    // Префикс L"..." делает строку широкой (UTF-16 на MSVC).
    std::wstring alist_path =
        L"E:\\1111\\filename.alist";
    int list_size = 16;   // 1 (plain BP), 2, 4, 8, 16, 32

    // ── Параметры кода ───────────────────────────────────────────────────
    int   K         = 144;    // информационных бит  (= N - M = 271 - 127)
    int   N         = 271;    // длина кодового слова
    int   n_ite     = 10;     // итераций BP
    int   fe        = 100;    // кадров с ошибками до остановки
    int   seed      = 0;      // зерно ГПСЧ канала
    float ebn0_min  = 0.00f;
    float ebn0_max  = 15.00f;
    float ebn0_step = 0.50f;
    float R;                  // кодовая скорость (вычисляется в init_params)

    tools::Sparse_matrix   H;             // матрица проверок чётности

    // Максимальная степень check-узла для SPA (6 — как в LDPC_old.cpp)
    tools::Update_rule_SPA<> update_rule = tools::Update_rule_SPA<>(6);
};

void init_params(params& p);

// ─────────────────────────────────────────────────────────────────────────────
struct modules
{
    std::unique_ptr<module::Source_random<>>         source;
    std::unique_ptr<module::Encoder_LDPC<>>         encoder;
    std::unique_ptr<module::Modem_BPSK<>>            modem;
    std::unique_ptr<module::Channel_AWGN_LLR<>>      channel;
    std::unique_ptr<module::Decoder_LDPC_BP_list<>>  decoder;
    std::unique_ptr<module::Monitor_BFER<>>           monitor;
};

void init_modules(const params& p, modules& m);

// ─────────────────────────────────────────────────────────────────────────────
struct buffers
{
    std::vector<int  > ref_bits;
    std::vector<int  > enc_bits;
    std::vector<float> symbols;
    std::vector<float> sigma;
    std::vector<float> noisy_symbols;
    std::vector<float> LLRs;
    std::vector<int  > dec_bits;
};

void init_buffers(const params& p, buffers& b);

// ─────────────────────────────────────────────────────────────────────────────
struct utils
{
    std::unique_ptr<tools::Sigma<>>               noise;
    std::vector<std::unique_ptr<tools::Reporter>> reporters;
    std::unique_ptr<tools::Terminal_std>          terminal;
};

void init_utils(const modules& m, utils& u);

// ─────────────────────────────────────────────────────────────────────────────
int main()
{
    const std::string v = "v" + std::to_string(tools::version_major()) + "." +
                               std::to_string(tools::version_minor())   + "." +
                               std::to_string(tools::version_release());

    std::cout << "#----------------------------------------------------------\n"
              << "# AFF3CT (" << v << ")  —  LDPC List-BP decoder demo\n"
              << "#----------------------------------------------------------\n"
              << "#\n";

    params  p; init_params (p);
    modules m; init_modules(p, m);
    buffers b; init_buffers(p, b);
    utils   u; init_utils  (m, u);

    u.terminal->legend();

    for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
    {
        const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
        std::fill(b.sigma.begin(), b.sigma.end(), tools::esn0_to_sigma(esn0));

        u.noise->set_values(b.sigma[0], ebn0, esn0);
        u.terminal->start_temp_report();

        while (!m.monitor->fe_limit_achieved() && !u.terminal->is_interrupt())
        {
            m.source ->generate   (b.ref_bits);
            m.encoder->encode     (b.ref_bits,        b.enc_bits);
            m.modem  ->modulate   (b.enc_bits,         b.symbols);
            m.channel->add_noise  (b.sigma, b.symbols, b.noisy_symbols);
            m.modem  ->demodulate (b.sigma, b.noisy_symbols, b.LLRs);
            m.decoder->decode_siho(b.LLRs,             b.dec_bits);
            m.monitor->check_errors(b.dec_bits,        b.ref_bits);
        }

        u.terminal->final_report();
        m.monitor->reset();
        u.terminal->reset();

        if (u.terminal->is_over()) break;
    }

    std::cout << "# End of the simulation\n";
    return 0;
}

// ─────────────────────────────────────────────────────────────────────────────
void init_params(params& p)
{
    p.R = static_cast<float>(p.K) / static_cast<float>(p.N);

    // Загрузка H. p.alist_path — wstring, std::ifstream на MSVC
    // принимает wstring напрямую, что корректно обрабатывает кириллицу.
#ifdef _WIN32
    std::ifstream h_file(p.alist_path); // MSVC: ifstream(wstring) — OK
#else
    std::ifstream h_file(std::string(p.alist_path.begin(), p.alist_path.end()));
#endif
    if (!h_file.is_open())
        throw std::runtime_error("Cannot open alist file: " +
            std::string(p.alist_path.begin(), p.alist_path.end()));
    p.H = tools::AList::read(h_file);
    h_file.close();

    const int n_amb = (p.list_size >= 2)
        ? static_cast<int>(std::floor(std::log2(static_cast<double>(p.list_size))))
        : 0;

    std::cout << "# * Simulation parameters:\n"
              << "#    ** H matrix path      = " << std::string(p.alist_path.begin(), p.alist_path.end()) << "\n"
              << "#    ** Info. bits  (K)    = " << p.K           << "\n"
              << "#    ** Frame size  (N)    = " << p.N           << "\n"
              << "#    ** Code rate   (R)    = " << p.R           << "\n"
              << "#    ** BP iterations      = " << p.n_ite       << "\n"
              << "#    ** List size          = " << p.list_size   << "\n"
              << "#    ** Ambiguous bits     = " << n_amb
                                                 << "  (2^" << n_amb
                                                 << " = " << (1 << n_amb)
                                                 << " candidates)\n"
              << "#    ** Frame errors       = " << p.fe          << "\n"
              << "#    ** Noise seed         = " << p.seed        << "\n"
              << "#    ** SNR min  (dB)      = " << p.ebn0_min    << "\n"
              << "#    ** SNR max  (dB)      = " << p.ebn0_max    << "\n"
              << "#    ** SNR step (dB)      = " << p.ebn0_step   << "\n"
              << "#\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Аналитика: статистика по матрице H и параметрам симуляции → analytics.txt
// ─────────────────────────────────────────────────────────────────────────────
static void write_analytics(const params& p,
                             const std::vector<uint32_t>& info_bits_pos)
{
    const std::string fname = "analytics.txt";
    std::ofstream f(fname);
    if (!f.is_open()) { std::cerr << "Cannot write " << fname << "\n"; return; }

    const int N = p.N, M = p.N - p.K, K = p.K;

    // Степени переменных узлов
    std::vector<int> var_deg(static_cast<size_t>(N));
    for (int v = 0; v < N; ++v)
        var_deg[static_cast<size_t>(v)] =
            static_cast<int>(p.H.get_cols_from_row(static_cast<size_t>(v)).size());

    // Степени check-узлов
    std::vector<int> chk_deg(static_cast<size_t>(M));
    for (int c = 0; c < M; ++c)
        chk_deg[static_cast<size_t>(c)] =
            static_cast<int>(p.H.get_rows_from_col(static_cast<size_t>(c)).size());

    const int var_deg_min = *std::min_element(var_deg.begin(), var_deg.end());
    const int var_deg_max = *std::max_element(var_deg.begin(), var_deg.end());
    const int chk_deg_min = *std::min_element(chk_deg.begin(), chk_deg.end());
    const int chk_deg_max = *std::max_element(chk_deg.begin(), chk_deg.end());
    const double var_deg_avg = static_cast<double>(p.H.get_n_connections()) / N;
    const double chk_deg_avg = static_cast<double>(p.H.get_n_connections()) / M;

    const int n_amb = (p.list_size >= 2)
        ? static_cast<int>(std::floor(std::log2(static_cast<double>(p.list_size))))
        : 0;

    f << "====================================================\n"
      << " LDPC List-BP — аналитика\n"
      << "====================================================\n\n"
      << "[Параметры кода]\n"
      << "  N (codeword length)   = " << N << "\n"
      << "  M (check nodes)       = " << M << "\n"
      << "  K (info bits)         = " << K << "\n"
      << "  Rate R = K/N          = " << p.R << "\n"
      << "  H connections total   = " << p.H.get_n_connections() << "\n\n"
      << "[Граф Таннера]\n"
      << "  Variable node degree  min=" << var_deg_min
                              << " max=" << var_deg_max
                              << " avg=" << var_deg_avg << "\n"
      << "  Check node degree     min=" << chk_deg_min
                              << " max=" << chk_deg_max
                              << " avg=" << chk_deg_avg << "\n\n"
      << "[List-BP]\n"
      << "  list_size             = " << p.list_size << "\n"
      << "  n_ambiguous bits      = " << n_amb << "\n"
      << "  candidates per frame  = " << (1 << n_amb) << "\n"
      << "  BP iterations (max)   = " << p.n_ite << "\n\n"
      << "[info_bits_pos  (первые 10 / последние 10)]\n  ";
    for (int i = 0; i < std::min(10, K); ++i)
        f << info_bits_pos[static_cast<size_t>(i)] << " ";
    f << "... ";
    for (int i = std::max(0, K - 10); i < K; ++i)
        f << info_bits_pos[static_cast<size_t>(i)] << " ";
    f << "\n\n"
      << "[Симуляция]\n"
      << "  Eb/N0 range           = " << p.ebn0_min << " .. "
                                       << p.ebn0_max << " dB"
                                       << " (step " << p.ebn0_step << ")\n"
      << "  frame errors limit    = " << p.fe  << "\n"
      << "  PRNG seed             = " << p.seed << "\n"
      << "====================================================\n";

    f.close();
    std::cout << "# Analytics written to: " << fname << "\n#\n";
}


// ─────────────────────────────────────────────────────────────────────────────
// GF(2) Gaussian elimination over H (M×N) to find:
//   • pivot_cols  — M parity-bit positions
//   • info_cols   — K=N-M information-bit positions
// Returns G as a Sparse_matrix (vertical, K rows × N cols) and fills info_bits_pos.
// This replaces Encoder_LDPC_from_H which crashes on matrices with no unit columns.
// ─────────────────────────────────────────────────────────────────────────────
static tools::Sparse_matrix
build_G_from_H(const tools::Sparse_matrix &H_sparse,
               int K, int N,
               std::vector<uint32_t> &info_bits_pos)
{
    const int M = N - K;

    // ── Dense H (M × N) ──────────────────────────────────────────────────
    std::vector<std::vector<uint8_t>> mat(
        static_cast<size_t>(M),
        std::vector<uint8_t>(static_cast<size_t>(N), 0));

    // H_sparse из AList::read: rows=variable nodes (N=271), cols=check nodes (M=127)
    //   get_cols_from_row(v) → список check-узлов, смежных с переменным узлом v
    // Нам нужна H[check_node][variable_node], поэтому:
    for (int v = 0; v < N; ++v)
        for (uint32_t c : H_sparse.get_cols_from_row(static_cast<size_t>(v)))
            mat[static_cast<size_t>(c)][static_cast<size_t>(v)] = 1;

    // ── GF(2) row reduction to find M pivot columns ───────────────────────
    std::vector<int> pivot_cols;
    pivot_cols.reserve(static_cast<size_t>(M));
    int cur_row = 0;
    for (int col = 0; col < N && cur_row < M; ++col)
    {
        // Find pivot
        int piv = -1;
        for (int r = cur_row; r < M; ++r)
            if (mat[static_cast<size_t>(r)][static_cast<size_t>(col)]) { piv = r; break; }
        if (piv < 0) continue;

        // Swap
        std::swap(mat[static_cast<size_t>(cur_row)], mat[static_cast<size_t>(piv)]);

        // Eliminate column in all other rows
        for (int r = 0; r < M; ++r)
        {
            if (r == cur_row) continue;
            if (!mat[static_cast<size_t>(r)][static_cast<size_t>(col)]) continue;
            for (int c2 = 0; c2 < N; ++c2)
                mat[static_cast<size_t>(r)][static_cast<size_t>(c2)] ^=
                    mat[static_cast<size_t>(cur_row)][static_cast<size_t>(c2)];
        }
        pivot_cols.push_back(col);
        ++cur_row;
    }

    // ── Info-bit positions = non-pivot columns ────────────────────────────
    std::vector<bool> is_pivot(static_cast<size_t>(N), false);
    for (int c : pivot_cols) is_pivot[static_cast<size_t>(c)] = true;

    info_bits_pos.clear();
    info_bits_pos.reserve(static_cast<size_t>(K));
    for (int c = 0; c < N; ++c)
        if (!is_pivot[static_cast<size_t>(c)])
            info_bits_pos.push_back(static_cast<uint32_t>(c));

    // ── Build G (K × N sparse, vertical) from reduced H ──────────────────
    // For each info bit k (column info_bits_pos[k]):
    //   G_row[k][info_bits_pos[k]] = 1                  (systematic part)
    //   G_row[k][pivot_cols[r]]    = mat[r][info_bits_pos[k]]  (parity part)
    tools::Sparse_matrix G(static_cast<size_t>(N),   // rows = N (vertical)
                           static_cast<size_t>(K));  // cols = K

    for (int k = 0; k < K; ++k)
    {
        const int info_col = static_cast<int>(info_bits_pos[static_cast<size_t>(k)]);
        // Systematic 1
        G.add_connection(static_cast<size_t>(info_col), static_cast<size_t>(k));
        // Parity entries
        for (int r = 0; r < M; ++r)
            if (mat[static_cast<size_t>(r)][static_cast<size_t>(info_col)])
                G.add_connection(static_cast<size_t>(pivot_cols[static_cast<size_t>(r)]),
                                 static_cast<size_t>(k));
    }

    return G;
}

void init_modules(const params& p, modules& m)
{
    m.source  = std::make_unique<module::Source_random<>>(p.K);
    m.modem   = std::make_unique<module::Modem_BPSK<>>(p.N);
    m.channel = std::make_unique<module::Channel_AWGN_LLR<>>(p.N);
    m.channel->set_seed(p.seed);

    // Строим G и info_bits_pos через GF(2) исключение —
    // оба метода AFF3CT (IDENTITY и LU_DEC) падают на данной матрице,
    // т.к. в ней нет ни одного единичного столбца (все переменные узлы степени 3).
    std::vector<uint32_t> info_bits_pos;
    tools::Sparse_matrix G = build_G_from_H(p.H, p.K, p.N, info_bits_pos);

    std::cout << "# info_bits_pos[0..3] = "
              << info_bits_pos[0] << " " << info_bits_pos[1] << " "
              << info_bits_pos[2] << " " << info_bits_pos[3] << " ...#";
    write_analytics(p, info_bits_pos);

    // Encoder_LDPC<> принимает готовую G (вертикальная ориентация)
    m.encoder = std::make_unique<module::Encoder_LDPC<>>(p.K, p.N, G);

    m.decoder = std::make_unique<module::Decoder_LDPC_BP_list<>>(
                    p.K, p.N, p.n_ite, p.H, info_bits_pos,
                    p.list_size,
                    true,  // enable_syndrome
                    1      // syndrome_depth
                );

    m.monitor = std::make_unique<module::Monitor_BFER<>>(p.K, p.fe);
}

void init_buffers(const params& p, buffers& b)
{
    b.ref_bits      = std::vector<int  >(static_cast<size_t>(p.K));
    b.enc_bits      = std::vector<int  >(static_cast<size_t>(p.N));
    b.symbols       = std::vector<float>(static_cast<size_t>(p.N));
    b.sigma         = std::vector<float>(1);
    b.noisy_symbols = std::vector<float>(static_cast<size_t>(p.N));
    b.LLRs          = std::vector<float>(static_cast<size_t>(p.N));
    b.dec_bits      = std::vector<int  >(static_cast<size_t>(p.K));
}

void init_utils(const modules& m, utils& u)
{
    u.noise = std::make_unique<tools::Sigma<>>();
    u.reporters.push_back(std::make_unique<tools::Reporter_noise<>>(*u.noise));
    u.reporters.push_back(std::make_unique<tools::Reporter_BFER <>>(*m.monitor));
    u.reporters.push_back(std::make_unique<tools::Reporter_throughput<>>(*m.monitor));
    u.terminal = std::make_unique<tools::Terminal_std>(u.reporters);
}
