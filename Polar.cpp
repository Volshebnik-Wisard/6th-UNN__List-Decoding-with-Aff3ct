#define _SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING 
#define _SILENCE_ALL_MS_EXT_DEPRECATION_WARNINGS
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include "codec_polar.hpp"
#include <aff3ct.hpp>
using namespace aff3ct;
float current_snr_for_frozen = 0.0f;
struct params
{
	int   K = 448;     // number of information bits
	int   N = 512;     // codeword size
	int   fe = 100;     // number of frame errors
	int   seed = 0;     // PRNG seed for the AWGN channel
	float ebn0_min = 0.00f; // minimum SNR value
	float ebn0_max = 10.01f; // maximum SNR value
	float ebn0_step = 0.50f; // SNR step
	float R;                   // code rate (R=K/N)
	int L = 16; //L — это размер списка (list size) в списочном декодировании.
	int crc_size = 16;
	int K_data;             // реальные данные (K - crc_size)
	std::vector<bool> frozen_bits;
};
void init_params(params& p);

struct modules
{
	std::unique_ptr<module::Source_random<>>          source;
	std::unique_ptr<module::Encoder_polar_sys<>>              encoder;

	std::unique_ptr<module::Modem_BPSK<>>             modem;
	std::unique_ptr<module::Channel_AWGN_LLR<>>       channel;
	//std::unique_ptr<module::Decoder_polar_SC_fast_sys<>>          decoder; //Обычное декодирование (базовый уровень)
	// Списочное декодирование (основной объект исследования)
	//std::unique_ptr<module::Decoder_polar_SCL_fast_sys<>>         decoder; //SCL без CA (CRC-Aided) - менее эффективные версии
	std::unique_ptr<module::Decoder_polar_SCL_fast_CA_sys<>>      decoder; //CA-SCL дает значительный выигрыш, поэтому сравнивать стоит именно с CA
	//std::unique_ptr<module::Decoder_polar_ASCL_fast_CA_sys<>>     decoder; //адаптивное списочное (для продвинутого исследования)
	std::unique_ptr<module::CRC_polynomial<>>                    crc;
	std::unique_ptr<module::Monitor_BFER<>>           monitor;
};
void init_modules(const params& p, modules& m);

struct buffers
{
	std::vector<int> src_bits;           // K_data бит (без CRC)
	std::vector<int> src_bits_with_crc;  // K бит (с CRC)
	std::vector<int> enc_bits;           // N бит
	std::vector<float> symbols;          // N символов
	std::vector<float> sigma;            // 1 значение
	std::vector<float> noisy_symbols;    // N символов
	std::vector<float> LLRs;             // N LLR
	std::vector<int> dec_bits;           // K бит (данные + CRC)
	std::vector<int> data_bits_only;     // K_data бит (только данные)
};
void init_buffers(const params& p, buffers& b);

struct utils
{
	std::unique_ptr<tools::Sigma<>>               noise;     // a sigma noise type
	std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std>          terminal;  // manage the output text in the terminal
};
void init_utils(const modules& m, utils& u);

int main(int argc, char** argv)
{
	// get the AFF3CT version
	const std::string v = "v" + std::to_string(tools::version_major()) + "." +
		std::to_string(tools::version_minor()) + "." +
		std::to_string(tools::version_release());

	std::cout << "#----------------------------------------------------------" << std::endl;
	std::cout << "# This is a basic program using the AFF3CT library (" << v << ")" << std::endl;
	std::cout << "# Feel free to improve it as you want to fit your needs." << std::endl;
	std::cout << "#----------------------------------------------------------" << std::endl;
	std::cout << "#" << std::endl;

	params p;  init_params(p); // create and initialize the parameters defined by the user
	modules m; init_modules(p, m); // create and initialize the modules
	buffers b; init_buffers(p, b); // create and initialize the buffers required by the modules
	utils u;   init_utils(m, u); // create and initialize the utils

	// display the legend in the terminal
	u.terminal->legend();

	// loop over the various SNRs
	for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
	{
		// compute the current sigma for the channel noise
		const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
		std::fill(b.sigma.begin(), b.sigma.end(), tools::esn0_to_sigma(esn0));

		u.noise->set_values(b.sigma[0], ebn0, esn0);

		// display the performance (BER and FER) in real time (in a separate thread)
		u.terminal->start_temp_report();

		// run the simulation chain
		while (!m.monitor->fe_limit_achieved() && !u.terminal->is_interrupt())
		{
// 1. Генерация K_data бит
            m.source->generate(b.src_bits);
            
            // 2. Добавление CRC к данным
            m.crc->build(b.src_bits, b.src_bits_with_crc);
            
            // 3. Полярное кодирование (K бит -> N бит)
            m.encoder->encode(b.src_bits_with_crc, b.enc_bits);
            
            // 4. BPSK модуляция
            m.modem->modulate(b.enc_bits, b.symbols);
            
            // 5. AWGN канал
            m.channel->add_noise(b.sigma, b.symbols, b.noisy_symbols);
            
            // 6. BPSK демодуляция -> LLR
            m.modem->demodulate(b.sigma, b.noisy_symbols, b.LLRs);
            
            // 7. CA-SCL декодирование (возвращает K бит: данные + CRC)
            m.decoder->decode_siho(b.LLRs, b.dec_bits);
            
            // 8. Извлекаем только данные (первые K_data бит)
            std::copy_n(b.dec_bits.begin(), p.K_data, b.data_bits_only.begin());
            
            // 9. Проверка ошибок (только данные)
            m.monitor->check_errors(b.data_bits_only, b.src_bits);
		}

		// display the performance (BER and FER) in the terminal
		u.terminal->final_report();

		// reset the monitor for the next SNR
		m.monitor->reset();
		u.terminal->reset();

		// if user pressed Ctrl+c twice, exit the SNRs loop
		if (u.terminal->is_over()) break;
	}

	std::cout << "# End of the simulation" << std::endl;

	return 0;
}

void init_params(params& p)
{
	p.R = (float)p.K / (float)p.N;
	p.K_data = p.K - p.crc_size;  // реальные информационные биты
	std::cout << "# * Simulation parameters: " << std::endl;
	std::cout << "#    ** Frame errors   = " << p.fe << std::endl;
	std::cout << "#    ** Noise seed     = " << p.seed << std::endl;
	std::cout << "#    ** Info. bits (K) = " << p.K << std::endl;
	std::cout << "#    ** Frame size (N) = " << p.N << std::endl;
	std::cout << "#    ** Code rate  (R) = " << p.R << std::endl;
	std::cout << "#    ** SNR min   (dB) = " << p.ebn0_min << std::endl;
	std::cout << "#    ** SNR max   (dB) = " << p.ebn0_max << std::endl;
	std::cout << "#    ** SNR step  (dB) = " << p.ebn0_step << std::endl << std::endl;
	std::cout << "#    ** Encoder_polar_sys " << std::endl;
	std::cout << "#    ** Decoder_polar_SCL_fast_CA_sys " << std::endl;
	std::cout << "#    ** List Size = " << p.L << std::endl;
	std::cout << "#    ** CRC Size = " << p.crc_size << std::endl;
	std::cout << "#    ** Polynomial = " << "16-5G; " "0x1023" << std::endl << std::endl;
	std::cout << "#" << std::endl;
	p.frozen_bits = generate_frozen_bits(p.K, p.N, p.ebn0_max * 0.4);
}

void init_modules(const params& p, modules& m)
{
	m.crc = std::unique_ptr<module::CRC_polynomial         <>>(new module::CRC_polynomial         <>(p.K_data,"0x1023", p.crc_size));
	//m.CRC->build(p.K, p.K_total);
	m.source = std::unique_ptr<module::Source_random         <>>(new module::Source_random         <>(p.K_data));
	m.encoder = std::unique_ptr<module::Encoder_polar_sys<>>(new module::Encoder_polar_sys<>(p.K, p.N, p.frozen_bits));
	m.modem = std::unique_ptr<module::Modem_BPSK            <>>(new module::Modem_BPSK            <>(p.N));
	m.channel = std::unique_ptr<module::Channel_AWGN_LLR      <>>(new module::Channel_AWGN_LLR      <>(p.N));
	m.decoder = std::unique_ptr<module::Decoder_polar_SCL_fast_CA_sys<>>(new module::Decoder_polar_SCL_fast_CA_sys<>(p.K, p.N, p.L, p.frozen_bits, *m.crc));
	m.monitor = std::unique_ptr<module::Monitor_BFER          <>>(new module::Monitor_BFER          <>(p.K_data, p.fe));
	m.channel->set_seed(p.seed);
};

void init_buffers(const params& p, buffers& b)
{
	b.src_bits = std::vector<int>(p.K_data);           // данные
	b.src_bits_with_crc = std::vector<int>(p.K);       // данные + CRC
	b.enc_bits = std::vector<int>(p.N);
	b.symbols = std::vector<float>(p.N);
	b.sigma = std::vector<float>(1);
	b.noisy_symbols = std::vector<float>(p.N);
	b.LLRs = std::vector<float>(p.N);
	b.dec_bits = std::vector<int>(p.K);                // данные + CRC
	b.data_bits_only = std::vector<int>(p.K_data);     // только данные
}

void init_utils(const modules& m, utils& u)
{
	// create a sigma noise type
	u.noise = std::unique_ptr<tools::Sigma<>>(new tools::Sigma<>());
	// report the noise values (Es/N0 and Eb/N0)
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_noise<>(*u.noise)));
	// report the bit/frame error rates
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*m.monitor)));
	// report the simulation throughputs
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*m.monitor)));
	// create a terminal that will display the collected data from the reporters
	u.terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(u.reporters));
}