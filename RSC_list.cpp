//#include <iostream>
//#include <memory>
//#include <vector>
//#include <string>
//#define NOMINMAX
////#include <streampu.hpp>
//#include <aff3ct.hpp>
//using namespace aff3ct;
//struct params
//{
//	int   K = 1024;     // number of information bits
//	int   n_ff = 0;
//	int   N = 0;
//	int   fe = 100;     // number of frame errors
//	int   seed = 0;     // PRNG seed for the AWGN channel
//	float ebn0_min = 0.00f; // minimum SNR value
//	float ebn0_max = 10.01f; // maximum SNR value
//	float ebn0_step = 0.50f; // SNR step
//	float R;                   // code rate (R=K/N)
//
//	int L = 16; //L — это размер списка (list size) в списочном декодировании.
//	int crc_size = 16;
//	int K_data = 0;             // реальные данные (K - crc_size)
//
//	std::vector<int> poly = { 0133, 0171 };
//	bool buffered_encoding = false;
//};
//void init_params(params& p);
//
//struct modules
//{
//	std::unique_ptr<spu::module::Source_random<>>          source;
//	//std::unique_ptr<     module::Encoder_repetition_sys<>> encoder;
//	std::unique_ptr<     module::Modem_BPSK<>>             modem;
//	std::unique_ptr<     module::Channel_AWGN_LLR<>>       channel;
//	//std::unique_ptr<     module::Decoder_repetition_std<>> decoder;
//	std::unique_ptr<     module::Monitor_BFER<>>           monitor;
//
//
//	std::unique_ptr<module::CRC_polynomial<>>                    crc;
//
//	//std::unique_ptr < module::Encoder_RSC_generic_json_sys<>> encoder;
//	std::unique_ptr < module::Encoder_RSC_generic_sys<>> encoder;
//	//std::unique_ptr < module::Encoder_RSC_sys<>> encoder;
//	std::unique_ptr < module::Decoder_Viterbi_list_parallel<int, float>> decoder;
//
//};
//void init_modules(const params& p, modules& m);
//
//struct buffers
//{
//	std::vector<int     > ref_bits;
//	std::vector<uint32_t> ref_count; // unused in this code, just to avoid run failure in the source
//	std::vector<int> src_bits;           // K_data бит (без CRC)
//	std::vector<int> src_bits_with_crc;  // K бит (с CRC)
//	std::vector<int> dec_bits;           // K бит (данные + CRC)
//	std::vector<int> data_bits_only;     // K_data бит (только данные)
//	std::vector<int     > enc_bits;
//	std::vector<float   > symbols;
//	std::vector<float   > sigma;
//	std::vector<float   > noisy_symbols;
//	std::vector<float   > LLRs;
//};
//void init_buffers(const params& p, buffers& b);
//
//struct utils
//{
//	std::unique_ptr<tools::Sigma<>>                    noise;     // a sigma noise type
//	std::vector<std::unique_ptr<spu::tools::Reporter>> reporters; // list of reporters dispayed in the terminal
//	std::unique_ptr<spu::tools::Terminal_std>          terminal;  // manage the output text in the terminal
//};
//void init_utils(const modules& m, utils& u);
//
//int main(int argc, char** argv)
//{
//	// StreamPU will catch and manage sigint
//	spu::tools::Signal_handler::init();
//
//	// get the AFF3CT version
//	const std::string v = "v" + std::to_string(tools::version_major()) + "." +
//		std::to_string(tools::version_minor()) + "." +
//		std::to_string(tools::version_release());
//
//	std::cout << "#----------------------------------------------------------" << std::endl;
//	std::cout << "# This is a basic program using the AFF3CT library (" << v << ")" << std::endl;
//	std::cout << "# Feel free to improve it as you want to fit your needs." << std::endl;
//	std::cout << "#----------------------------------------------------------" << std::endl;
//	std::cout << "#" << std::endl;
//
//	params p;  init_params(p); // create and initialize the parameters defined by the user
//	modules m; init_modules(p, m); // create and initialize the modules
//	buffers b; init_buffers(p, b); // create and initialize the buffers required by the modules
//	utils u;   init_utils(m, u); // create and initialize the utils
//
//	// display the legend in the terminal
//	u.terminal->legend();
//
//	// loop over the various SNRs
//	for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
//	{
//		// compute the current sigma for the channel noise
//		const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
//		std::fill(b.sigma.begin(), b.sigma.end(), tools::esn0_to_sigma(esn0));
//
//		u.noise->set_values(b.sigma[0], ebn0, esn0);
//
//		// display the performance (BER and FER) in real time (in a separate thread)
//		u.terminal->start_temp_report();
//
//		// run the simulation chain
//		while (!m.monitor->fe_limit_achieved() && !spu::tools::Signal_handler::is_sigint())
//		{
//
//			// 1. Генерация K_data бит
//			m.source->generate(b.src_bits, b.ref_count);
//
//			// 2. Добавление CRC: K_data → K_data + crc_size
//			m.crc->build(b.src_bits, b.src_bits_with_crc);
//
//			// 3. RSC кодирование: K → N
//			m.encoder->encode(b.src_bits_with_crc, b.enc_bits);
//
//			// 4. BPSK модуляция
//			m.modem->modulate(b.enc_bits, b.symbols);
//
//			// 5. AWGN канал
//			m.channel->add_noise(b.sigma, b.symbols, b.noisy_symbols);
//
//			// 6. BPSK демодуляция → LLR
//			m.modem->demodulate(b.sigma, b.noisy_symbols, b.LLRs);
//
//			// 7. Списочное декодирование Витерби (возвращает K бит: данные + CRC)
//			m.decoder->decode_siho(b.LLRs, b.dec_bits);
//
//			// 8. Извлекаем только K_data бит (без CRC-хвоста)
//			std::copy_n(b.dec_bits.begin(), p.K_data, b.data_bits_only.begin());
//
//			// 9. Проверка ошибок
//			m.monitor->check_errors(b.data_bits_only, b.src_bits);
//		}
//
//		// reset sigint if previously triggered
//		if (spu::tools::Signal_handler::is_sigint())
//			spu::tools::Signal_handler::reset_sigint();
//
//		// display the performance (BER and FER) in the terminal
//		u.terminal->final_report();
//
//		// reset the monitor for the next SNR
//		m.monitor->reset();
//	}
//
//	std::cout << "# End of the simulation" << std::endl;
//
//	return 0;
//}
//
//void init_params(params& p)
//{
//
//
//
//	p.K_data = p.K - p.crc_size;  // реальные информационные биты
//
//	p.n_ff = (int)(std::floor(std::log2(
//		*std::max_element(p.poly.begin(), p.poly.end()))));
//
//	// N = 2*(K + n_ff): K сист. + K парит. + 2*n_ff хвост
//	p.N = 2 * (p.K + p.n_ff);
//
//	p.R = (float)p.K / (float)p.N;
//
//	std::cout << "# * Simulation parameters: " << std::endl;
//	std::cout << "#    ** Frame errors   = " << p.fe << std::endl;
//	std::cout << "#    ** Noise seed     = " << p.seed << std::endl;
//	std::cout << "#    ** Info. bits    (K) = " << p.K << "\n";
//	std::cout << "#    ** Frame size    (N) = " << p.N << "\n";
//	std::cout << "#    ** Flip-flops (n_ff) = " << p.n_ff << "\n";
//	std::cout << "#    ** Trellis states    = " << (1 << p.n_ff) << "\n";
//	std::cout << "#    ** Code rate     (R) = " << p.R << "\n";
//	std::cout << "#    ** Polynomials       = {"
//		<< std::oct << p.poly[0] << ", " << p.poly[1]
//		<< std::dec << "} (octal)\n";
//	std::cout << "#    ** SNR min   (dB) = " << p.ebn0_min << std::endl;
//	std::cout << "#    ** SNR max   (dB) = " << p.ebn0_max << std::endl;
//	std::cout << "#    ** SNR step  (dB) = " << p.ebn0_step << std::endl;
//	std::cout << "#    ** Encoder_RSC_generic_sys " << std::endl;
//	std::cout << "#    ** Decoder_Viterbi_list_parallel " << std::endl;
//	std::cout << "#    ** List Size = " << p.L << std::endl;
//	std::cout << "#    ** CRC Size = " << p.crc_size << std::endl;
//	std::cout << "#    ** Polynomial = " << "16-5G; " "0x1023" << std::endl << std::endl;
//	std::cout << "#" << std::endl;
//}
//
//void init_modules(const params& p, modules& m)
//{
//	m.crc = std::unique_ptr<module::CRC_polynomial         <>>(new module::CRC_polynomial         <>(p.K_data, "0x1023", p.crc_size));
//	//m.CRC->build(p.K, p.K_total);
//	m.source = std::unique_ptr<spu::module::Source_random         <>>(new spu::module::Source_random         <>(p.K_data));
//	m.encoder = std::unique_ptr<     module::Encoder_RSC_generic_sys<>>(new      module::Encoder_RSC_generic_sys<>(p.K, p.N, false, p.poly));
//	/*std::cout << "# DEBUG encoder: K=" << m.encoder->get_K()
//		<< " N=" << m.encoder->get_N() << std::endl;*/
//	m.modem = std::unique_ptr<     module::Modem_BPSK            <>>(new      module::Modem_BPSK            <>(p.N));
//	const auto& trellis = m.encoder->get_trellis();
//	m.channel = std::unique_ptr<     module::Channel_AWGN_LLR      <>>(new      module::Channel_AWGN_LLR      <>(p.N));
//	//m.decoder = std::unique_ptr<     module::Decoder_Viterbi_SIHO<int, float>>(new      module::Decoder_Viterbi_SIHO<int, float>(p.K, trellis, true));
//	m.decoder = std::unique_ptr<     module::Decoder_Viterbi_list_parallel<int, float>>(new      module::Decoder_Viterbi_list_parallel<int, float>(p.K, p.N, p.L, *m.crc, trellis, true));
//	//std::cout << "# DEBUG decoder: K=" << m.decoder->get_K()
//	//	<< " N=" << m.decoder->get_N() << std::endl;
//	m.monitor = std::unique_ptr<     module::Monitor_BFER          <>>(new      module::Monitor_BFER          <>(p.K_data, p.fe));
//	m.channel->set_seed(p.seed);
//};
//
//void init_buffers(const params& p, buffers& b)
//{
//	b.src_bits = std::vector<int>(p.K_data);           // данные
//	b.src_bits_with_crc = std::vector<int>(p.K);       // данные + CRC
//	b.dec_bits = std::vector<int>(p.K);                // данные + CRC
//	b.data_bits_only = std::vector<int>(p.K_data);     // только данные
//	b.ref_bits = std::vector<int     >(p.K);
//	b.ref_count = std::vector<uint32_t>(1);
//	b.enc_bits = std::vector<int     >(p.N);
//	b.symbols = std::vector<float   >(p.N);
//	b.sigma = std::vector<float   >(1);
//	b.noisy_symbols = std::vector<float   >(p.N);
//	b.LLRs = std::vector<float   >(p.N);
//}
//
//void init_utils(const modules& m, utils& u)
//{
//	// create a sigma noise type
//	u.noise = std::unique_ptr<tools::Sigma<>>(new tools::Sigma<>());
//	// report the noise values (Es/N0 and Eb/N0)
//	u.reporters.push_back(std::unique_ptr<spu::tools::Reporter>(new tools::Reporter_noise<>(*u.noise)));
//	// report the bit/frame error rates
//	u.reporters.push_back(std::unique_ptr<spu::tools::Reporter>(new tools::Reporter_BFER<>(*m.monitor)));
//	// report the simulation throughput
//	u.reporters.push_back(std::unique_ptr<spu::tools::Reporter>(new tools::Reporter_throughput<>(*m.monitor)));
//	// create a terminal that will display the collected data from the reporters
//	u.terminal = std::unique_ptr<spu::tools::Terminal_std>(new spu::tools::Terminal_std(u.reporters));
//}