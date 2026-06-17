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
//	int   K = 128;     // number of information bits
//	int   n_ff = 0;
//	int   N = 0;
//	int   fe = 100;     // number of frame errors
//	int   seed = 0;     // PRNG seed for the AWGN channel
//	float ebn0_min = 0.00f; // minimum SNR value
//	float ebn0_max = 10.01f; // maximum SNR value
//	float ebn0_step = 0.50f; // SNR step
//	float R;                   // code rate (R=K/N)
//	std::vector<int> poly = { 011, 013 };
//
//
//	//А) По длине кодового ограничения(K)
//	//	K = 3 (5, 7) : Простой, быстрый, малая память
//
//	//	K = 4 (11, 13) : Лучше на 0.5 - 1 дБ, стандарт LTE
//
//	//	K = 7 (91, 121) : Еще лучше на 1 - 2 дБ, но сложнее
//
//	//	Б) По порождающим полиномам
//	//	Для одного K = 4 сравните:
//
//	//{13, 15} (стандартный)
//
//	//{
//	//	13, 17
//	//} (другой)
//
//	//	{
//	//		15, 17
//	//	} (симметричный)
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
//	//std::unique_ptr < module::Encoder_RSC_generic_json_sys<>> encoder;
//	std::unique_ptr < module::Encoder_RSC_generic_sys<>> encoder;
//	//std::unique_ptr < module::Encoder_RSC_sys<>> encoder;
//	//std::unique_ptr < module::Decoder_Viterbi_list_parallel<>> decoder;
//	std::unique_ptr < module::Decoder_Viterbi_SIHO<int, float>> decoder;
//	// 
//	//std::unique_ptr < module::Decoder_RSC_BCJR<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_intra<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_intra_fast_x2_AVX<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_intra_fast_x2_SSE<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_std<>> decoder;
//	// 
//	//std::unique_ptr < module::Decoder_RSC_BCJR_inter_very_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_intra<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_intra_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_intra_std<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_generic<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_generic_std<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_generic_std_json<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_scan<>> decoder;
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_std<>> decoder;
//	// 
//	//std::unique_ptr < module::Decoder_RSC_BCJR_seq_very_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_chase_pyndiah<>> decoder;
//	//std::unique_ptr < module::Decoder_chase_pyndiah_fast<>> decoder;
//	//std::unique_ptr < module::Decoder_chase_std<>> decoder;
//	//std::unique_ptr < module::Decoder_ML_naive<>> decoder;
//	//std::unique_ptr < module::Decoder_ML_std<>> decoder;
//
//};
//void init_modules(const params& p, modules& m);
//
//struct buffers
//{
//	std::vector<int     > ref_bits;
//	std::vector<uint32_t> ref_count; // unused in this code, just to avoid run failure in the source
//	std::vector<int     > enc_bits;
//	std::vector<float   > symbols;
//	std::vector<float   > sigma;
//	std::vector<float   > noisy_symbols;
//	std::vector<float   > LLRs;
//	std::vector<int     > dec_bits;
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
//			m.source->generate(b.ref_bits, b.ref_count);
//			m.encoder->encode(b.ref_bits, b.enc_bits);
//			m.modem->modulate(b.enc_bits, b.symbols);
//			m.channel->add_noise(b.sigma, b.symbols, b.noisy_symbols);
//			m.modem->demodulate(b.sigma, b.noisy_symbols, b.LLRs);
//			m.decoder->decode_siho(b.LLRs, b.dec_bits);
//			m.monitor->check_errors(b.dec_bits, b.ref_bits);
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
//	p.n_ff = (int)(std::floor(std::log2(
//		*std::max_element(p.poly.begin(), p.poly.end()))));
//
//	// N = 2*(K + n_ff): K сист. + K парит. + 2*n_ff хвост
//	p.N = 2 * (p.K + p.n_ff);
//	p.R = (float)p.K / (float)p.N;
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
//	std::cout << "#    ** Decoder_Viterbi_SIHO " << std::endl;
//	std::cout << "#" << std::endl;
//}
//
//void init_modules(const params& p, modules& m)
//{
//	m.source = std::unique_ptr<spu::module::Source_random         <>>(new spu::module::Source_random         <>(p.K));
//	m.encoder = std::unique_ptr<     module::Encoder_RSC_generic_sys<>>(new      module::Encoder_RSC_generic_sys<>(p.K, p.N, p.buffered_encoding, p.poly));
//	m.modem = std::unique_ptr<     module::Modem_BPSK            <>>(new      module::Modem_BPSK            <>(p.N));
//	const auto& trellis = m.encoder->get_trellis();
//	m.channel = std::unique_ptr<     module::Channel_AWGN_LLR      <>>(new      module::Channel_AWGN_LLR      <>(p.N));
//	m.decoder = std::unique_ptr<     module::Decoder_Viterbi_SIHO<int, float>>(new      module::Decoder_Viterbi_SIHO<int, float>(p.K, trellis, p.buffered_encoding));
//	m.monitor = std::unique_ptr<     module::Monitor_BFER          <>>(new      module::Monitor_BFER          <>(p.K, p.fe));
//	m.channel->set_seed(p.seed);
//};
//
//void init_buffers(const params& p, buffers& b)
//{
//	b.ref_bits = std::vector<int     >(p.K);
//	b.ref_count = std::vector<uint32_t>(1);
//	b.enc_bits = std::vector<int     >(p.N);
//	b.symbols = std::vector<float   >(p.N);
//	b.sigma = std::vector<float   >(1);
//	b.noisy_symbols = std::vector<float   >(p.N);
//	b.LLRs = std::vector<float   >(p.N);
//	b.dec_bits = std::vector<int     >(p.K);
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