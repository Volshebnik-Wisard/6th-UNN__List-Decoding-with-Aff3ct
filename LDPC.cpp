#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
using namespace aff3ct;

std::string LDPC_H_MATRIX_NAME("E:/1111/filename.alist");

auto get_h_alist(std::string file_name) {
	std::ifstream h_file(file_name);
	auto h_mat = aff3ct::tools::AList::read(h_file);
	h_file.close();
	return h_mat;
}

struct params
{
	int   K = 144;     // number of information bits
	int   N = 271;     // codeword size
	int n_ite = 10;        // number of decoder iterations
	int   t = (int)(N - K) / 2.;   // correction power
	int   fe = 100;     // number of frame errors
	int   seed = 0;     // PRNG seed for the AWGN channel
	float ebn0_min = 0.00f; // minimum SNR value
	float ebn0_max = 15.00f; // maximum SNR value
	float ebn0_step = 0.5f; // SNR step
	float R;                   // code rate (R=K/N)
	tools::Sparse_matrix H; // parity check matrix
	// The maximum check-node degree for the AR4JA code used is 6
	Update_rule_SPA<> update_rule = tools::Update_rule_SPA<>(6);
	std::vector<uint32_t> info_bits_pos;
	RS_polynomial_generator Pol_Gen{ (1 << (int)std::ceil(std::log2(N))) - 1, t };
	//std::unique_ptr<dvbs2_values> DVBS2 = tools::build_dvbs2(K,N);
};
void init_params(params& p);

struct modules
{
	std::unique_ptr<module::Source_random<>>          source;
	//std::unique_ptr<module::Encoder_repetition_sys<>> encoder;
	//std::unique_ptr<module::Encoder_RS<>>             encoder;
	//std::unique_ptr<module::Encoder_LDPC<>>             encoder;
	std::unique_ptr<module::Encoder_LDPC_from_H<>>             encoder;
	//std::unique_ptr<module::Encoder_LDPC_DVBS2<>>             encoder;
	//std::unique_ptr<module::Encoder_LDPC_from_IRA<>>             encoder;
	//std::unique_ptr<module::Encoder_LDPC_from_QC<>>             encoder;
	std::unique_ptr<module::Modem_BPSK<>>             modem;
	std::unique_ptr<module::Channel_AWGN_LLR<>>       channel;
	//std::unique_ptr<module::Decoder_repetition_std<>> decoder;
	//std::unique_ptr<module::Decoder_RS_std<>>     decoder_std;
	//std::unique_ptr<module::Decoder_LDPC_bit_flipping<>>     decoder;           //abstact
	//std::unique_ptr<module::Decoder_LDPC_bit_flipping_hard<>>     decoder;      //abstact
	//std::unique_ptr<module::Decoder_LDPC_bit_flipping_OMWBF<>>     decoder;
	std::unique_ptr<module::Decoder_LDPC_BP_flooding<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_GALA<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_GALB<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_GALE<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_Gallager_A<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_Gallager_B<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_Gallager_E<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_flooding_SPA<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_horizontal_layered<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_peeling<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_BP_vertical_layered<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_PPBF<>>     decoder;
	//std::unique_ptr<module::Decoder_LDPC_probabilistic_parallel_bit_flipping<>>     decoder;
	std::unique_ptr<module::Monitor_BFER<>>           monitor;
};
void init_modules(const params& p, modules& m);

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
			m.source->generate(b.ref_bits);
			m.encoder->encode(b.ref_bits, b.enc_bits);
			m.modem->modulate(b.enc_bits, b.symbols);
			m.channel->add_noise(b.sigma, b.symbols, b.noisy_symbols);
			m.modem->demodulate(b.sigma, b.noisy_symbols, b.LLRs);
			m.decoder->decode_siho(b.LLRs, b.dec_bits);
			//m.decoder_std->decode_siho(b.LLRs, b.dec_bits);
			m.monitor->check_errors(b.dec_bits, b.ref_bits);
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
	//p.info_bits_pos.reserve(p.K);
	std::cout << "# * Simulation parameters: " << std::endl;
	std::cout << "#    ** Correction power = " << p.t << std::endl;
	std::cout << "#    ** Frame errors   = " << p.fe << std::endl;
	std::cout << "#    ** Noise seed     = " << p.seed << std::endl;
	std::cout << "#    ** Info. bits (K) = " << p.K << std::endl;
	std::cout << "#    ** Frame size (N) = " << p.N << std::endl;
	std::cout << "#    ** Code rate  (R) = " << p.R << std::endl;
	std::cout << "#    ** SNR min   (dB) = " << p.ebn0_min << std::endl;
	std::cout << "#    ** SNR max   (dB) = " << p.ebn0_max << std::endl;
	std::cout << "#    ** SNR step  (dB) = " << p.ebn0_step << std::endl;
	std::cout << "#" << std::endl;
	/*std::ifstream alist_file;
	alist_file.open("E:/1111/filename.alist", std::ifstream::in);
	p.H = tools::AList::read(alist_file);*/

	for (unsigned j = 0; j < p.K; ++j) {
		p.info_bits_pos.push_back(j);
	}
}

void init_modules(const params& p, modules& m)
{
	auto h_mat = get_h_alist(LDPC_H_MATRIX_NAME);
	m.source = std::unique_ptr<module::Source_random        <>>(new module::Source_random         <>(p.K));
	//m.encoder = std::unique_ptr<module::Encoder_repetition_sys<>>(new module::Encoder_repetition_sys<>(p.K, p.N));
	//m.encoder = std::unique_ptr<module::Encoder_RS        <>>(new module::Encoder_RS       <>(p.K, p.N, p.Pol_Gen));
	m.encoder = std::unique_ptr<module::Encoder_LDPC_from_H  <>>(new module::Encoder_LDPC_from_H       <>(p.K, p.N, h_mat));
	//m.encoder = std::unique_ptr<module::Encoder_LDPC_DVBS2  <>>(new module::Encoder_LDPC_DVBS2       <>(*p.DVBS2));
	m.modem = std::unique_ptr<module::Modem_BPSK            <>>(new module::Modem_BPSK            <>(p.N));
	m.channel = std::unique_ptr<module::Channel_AWGN_LLR    <>>(new module::Channel_AWGN_LLR      <>(p.N));
	//m.decoder = std::unique_ptr<module::Decoder_repetition_std<>>(new module::Decoder_repetition_std<>(p.K, p.N));
	//m.decoder_std = std::unique_ptr<module::Decoder_RS_std  <>>(new module::Decoder_RS_std   <>(p.K, p.N, p.Pol_Gen));
	//m.decoder = std::unique_ptr<module::Decoder_LDPC_bit_flipping_OMWBF<>>(new module::Decoder_LDPC_bit_flipping_OMWBF<>(p.K, p.N, 30, p.H, p.info_bits_pos));
	m.decoder = std::unique_ptr<module::Decoder_LDPC_BP_flooding<>>(new module::Decoder_LDPC_BP_flooding<>(p.K, p.N, p.n_ite, p.H, p.info_bits_pos, p.update_rule));
	m.monitor = std::unique_ptr<module::Monitor_BFER        <>>(new module::Monitor_BFER          <>(p.K, p.fe));
	m.channel->set_seed(p.seed);
};

void init_buffers(const params& p, buffers& b)
{
	b.ref_bits = std::vector<int  >(p.K);
	b.enc_bits = std::vector<int  >(p.N);
	b.symbols = std::vector<float>(p.N);
	b.sigma = std::vector<float>(1);
	b.noisy_symbols = std::vector<float>(p.N);
	b.LLRs = std::vector<float>(p.N);
	b.dec_bits = std::vector<int  >(p.K);
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