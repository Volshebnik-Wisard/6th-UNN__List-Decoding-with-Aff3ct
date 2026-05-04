#ifndef DECODER_LDPC_HPP
#define DECODER_LDPC_HPP

#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>

#include "Module/Decoder/Decoder_SISO.hpp"

namespace aff3ct
{
	namespace module
	{

		template <typename B = int, typename R = float>
		class Decoder_LDPC : public Decoder_SISO<B, R>
		{
		private:
			struct Edge {
				int variableNodeIndex;
				int checkNodeIndex;
			};

			class CheckNode {
			private:
				std::map<int, double> receivedMessage;
			public:
				double calcMessage(int to);
				void receiveMessage(int from, double message);
				void clear();
			};

			class VariableNode {
			private:
				double channelLLR = 0.0;
				std::map<int, double> receivedMessage;
				bool isFrozen = false;
			public:
				void setIsFrozen(bool frozen);
				double calcInitialMessage();
				double calcMessage(int to);
				void receiveMessage(int from, double message);
				double marginalize();
				int estimateSendBit();
				void clear();
				void setChannelLLR(double llr);
			};

			int codeLength;
			std::vector<Edge> edges;
			std::vector<int> informationBitIndexes;
			std::vector<int> frozenBitIndexes;
			std::vector<VariableNode> variableNodes;
			std::vector<CheckNode> checkNodes;
			const int decodeIteration = 40;

			bool isSatisfyAllChecks();
			void executeMessagePassing(const std::vector<double>& channelOutputs);
			std::vector<int> selectBestCandidate(const std::vector<std::vector<int>>& candidates, const std::vector<double>& channelOutputs);

		public:
			Decoder_LDPC(const int K, const int N, const std::vector<int>& frozen_bits);
			virtual ~Decoder_LDPC() = default;
			virtual Decoder_LDPC<B, R>* clone() const;

		protected:
			virtual int _decode_siso(const R* Y_N1, R* Y_N2, const size_t frame_id);
			virtual int _decode_siho(const R* Y_N, B* V_K, const size_t frame_id);
			virtual int _decode_siho_cw(const R* Y_N, B* V_N, const size_t frame_id);

			std::vector<int> decode(const std::vector<double>& channelOutputs);
			std::vector<std::vector<int>> listDecode(const std::vector<double>& channelOutputs, int listSize);
			double getRate() const;
			double getListRate(int listSize) const;
			int getRealCodeLength() const;
		};

	}
}

#endif /* DECODER_LDPC_HPP */