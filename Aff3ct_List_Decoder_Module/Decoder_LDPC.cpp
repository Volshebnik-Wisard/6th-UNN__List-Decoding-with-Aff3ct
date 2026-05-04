#include "Decoder_LDPC.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B, typename R>
double Decoder_LDPC<B,R>::CheckNode::calcMessage(int to)
{
    double product = 1.0;
    for (const auto& [index, message] : receivedMessage) {
        if (index == to) continue;
        if (std::isinf(message)) {
            if (message < 0) product *= -1;
            continue;
        }
        product *= std::tanh(message / 2);
    }
    return 2 * std::atanh(product);
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::CheckNode::receiveMessage(int from, double message)
{
    receivedMessage[from] = message;
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::CheckNode::clear()
{
    receivedMessage.clear();
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::VariableNode::setIsFrozen(bool frozen)
{
    isFrozen = frozen;
}

template <typename B, typename R>
double Decoder_LDPC<B,R>::VariableNode::calcInitialMessage()
{
    if (isFrozen) return std::numeric_limits<double>::infinity();
    return channelLLR;
}

template <typename B, typename R>
double Decoder_LDPC<B,R>::VariableNode::calcMessage(int to)
{
    if (isFrozen) return std::numeric_limits<double>::infinity();

    double sum = channelLLR;
    for (const auto& [index, message] : receivedMessage) {
        if (index == to) continue;
        if (std::isinf(message)) return message;
        sum += message;
    }
    return sum;
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::VariableNode::receiveMessage(int from, double message)
{
    receivedMessage[from] = message;
}

template <typename B, typename R>
double Decoder_LDPC<B,R>::VariableNode::marginalize()
{
    return calcMessage(-1);
}

template <typename B, typename R>
int Decoder_LDPC<B,R>::VariableNode::estimateSendBit()
{
    double llr = marginalize();
    if (llr > 0) return 0;
    if (llr < 0) return 1;

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    return (int)dis(gen);
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::VariableNode::clear()
{
    channelLLR = 0.0;
    receivedMessage.clear();
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::VariableNode::setChannelLLR(double llr)
{
    channelLLR = llr;
}

template <typename B, typename R>
Decoder_LDPC<B,R>::Decoder_LDPC(const int K, const int N, const std::vector<int>& frozen_bits)
{
    const std::string name = "Decoder_LDPC";
    this->set_name(name);

    // Инициализировать замороженные биты
    frozenBitIndexes = frozen_bits;
    informationBitIndexes.reserve(K);
    for (int i = 0; i < N; i++) {
        if (std::find(frozen_bits.begin(), frozen_bits.end(), i) == frozen_bits.end()) {
            informationBitIndexes.push_back(i);
        }
    }

    // Инициализация узлов
    variableNodes.resize(N);
    for (int index : frozenBitIndexes) {
        variableNodes[index].setIsFrozen(true);
    }
}

template <typename B, typename R>
Decoder_LDPC<B,R>* Decoder_LDPC<B,R>::clone() const
{
    auto m = new Decoder_LDPC(*this);
    m->deep_copy(*this);
    return m;
}

template <typename B, typename R>
bool Decoder_LDPC<B,R>::isSatisfyAllChecks()
{
    std::vector<B> estimates(codeLength);
    for (int i = 0; i < codeLength; i++) {
        double llr = variableNodes[i].marginalize();
        if (llr == 0) return false;
        estimates[i] = (llr < 0) ? 1 : 0;
    }

    std::vector<int> checks(checkNodes.size(), 0);
    for (const auto& edge : edges) {
        checks[edge.checkNodeIndex] += estimates[edge.variableNodeIndex];
    }

    for (int c : checks) {
        if (c % 2 != 0) return false;
    }
    return true;
}

template <typename B, typename R>
void Decoder_LDPC<B,R>::executeMessagePassing(const std::vector<double>& channelOutputs)
{
    // Инициализация переменных узлов с помощью chanelLLR 
    for (int i = 0; i < codeLength; i++) {
        variableNodes[i].setChannelLLR(channelOutputs[i]);
    }

    // Первая половина итерации: переменные узлы отправляют начальные сообщения для проверки узлов
    for (int i = 0; i < (int)edges.size(); i++) {
        const auto& edge = edges[i];
        double message = variableNodes[edge.variableNodeIndex].calcInitialMessage();
        checkNodes[edge.checkNodeIndex].receiveMessage(i, message);
    }

    // Основной цикл передачи сообщений
    for (int iter = 0; iter < decodeIteration; iter++) {
        // Проверяемые узлы вычисляют и отправляют сообщения обратно в переменные узлы
        for (int i = 0; i < (int)edges.size(); i++) {
            const auto& edge = edges[i];
            double message = checkNodes[edge.checkNodeIndex].calcMessage(i);
            variableNodes[edge.variableNodeIndex].receiveMessage(i, message);
        }

        // Переменные узлы обновляются и отправляют новые сообщения
        for (int i = 0; i < (int)edges.size(); i++) {
            const auto& edge = edges[i];
            double message = variableNodes[edge.variableNodeIndex].calcMessage(i);
            checkNodes[edge.checkNodeIndex].receiveMessage(i, message);
        }

        // Досрочное расторжение, если все проверки удовлетворены
        if (isSatisfyAllChecks()) break;
    }
}

template <typename B, typename R>
std::vector<int> Decoder_LDPC<B,R>::decode(const std::vector<double>& channelOutputs)
{
    executeMessagePassing(channelOutputs);
    std::vector<int> decoded;
    for (int index : informationBitIndexes) {
        decoded.push_back(variableNodes[index].estimateSendBit());
    }
    return decoded;
}

template <typename B, typename R>
std::vector<int> Decoder_LDPC<B, R>::selectBestCandidate(
    const std::vector<std::vector<int>>& candidates,
    const std::vector<double>& channelOutputs)
{
    if (candidates.empty()) return {};

    // Критерий 1: Выбор первого валидного (удовлетворяющего проверкам)
    for (const auto& candidate : candidates) {
        // Временная подстановка candidate в variableNodes для проверки isSatisfyAllChecks()
        for (size_t i = 0; i < informationBitIndexes.size(); i++) {
            int index = informationBitIndexes[i];
            variableNodes[index].setChannelLLR(candidate[i] ? -1.0 : 1.0); // Эмулируем LLR для бита
        }
        if (isSatisfyAllChecks()) {
            return candidate;
        }
    }

    // Критерий 2: Выбор по максимальной сумме |LLR| (наиболее надежный)
    std::vector<int> bestCandidate;
    double maxLLRSum = -1.0;
    for (const auto& candidate : candidates) {
        double llrSum = 0.0;
        for (size_t i = 0; i < candidate.size(); i++) {
            int index = informationBitIndexes[i];
            llrSum += std::abs(variableNodes[index].marginalize());
        }
        if (llrSum > maxLLRSum) {
            maxLLRSum = llrSum;
            bestCandidate = candidate;
        }
    }
    return bestCandidate;
}

template <typename B, typename R>
std::vector<std::vector<int>> Decoder_LDPC<B,R>::listDecode(const std::vector<double>& channelOutputs, int listSize)
{
    executeMessagePassing(channelOutputs);

    // Найти наименее надежные биты
    std::vector<std::pair<int, double>> llrs;
    for (int index : informationBitIndexes) {
        double llr = variableNodes[index].marginalize();
        llrs.emplace_back(index, llr);
    }

    // Сортировать по надежности (сначала ближе к 0)
    std::sort(llrs.begin(), llrs.end(), [](const auto& a, const auto& b) {
        return std::abs(a.second) < std::abs(b.second);
        });

    // Генерация уникального декодированного вектора
    std::vector<int> uniqueDecoded;
    for (int index : informationBitIndexes) {
        uniqueDecoded.push_back(variableNodes[index].estimateSendBit());
    }

    std::vector<std::vector<int>> listDecoded = { uniqueDecoded };

    // Сгенерировать список кандидатов, перевернув наименее надежные биты
    int ambiguousBitCount = std::floor(std::log2(listSize));
    for (int i = 0; i < ambiguousBitCount; i++) {
        auto [index, _] = llrs[i];
        std::vector<std::vector<int>> temp;
        for (const auto& v : listDecoded) {
            std::vector<int> inverted = v;
            inverted[index] = 1 - inverted[index];
            temp.push_back(inverted);
        }
        listDecoded.insert(listDecoded.end(), temp.begin(), temp.end());
    }

    // Выбор лучшего кандидата
    return selectBestCandidate(listDecoded, channelOutputs);
}

template <typename B, typename R>
int Decoder_LDPC<B,R>::_decode_siso(const R *Y_N1, R *Y_N2, const size_t frame_id)
{
    // В этом примере не реализовано для LDPC
    return 0;
}

template <typename B, typename R>
int Decoder_LDPC<B,R>::_decode_siho(const R *Y_N, B *V_K, const size_t frame_id)
{
    std::vector<R> Y_N_vec(Y_N, Y_N + this->N);
    auto decoded = decode(Y_N_vec);
    std::copy(decoded.begin(), decoded.end(), V_K);
    return 0;
}

template <typename B, typename R>
int Decoder_LDPC<B,R>::_decode_siho_cw(const R *Y_N, B *V_N, const size_t frame_id)
{
    std::vector<R> Y_N_vec(Y_N, Y_N + this->N);
    executeMessagePassing(Y_N_vec);
    for (int i = 0; i < this->N; i++) {
        V_N[i] = variableNodes[i].estimateSendBit();
    }
    return 0;
}

template <typename B, typename R>
double Decoder_LDPC<B,R>::getRate() const
{
    return static_cast<double>(informationBitIndexes.size()) / getRealCodeLength();
}

template <typename B, typename R>
double Decoder_LDPC<B,R>::getListRate(int listSize) const
{
    int ambiguousBitCount = std::floor(std::log2(listSize));
    return static_cast<double>(informationBitIndexes.size() - ambiguousBitCount) / getRealCodeLength();
}

template <typename B, typename R>
int Decoder_LDPC<B,R>::getRealCodeLength() const
{
    return codeLength - frozenBitIndexes.size();
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Decoder_LDPC<B_8, Q_8>;
template class aff3ct::module::Decoder_LDPC<B_16, Q_16>;
template class aff3ct::module::Decoder_LDPC<B_32, Q_32>;
template class aff3ct::module::Decoder_LDPC<B_64, Q_64>;
#else
template class aff3ct::module::Decoder_LDPC<B, Q>;
#endif
// ==================================================================================== explicit template instantiation