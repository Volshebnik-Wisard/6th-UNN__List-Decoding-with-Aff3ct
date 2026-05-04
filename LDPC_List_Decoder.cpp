#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <algorithm>
#include <stdexcept>

using namespace std;


// Представляет собой контрольный узел в графе Таннера кода LDPC
class CheckNode {
private:
	map<int, double> receivedMessage; // Карта, хранящая сообщения из подключенных узлов переменных

public:
	double calcMessage(int to) // Вычисляет исходящее сообщение с использованием операций tanh/atanh
	{
		double product = 1.0;
		for (const auto& [index, message] : receivedMessage) {
			if (index == to) continue;
			if (isinf(message) {
				if (message < 0) product *= -1;
				continue;
			}
			product *= tanh(message / 2);
		}
		return 2 * atanh(product);
	}

	// Получение сообщения (Сохраняет входящее сообщение (узел))
	void receiveMessage(int from, double message) {
		receivedMessage[from] = message;
	}

	// Очистить все полученные сообщения
	void clear() {
		receivedMessage.clear();
	}
};

// Представляет собой соединение между переменным и контрольным узлами
struct Edge {
	int variableNodeIndex; // Индекс узла переменной, соединенного этим ребром.
	int checkNodeIndex; // Индекс узла проверки (контрольные), соединенного этим ребром.
};

// представляет собой узел в графе Таннера LDPC
class VariableNode {
private:
	double channelLLR = 0.0; // Начальное отношение правдоподобия логарифма из канала
	map<int, double> receivedMessage; // Карта, хранящая сообщения от подключенных контрольных узлов (переменных узлов)
	bool isFrozen = false; // Флаг, указывающий, заморожен ли бит (исправлен)

public:
	//-Замороженные биты : Возвращает `+Inf` (устанавливает бит `0`).
	//- Незамороженные: Возвращает `ChannelLLR`.

	// Устанавливает замороженный статус (как в полярных кодах)
	void setIsFrozen(bool frozen) {
		isFrozen = frozen;
	}

	// Возвращает начальный LLR (Inf для замороженных битов)
	double calcInitialMessage() {
		if (isFrozen) return INFINITY;
		return channelLLR;
	}

	// Вычисляет исходящее сообщение путем суммирования значений LLR
	/*
	1. Суммирует `ChannelLLR` и все входящие LLR (исключая целевой узел).
    2. Если какое-либо входящее сообщение равно `±Inf`, распространяет его напрямую (раннее завершение).
       Замороженные биты: Всегда возвращает `+Inf` (переопределяет вычисления).
     */
	double calcMessage(int to) {
		if (isFrozen) return INFINITY;

		double sum = channelLLR;
		for (const auto& [index, message] : receivedMessage) {
			if (index == to) continue;
			if (isinf(message)) return message;
			sum += message;
		}
		return sum;
	}

	// Сохраняет входящее сообщение
	void receiveMessage(int from, double message) {
		receivedMessage[from] = message;
	}

	// Вычисляет окончательное LLR для оценки битов путем суммирования всех входящих сообщений и `ChannelLLR`.
	double marginalize() {
		return calcMessage(-1); // -1 указывает на отсутствие целевого узла
	}

	// Оценивает значение бита на основе LLR (Преобразует значение LLR в битовое значение (0/1))
	// Жесткое решение (`0` или `1`)
	int estimateSendBit() {
		double llr = marginalize();
		if (llr > 0) return 0; //предпочтение биту `0`
		if (llr < 0) return 1; //предпочтение биту `1`

		// Случайное значение при нулевом LLR
		static random_device rd;
		static mt19937 gen(rd());
		uniform_int_distribution<> dis(0, 1);
		return dis(gen);
	}

	// Очистить состояние узла
	void clear() {
		channelLLR = 0.0;
		receivedMessage.clear();
	}

	// Установить канал LLR
	void setChannelLLR(double llr) {
		channelLLR = llr;
	}
};

// Код LDPC представляет собой полную структуру кода LDPC
class LDPCCode {
private:
	int codeLength; // Общая длина кодового слова (количество переменных узлов)
	vector<Edge> edges; // Связи между узлами (Список всех ребер, определяющих связи между переменными и контрольными узлами.)
	vector<int> informationBitIndexes; // Индексы информационных битов (переменных узлов, переносящих информационные биты (незамороженные))
	vector<int> frozenBitIndexes; // Индексы замороженных битов
	vector<VariableNode> variableNodes; // Массив переменных узлов
	vector<CheckNode> checkNodes; // Массив контрольных узлов
	const int decodeIteration = 40; // Кол-во итераций

	// Проверяет, выполнены ли все проверки на четность (удовлетворяет ли текущее кодовое слово всем ограничениям четности)
	/*
	1. Оценивает биты для всех узлов переменных с помощью `Marginalize()`.
    2. Суммирует биты, подключенные к каждому контрольному узлу (по модулю 2).
    3. Возвращает `true`, только если все суммы равны `0` (четность).
	*/
	bool isSatisfyAllChecks() {
		vector<int> estimates(codeLength);
		for (int i = 0; i < codeLength; i++) {
			double llr = variableNodes[i].marginalize();
			if (llr == 0) return false;
			estimates[i] = (llr < 0) ? 1 : 0;
		}

		vector<int> checks(checkNodes.size(), 0);
		for (const auto& edge : edges) {
			checks[edge.checkNodeIndex] += estimates[edge.variableNodeIndex];
		}

		for (int c : checks) {
			if (c % 2 != 0) return false;
		}
		return true;
	}

	// Запускает алгоритм передачи сообщений
	void executeMessagePassing(const vector<double>& channelOutputs) {
		// Инициализируйте переменные узлы с помощью `ChannelLLR`
		for (int i = 0; i < codeLength; i++) {
			variableNodes[i].setChannelLLR(channelOutputs[i]);
		}

		// Первая половина итерации: узлы переменные отправляют начальные сообщения узлам проверки.
		for (int i = 0; i < edges.size(); i++) {
			const auto& edge = edges[i];
			double message = variableNodes[edge.variableNodeIndex].calcInitialMessage();
			checkNodes[edge.checkNodeIndex].receiveMessage(i, message);
		}

		// Основной цикл передачи сообщений
		for (int iter = 0; iter < decodeIteration; iter++) {
			// Узлы проверки вычисляют и отправляют сообщения обратно узлам переменной.
			for (int i = 0; i < edges.size(); i++) {
				const auto& edge = edges[i];
				double message = checkNodes[edge.checkNodeIndex].calcMessage(i);
				variableNodes[edge.variableNodeIndex].receiveMessage(i, message);
			}

			// Узлы переменной обновляют и отправляют новые сообщения.
			for (int i = 0; i < edges.size(); i++) {
				const auto& edge = edges[i];
				double message = variableNodes[edge.variableNodeIndex].calcMessage(i);
				checkNodes[edge.checkNodeIndex].receiveMessage(i, message);
			}

			// Досрочное завершение, если все проверки будут выполнены
			if (isSatisfyAllChecks()) break;
		}
	}

public:
	// Декодирование выходных данных принимаемого канала (Извлекает декодированные биты (используя `EstimateSendBit()`) только для информационных битов.)
	vector<int> decode(const vector<double>& channelOutputs) {
		executeMessagePassing(channelOutputs);
		vector<int> decoded;
		for (int index : informationBitIndexes) {
			decoded.push_back(variableNodes[index].estimateSendBit());
		}
		return decoded;
	}

	// Декодирование списка для неоднозначных битов (Генерирует несколько кодовых слов-кандидатов для неоднозначных битов.)
	// Ранжирует ненадежные биты по величине LLR и переворачивает их для изучения альтернатив.
	vector<vector<int>> listDecode(const vector<double>& channelOutputs, int listSize) {
		// Определяет наименее надежные биты (ближайшие к `0` LLR) среди информационных битов.
		executeMessagePassing(channelOutputs);
		// Генерирует кандидатов `listSize`, комбинаторно переворачивая эти биты.
		int ambiguousBitCount = floor(log2(listSize));
		vector<pair<int, double>> llrs;
		for (int index : informationBitIndexes) {
			double llr = variableNodes[index].marginalize();
			llrs.emplace_back(index, llr);
		}
		// Сортировка по величине LLR (сначала наименее надежный)
		sort(llrs.begin(), llrs.end(), [](const auto& a, const auto& b) {
			return abs(a.second) < abs(b.second);
			});

		// Получение уникального декодированного вектора
		vector<int> uniqueDecoded;
		for (int index : informationBitIndexes) {
			uniqueDecoded.push_back(variableNodes[index].estimateSendBit());
		}

		vector<vector<int>> listDecoded = { uniqueDecoded };

		// Возвращает список кодовых слов, включая исходные и перевернутые варианты.
		for (int i = 0; i < ambiguousBitCount; i++) {
			auto [index, _] = llrs[i];
			vector<vector<int>> temp;
			for (const auto& v : listDecoded) {
				vector<int> inverted = v;
				inverted[index] = 1 - inverted[index];
				temp.push_back(inverted);
			}
			listDecoded.insert(listDecoded.end(), temp.begin(), temp.end());
		}

		return listDecoded;
	}

	// Возвращает скорость кодового слова
	double getRate() const {
		return static_cast<double>(informationBitIndexes.size()) / getRealCodeLength();
	}

	// Возвращает эффективную скорость для декодирования списка
	double getListRate(int listSize) const {
		int ambiguousBitCount = floor(log2(listSize));
		return static_cast<double>(informationBitIndexes.size() - ambiguousBitCount) / getRealCodeLength();
	}

	// Возвращает эффективную длину кодового слова
	int getRealCodeLength() const {
		return codeLength - frozenBitIndexes.size();
	}

	// Создание случайныго обычного LDPC-кода
		static LDPCCode constructCode(
		int originalCodeLength,
		int informationBitSize,
		int variableNodeDegree,
		int checkNodeDegree) {

		LDPCCode code;
		code.codeLength = originalCodeLength;

		// Устанавливает информационные и замороженные битовые индексы
		for (int i = 0; i < informationBitSize; i++) {
			code.informationBitIndexes.push_back(i);
		}

		int originalInfoBitSize = originalCodeLength -
			(originalCodeLength * variableNodeDegree) / checkNodeDegree;

		for (int i = informationBitSize; i < originalInfoBitSize; i++) {
			code.frozenBitIndexes.push_back(i);
		}

		// Создание ребер
		code.edges = createRandomEdges(originalCodeLength, variableNodeDegree, checkNodeDegree);

		// Инициализировать узлы
		code.variableNodes.resize(originalCodeLength);
		code.checkNodes.resize((originalCodeLength * variableNodeDegree) / checkNodeDegree);

		// Набор замороженных бит
		for (int index : code.frozenBitIndexes) {
			code.variableNodes[index].setIsFrozen(true);
		}

		return code;
	}

private:
	// Создаёт случайные ребра для графа Таннера
	// граф Таннера является случайным, но регулярным (фиксированные степени).
	static vector<Edge> createRandomEdges(int length, int variableNodeDegree, int checkNodeDegree) {
		if (length * variableNodeDegree % checkNodeDegree != 0) {
			throw invalid_argument("Invalid length and degree combination");
		}

		random_device rd;
		mt19937 gen(rd());

		vector<Edge> edges(length * variableNodeDegree);
		vector<int> temp(length * variableNodeDegree);

		// Initialize with variable node indices
		for (int i = 0; i < temp.size(); i++) {
			temp[i] = i / variableNodeDegree;
		}

		// Shuffle the connections
		shuffle(temp.begin(), temp.end(), gen);

		// Create edges
		for (int k = 0; k < temp.size(); k++) {
			edges[k] = Edge{ temp[k], k / checkNodeDegree };
		}

		return edges;
	}
};