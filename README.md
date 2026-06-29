## Содержание
1. [Результат работы системы](#1-результат-работы-системы)
2. [Содержание веток репозитория](#2-содержание-веток-репозитория)
3. [Сборка библиотеки aff3ct](#3-сборка-библиотеки-aff3ct)
4. [Иерархическая структура проекта](#4-иерархическая-структура-проекта)
5. [Источники](#5-источники)

## 1. Результат работы системы

Идея алгоритма списочного декодирования заключается в следующем:

Шаг 1: Запуск BP, получение LLR.

•	Инициализация LLR в переменных узлах (матрица значений LLR, повторяющая структуру матрицы H).
•	Циклический обмен сообщениями между Variable и Check узлами (горизонтальный шаг: обработка сообщения (V2C) в области вероятности для перехода от LLR к вероятностям. Исключается i-ый узел переменных из рассмотрения):

1.	Variable → Check: Переменные узлы отправляют начальные сообщения в проверяемые узлы. Отправка начальных сообщений. 
2.	Check → Variable: Проверяемые узлы вычисляют и отправляют сообщения обратно в переменные узлы. Обновление сообщений.
3.	Variable → Check: Переменные узлы обновляются и отправляют новые сообщения. Обновление LLR
В ходе итераций SPA проверочные узлы отправляют переменным узлам (VN) "корректирующие" LLR, основанные на проверочных уравнениях. LLR для каждого бита — это сумма:

1.	Начального LLR из канала (прямая информация о бите).
2.	Сообщений от всех проверочных узлов (корректирующая информация).

•	Проверка, удовлетворяет ли текущее декодированное слово всем проверочным уравнениям LDPC-кода (после каждой итерации сообщений априорные вероятности становятся апостериорными):

1.	Получение текущих оценок битов (LLR > 0, то бит = 0; LLR < 0, то бит = 1; если LLR = 0, условие сразу возвращает false (неопределённость)).
2.	Проверка чётности для всех проверочных узлов (если сумма по модулю 2 нечётна, проверка не пройдена).
	То есть если условие выполняется, алгоритм останавливается досрочно. По итогу итеративного декодирования LDPC возможны два исхода: либо синдром зануляется, либо проходит заданное количество итераций. В первом случае первое слово в списке автоматически проходит проверку на синдром, а все остальные — с меньшей вероятностью, так как в них псевдослучайно перевёрнуты биты. Во втором случае первое слово в списке точно не может пройти синдромную проверку, а все остальные - маловероятно. 

Шаг 2: если BP не сошелся:

•	Поиск наименее надежных битов
•	Сортировка по надежности, то есть по возрастанию LLR (чем ближе к 0, тем менее надежен бит)
•	Создание базового кандидата (первый кандидат — текущая оценка битов)
•	Генерация списка кандидатов, путём переворота (инвертирования) наименее надёжных битов. Каждая инверсия создает нового кандидата.

Когда все биты обработаны производится выбор лучшего кандидата:

•	Критерий 1: Выбор первого валидного (удовлетворяющего проверкам). Временная подстановка кандидата в переменный узел для проверки (на основе минимального синдрома). Если таких нет, то:
•	Критерий 2: Выбор по максимальной сумме |LLR| (наиболее надежный кандидат)

Итог: списочный декодер возвращает оптимальный декодированный вектор.

Для просмотра графиков (.grf) используется программа Graph 4.4.2


## 2. Содержание веток репозитория

| ОК | Ветки | Версия библиотеки |
|:------------------------:|------------------------|:----------------------:|
|:construction:|**4th-UNN__Folded-Reed-Solomon-(old)** — для пересборки библиотеки.|(aff3ct-3.0.2)|
|:x:|**4th-UNN__Folded-Reed-Solomon** — .lib + сборка с помощью кроссплатформенной системы CMake.|(aff3ct-3.0.2)|
|:construction:|**5th-UNN__LDPC-(old)** — для пересборки библиотеки.|(aff3ct-3.0.2)|
|:x:|**5th-UNN__LDPC** — .lib + сборка с помощью кроссплатформенной системы CMake.|(aff3ct-3.0.2)|
|:white_check_mark:|**6th-UNN__Polar** — просто собрать static lib (.lib).|(aff3ct-3.0.2)|
|:white_check_mark:|**6th-UNN__conv_RSС** — просто собрать static lib (.lib).|(aff3ct-4.4.0)|


## 3. Сборка библиотеки aff3ct

### Необходимые компоненты

- [CMake](https://cmake.org/)
- [Git](https://git-scm.com/)
- C++11-compatible compiler (MSVC 2022, GCC, Clang)

### 1. Клонирование репозитория

```bash
git clone https://github.com/aff3ct/my_project_with_aff3ct.git
git clone --recursive https://github.com/aff3ct/aff3ct.git
cd aff3ct
git submodule update --init --recursive
```

### 2. Конфигурирование CMakeLists.txt

### 2.1. Статическая и динамическая библиотека.

```cmake
set(AFF3CT_COMPILE_STATIC_LIB ON)
set(AFF3CT_COMPILE_SHARED_LIB ON)
```

### 2.2 \bigobj для Debug (MSVC)

```cmake
if(MSVC)
    add_compile_options(/bigobj)
else()
    add_compile_options(-Wa,-mbig-obj)
endif()
```

### 2.3 Отключение `addr2line`

Чтобы избежать проблем с разрешением символов в Windows, явно отключите `addr2line`:

```cmake
set(CPPTRACE_GET_SYMBOLS_WITH_ADD2LINE OFF)
```

<img width="660" alt="1" src="https://github.com/user-attachments/assets/70be9ff0-c73d-421a-9f42-87cc0ef5188c" />

### 3. Настройка JSON

| Переменная | Значение |
|---|---|
| `AFF3CT_COMPILE_SHARED_LIB` | ✅ |
| `AFF3CT_COMPILE_STATIC_LIB` | ✅ |
| `CPPTRACE_GET_SYMBOLS_WITH_ADD2LINE` | ❌ |
| `SPU_TESTS` | ❌ |
| Генератор CMake | Visual Studio 17 2022 (x64) |


### 4. Выбрать цель для сборки собрать и извлечь из out/install/ содержимое include, .lib, .dll

```
install/
├── bin/
├── include/
├── lib/
└── share/
```

### 5. Сборка 
```
cd /D E:\1-ННГУ\6-ВКР\2-Практика\2_LDPC_old

cmake -B build -G "Visual Studio 17 2022" -A x64 ^
  -DCMAKE_PREFIX_PATH="E:/1-ННГУ/6-ВКР/2-Практика/aff3ct-3.0.2/lib/cmake/aff3ct-3.0.2" ^
  -DCMAKE_CXX_FLAGS="/EHsc /MP4 /D_CRT_SECURE_NO_DEPRECATE"

cmake --build build --config Debug --target frs1
```

## 4. Иерархическая структура проекта

```
6th-UNN__List-Decoding-with-Aff3ct/ (для всех веток)
│
├── List_Decoding_with_Aff3ct.sln
├── 1_RS_old
│   ├── 0-RS_old.vcxproj
│   ├── 0-RS_old.vcxproj.filters
│   ├── 0-RS_old.vcxproj.user
│   ├── 1-RS.vcxproj
│   ├── 1-RS.vcxproj.filters
│   ├── 1-RS.vcxproj.user
│   ├── 123.bat
│   ├── CMakeLists.txt
│   ├── Decoder_FRS1.cpp
│   ├── Decoder_FRS1.hpp
│   ├── Decoder_FRS2.cpp
│   ├── Decoder_FRS2.hpp
│   ├── Encoder_FRS.cpp
│   ├── Encoder_FRS.hpp
│   ├── FRS_List_Selector.hpp
│   ├── gdiplus.h
│   ├── main_frs1.cpp
│   ├── main_frs2.cpp
│   ├── win_guard.h
│   ├── RS_Result_old
│   │   ├── FRS 0.grf
│   │   ├── FRS 1.grf
│   │   ├── FRS FInal.grf
│   │   ├── FRS Theory.grf
│   │   ├── Interleaver.grf
│   │   ├── Interleaver2.grf
│   │   ├── Packege Error.grf
│   │   └── RS + BPSK2.grf
│   ├── RS_Tests_old
│   │   ├── m=1 (2).txt
│   │   ├── m=1.txt
│   │   ├── m=2 (2).txt
│   │   ├── m=2.txt
│   │   ├── m=3 (2).txt
│   │   ├── m=3.txt
│   │   ├── m=4 (2).txt
│   │   ├── m=4.txt
│   │   ├── Package Error.txt
│   │   ├── QAM+Interveaver_test.txt
│   │   └── RS_test.txt
│   └── src_old
│       ├── FRS_Decoder_old.cpp
│       ├── FRS_Decoder_old.hpp
│       ├── FRS_old.cpp
│       └── RS_old.cpp
├── 2_LDPC_old
│   ├── 0-LDPC_old.vcxproj
│   ├── 0-LDPC_old.vcxproj.filters
│   ├── 0-LDPC_old.vcxproj.user
│   ├── 2-LDPC.vcxproj
│   ├── 2-LDPC.vcxproj.filters
│   ├── 2-LDPC.vcxproj.user
│   ├── LDPC_Result_old
│   │   ├── + (N=1008 ,K=504,M=504,R= 0.5).gz
│   │   ├── + (N=271,K=144,M=127,R= 0.5314).txt
│   │   ├── + (N=4095,K=3358,M=737,R= 0.82).txt
│   │   ├── + (N=4376,K=4094,M=282,R= 0.9356).gz
│   │   ├── + (N=495,K=433,M=62,R= 0.8747).gz
│   │   ├── + (N=504 ,K=252,M=252,R= 0.5).gz
│   │   ├── + (N=9972 ,K=4986,M=4986,R= 0.5).gz
│   │   ├── 123.txt
│   │   ├── aff3ct-3.0.2.exe
│   │   ├── filename.alist
│   │   └── LDPC + RS.grf
│   ├── LDPC_Tests_old
│   │   ├── LDPC_test.txt
│   │   └── RS.txt
│   └── src_old
│       ├── Decoder_LDPC_old.cpp
│       ├── Decoder_LDPC_old.hpp
│       ├── LDPC_List_Decoder_old.cpp
│       └── LDPC_old.cpp
├── 3_Polar
│   ├── 3-Polar.vcxproj
│   ├── 3-Polar.vcxproj.filters
│   ├── 3-Polar.vcxproj.user
│   ├── codec_polar.hpp
│   ├── Polar CRC.cpp
│   ├── Polar.cpp
│   ├── Polar_Result
│   │   ├── CMakeLists.txt
│   │   ├── Deepseek.txt
│   │   ├── LDPC + RS - Only.grf
│   │   ├── LDPC + RS - R.grf
│   │   └── LDPC + RS — Full.grf
│   └── Polar_Tests
│       ├── 123.txt
│       ├── Decoder_polar_ASCL_fast_CA_sys.txt
│       ├── Decoder_polar_SCL_fast_CA_sys.txt
│       ├── Decoder_polar_SCL_fast_sys.txt
│       ├── Decoder_polar_SCL_MEM_fast_sys.txt
│       ├── Decoder_polar_SC_fast_sys.txt
│       ├── Frozen.txt
│       └── Template.txt
├── 4_RSС
│   ├── 4-RSC.vcxproj
│   ├── 4-RSC.vcxproj.filters
│   ├── 4-RSC.vcxproj.user
│   ├── Conv.cpp
│   ├── Conv_list.cpp
│   ├── RSC.cpp
│   ├── RSC_list.cpp
│   ├── conv_RSC_Result
│   │   ├── conv + RSC - full.grf
│   │   └── conv + RSC.grf
│   ├── conv_Tests
│   │   ├── Viterbi
│   │   │   ├── 13-15.txt
│   │   │   ├── 133-171.txt
│   │   │   └── 5-7.txt
│   │   └── Viterbi_Lists
│   │       ├── 11-13.txt
│   │       ├── 13-11.txt
│   │       ├── 13-15.txt
│   │       ├── 13-17.txt
│   │       ├── 133-171.txt
│   │       ├── 15-17.txt
│   │       └── 5-7.txt
│   └── RSC_Tests
│       ├── Viterbi
│       │   ├── 11-13.txt
│       │   ├── 13-11.txt
│       │   ├── 13-15.txt
│       │   ├── 13-17.txt
│       │   ├── 133-171.txt
│       │   ├── 15-17.txt
│       │   └── 5-7.txt
│       └── Viterbi_Lists
│           ├── 13-15.txt
│           ├── 133-171.txt
│           └── 5-7.txt
├── aff3ct-3.0.2
└── aff3ct-4.4.0
```

## 5. Источники
- Aff3ct - https://aff3ct.github.io/
- Aff3ct Documentation - https://aff3ct.readthedocs.io/en/latest/
- Using AFF3CT as a library for your codes - https://github.com/aff3ct/my_project_with_aff3ct
- Graph - https://www.padowan.dk/
- Essential Coding Theory - https://cse.buffalo.edu/faculty/atri/courses/coding-theory/book/
