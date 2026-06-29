## Содержание
1. [Результат работы системы](#1-результат-работы-системы)
2. [Содержание веток репозитория](#2-содержание-веток-репозитория)
3. [Сборка библиотеки aff3ct](#3-сборка-библиотеки-aff3ct)
4. [Иерархическая структура проекта](#4-иерархическая-структура-проекта)
5. [Источники](#5-источники)

## 1. Результат работы системы

Для просмотра графиков (.grf) используется программа Graph 4.4.2

**Алгоритм**

<img width="628" height="274" alt="11" src="https://github.com/user-attachments/assets/1fdc9afc-ba78-4402-a1b6-072880217a8f" />


**Компромисс между скоростью кода R и долей ошибок ρ, которые могут быть исправлены алгоритмом 1 для m = 1, 2, 3, 4**

<img width="497" height="365" alt="2" src="https://github.com/user-attachments/assets/657aa264-a580-449f-aaf9-3d5c2be10e38" />


**Графики зависимостей BER от SNR при различных скоростях кода R для BPSK и FRS (красным m=1, синим m=2, зелёным m=3, розовым m=4)**

<img width="1798" height="956" alt="3" src="https://github.com/user-attachments/assets/5fcc4569-464f-4f18-9481-c64405d1edb5" />


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
