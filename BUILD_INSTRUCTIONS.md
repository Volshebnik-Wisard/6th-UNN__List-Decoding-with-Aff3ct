# Сборка и запуск: LDPC List-BP decoder (AFF3CT)

## Структура файлов проекта

Положите все файлы рядом в одну папку:

```
2_LDPC_old/
├── CMakeLists.txt
├── LDPC_list_main.cpp
├── Decoder_LDPC_BP_list.hpp
├── Decoder_LDPC_BP_list.hxx    ← включается автоматически из .hpp
└── filename.alist              ← ваша матрица H
```

---

## Требования

| Компонент | Версия |
|---|---|
| CMake | ≥ 3.1 |
| MSVC (Visual Studio) | 2019 или 2022, x64 |
| AFF3CT | 3.0.2 (уже есть в `E:\1-ННГУ\6-ВКР\2-Практика\aff3ct-3.0.2\`) |

---

## Сборка через CMake (Developer PowerShell / cmd)

### Шаг 1. Перейти в папку проекта

```bat
cd /D E:\1-ННГУ\6-ВКР\2-Практика\2_LDPC_old
```

### Шаг 2. Сконфигурировать

```bat
cmake -B build ^
      -DAFF3CT_DIR="E:/1-ННГУ/6-ВКР/2-Практика/aff3ct-3.0.2/lib/cmake/aff3ct-3.0.2"
```

Если CMake не находит AFF3CT — добавьте генератор явно:
```bat
cmake -B build -G "Visual Studio 17 2022" -A x64 ^
      -DAFF3CT_DIR="E:/1-ННГУ/6-ВКР/2-Практика/aff3ct-3.0.2/lib/cmake/aff3ct-3.0.2"
```

### Шаг 3. Собрать

```bat
cmake --build build --config Release
```

Исполняемый файл: `build\Release\aff3ct-ldpc-list-run.exe`

---

## Запуск

```bat
build\Release\aff3ct-ldpc-list-run.exe  filename.alist  [list_size]
```

| Аргумент | Значение |
|---|---|
| `filename.alist` | Путь к матрице H в формате AList |
| `list_size` | 1 (обычный BP), 2, 4, **8**, 16, 32 (по умолчанию 8) |

Количество кандидатов = 2^floor(log2(list_size)).

### Примеры

```bat
REM Обычный BP:
build\Release\aff3ct-ldpc-list-run.exe  filename.alist  1

REM List-BP, 8 кандидатов (3 неуверенных бита):
build\Release\aff3ct-ldpc-list-run.exe  filename.alist  8

REM List-BP, 32 кандидата (5 неуверенных бит):
build\Release\aff3ct-ldpc-list-run.exe  filename.alist  32
```

---

## Как устроен List-BP декодер

```
Decoder_LDPC_BP_list<B,R>
    └── Decoder_LDPC_BP_flooding_SPA<B,R      (из aff3ct-3.0.2.lib)
            └── Decoder_LDPC_BP_flooding<B,R,SPA>
                    ├── Decoder_SISO<B,R>
                    └── Decoder_LDPC_BP
```

`_decode_siho()` пошагово:

1. **BP-прогон** — родительский `_decode()` (flooding/SPA, ≤10 итераций,
   синдромный останов глубиной 1). Заполняет `this->post[]`.
2. **Ненадёжные биты** — среди K инфо-бит берём `n_amb = floor(log2(list_size))`
   с наименьшим `|post[i]|`.
3. **Перебор** — 2^n_amb паттернов инверсий над этими битами.
4. **Синдром** — `LDPC_syndrome::check_hard()`. Возвращаем первый валидный
   кандидат (минимальная сумма `|post|` по флипнутым позициям).
   Если ни одного — наименее плохой.

---

## Параметры (в LDPC_list_main.cpp)

```cpp
struct params {
    int   K         = 144;      // информационных бит
    int   N         = 271;      // длина кодового слова
    int   n_ite     = 10;       // итераций BP
    int   fe        = 100;      // кадров с ошибками до остановки
    int   seed      = 0;        // зерно ГПСЧ
    float ebn0_min  = 0.00f;
    float ebn0_max  = 15.00f;
    float ebn0_step = 0.50f;
};
```

После изменения — пересобрать: `cmake --build build --config Release`.

---

## Частые ошибки сборки

| Ошибка | Решение |
|---|---|
| `Could not find AFF3CT` | Проверьте путь `-DAFF3CT_DIR=...` |
| `LNK2019: unresolved external` | `aff3ct::aff3ct-static-lib` должен быть в `target_link_libraries` |
| `error C2039: 'post'` | `Decoder_LDPC_BP_list.hxx` должен лежать рядом с `.hpp` |
| `Cannot open alist file` | Передайте полный путь к `.alist`, например `E:\1111\filename.alist` |
