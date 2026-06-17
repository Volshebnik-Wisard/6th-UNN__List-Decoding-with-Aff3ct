cmake_minimum_required(VERSION 3.15)

project(
  aff3ct-frs
  LANGUAGES CXX
  DESCRIPTION "Folded Reed-Solomon codes with AFF3CT"
  VERSION 0.1.0
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

find_package(AFF3CT 3.0.2 REQUIRED
  HINTS "${CMAKE_PREFIX_PATH}"
)

# ── Пути ───────────────────────────────────────────────────────────────────
set(NOWIN_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/nowin")
set(WIN_GUARD  "${CMAKE_CURRENT_SOURCE_DIR}/win_guard.h")

# ── /FI win_guard.h — вставляется ПЕРЕД любым #include в каждом .cpp ──────
# Это гарантирует что windows.h подключится с WIN32_LEAN_AND_MEAN/NOMINMAX
# раньше чем это сделает aff3ct.hpp.
if(MSVC)
  add_compile_options(/FI "${WIN_GUARD}")
endif()

# WIN32_LEAN_AND_MEAN + NOMINMAX глобально
# NOGDI / NOUSER НЕ ставим глобально — rang.hpp нуждается в wincon.h
add_compile_definitions(
  WIN32_LEAN_AND_MEAN
  NOMINMAX
  NOMCX
  NOIME
)

# ── FRS static library ─────────────────────────────────────────────────────
# frs.lib компилирует только наш код (Encoder/Decoder).
# Здесь aff3ct.hpp НЕ включается напрямую, поэтому можно добавить NOGDI
# для защиты от GDI+ на уровне библиотеки.
add_library(frs STATIC
  Encoder_FRS.cpp
  Decoder_FRS1.cpp
  Decoder_FRS2.cpp
)

# nowin/ ПЕРЕД всеми системными путями — заглушает gdiplus.h и компанию
target_include_directories(frs BEFORE PUBLIC
  $<BUILD_INTERFACE:${NOWIN_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

target_link_libraries(frs PUBLIC
  aff3ct::aff3ct-static-lib
  Threads::Threads
)

# ── Executables ────────────────────────────────────────────────────────────
# main_frs1/2 включают <aff3ct.hpp> -> rang.hpp -> wincon.h.
# nowin/ всё равно нужен чтобы заглушить gdiplus.h,
# но NOGDI/NOUSER здесь НЕ добавляем.

add_executable(frs1 main_frs1.cpp)
target_include_directories(frs1 BEFORE PRIVATE
  $<BUILD_INTERFACE:${NOWIN_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(frs1 PRIVATE frs)

add_executable(frs2 main_frs2.cpp)
target_include_directories(frs2 BEFORE PRIVATE
  $<BUILD_INTERFACE:${NOWIN_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_link_libraries(frs2 PRIVATE frs)
