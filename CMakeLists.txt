cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(
  aff3ct-ldpc-list
  LANGUAGES CXX
  DESCRIPTION "LDPC List-BP decoder with AFF3CT"
  VERSION 0.0.1)

# --------------------------------------------------------------------------
# Threads
# --------------------------------------------------------------------------
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# --------------------------------------------------------------------------
# AFF3CT
# If not found automatically, pass the path explicitly:
#   cmake .. -DAFF3CT_DIR="C:/path/to/aff3ct/lib/cmake/aff3ct-3.0.2"
# --------------------------------------------------------------------------
find_package(AFF3CT 3.0.2 REQUIRED)

# --------------------------------------------------------------------------
# Executable
#
# Files in this project:
#   LDPC_list_main.cpp          – main / BER simulation
#   Decoder_LDPC_BP_list.hpp    – List-BP decoder declaration
#   Decoder_LDPC_BP_list.hxx    – List-BP decoder implementation (header-only)
#
# Decoder_LDPC_BP_list.hpp/.hxx are header-only (template class).
# They do NOT need to be listed as source files for compilation,
# but listing the .hpp makes IDEs (VS, CLion) show it in the project tree.
# --------------------------------------------------------------------------
add_executable(aff3ct-ldpc-list-run
    LDPC_list_main.cpp
    Decoder_LDPC_BP_list.hpp   # header-only – included for IDE visibility
)

target_link_libraries(aff3ct-ldpc-list-run
    PUBLIC
        aff3ct::aff3ct-static-lib
        Threads::Threads
)

# Source directory contains our decoder headers and the AFF3CT umbrella header
target_include_directories(aff3ct-ldpc-list-run
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# --------------------------------------------------------------------------
# Optional: nicer compiler warnings
# --------------------------------------------------------------------------
if(MSVC)
    target_compile_options(aff3ct-ldpc-list-run PRIVATE /W3)
else()
    target_compile_options(aff3ct-ldpc-list-run PRIVATE -Wall -Wextra -Wpedantic)
endif()
