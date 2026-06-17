cd /D E:\1-ННГУ\6-ВКР\2-Практика\1_RS_old

cmake -B build -G "Visual Studio 17 2022" -A x64 ^
  -DCMAKE_PREFIX_PATH="E:/1-ННГУ/6-ВКР/2-Практика/aff3ct-3.0.2/lib/cmake/aff3ct-3.0.2" ^
  -DCMAKE_CXX_FLAGS="/EHsc /MP4 /D_CRT_SECURE_NO_DEPRECATE"

cmake --build build --config Debug --target frs1

frs1.exe 63 255 3 0 llr