// win_guard.h — принудительно вставляется через /FI перед любым #include
// Цель: включить windows.h с нужными ограничениями ДО того как это сделает
// aff3ct.hpp, чтобы GDI+ не затянулся полностью.
//
// НЕ определяем NOUSER и NOGDI здесь — rang.hpp нуждается в консольных
// функциях (GetConsoleMode, FOREGROUND_RED и т.д.) из wincon.h / winbase.h,
// которые блокируются этими макросами.
// GDI+ блокируется через папку nowin/ с пустыми заглушками.

#pragma once

#if defined(_WIN32) || defined(_WIN64)

#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN   // убирает Multimedia, OLE, Crypto, Shell, WinSock, RPC
#  endif                          // но НЕ убирает wincon.h / winbase.h

#  ifndef NOMINMAX
#    define NOMINMAX              // убирает макросы min/max
#  endif

// NOMCX, NOIME — убираем ненужное, но оставляем консольные API
#  ifndef NOMCX
#    define NOMCX
#  endif
#  ifndef NOIME
#    define NOIME
#  endif

// НЕ определяем NOGDI и NOUSER — они нужны rang.hpp
// GDI+ заблокирован пустыми заглушками в папке nowin/

#  include <windows.h>

#endif // _WIN32
