//
// Created by jhewers on 11/07/22.
//

#ifndef PYVIEWSHED_LOGGING_H
#define PYVIEWSHED_LOGGING_H

#include <loguru.hpp>

// https://emilk.github.io/loguru/index.html

#define _LOG_F LOG_F
#define _ABORT_F ABORT_F
#define _VLOG_F VLOG_F
#define _DLOG_F DLOG_F

#define _LOG_SCOPE_FUNCTION LOG_SCOPE_FUNCTION
#define _LOG_SCOPE_F LOG_SCOPE_F
#define _VLOG_SCOPE_F VLOG_SCOPE_F

#define _LOG_WARNING(...) LOG_F(WARNING, __VA_ARGS__)
#define _LOG_INFO(...) LOG_F(INFO, __VA_ARGS__)
#define _LOG_VERBOSE(...) LOG_F(1, __VA_ARGS__)

#endif  // PYVIEWSHED_LOGGING_H
