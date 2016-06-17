#ifndef PREPROC_DEFS_H
#define PREPROC_DEFS_H

#if defined(_MSC_VER)
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#endif
#if defined (__GNUC__) && !defined (__clang__)
#define _XOPEN_SOURCE 500
#define _GNU_SOURCES
#endif


#endif // PREPROC_DEFS_H
