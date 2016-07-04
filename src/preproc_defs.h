#ifndef PREPROC_DEFS_H
#define PREPROC_DEFS_H

#if defined(_MSC_VER)
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

#if defined (__GNUC__) && !defined (__clang__)
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 500
#endif
#ifndef _GNU_SOURCES
#define _GNU_SOURCES
#endif
#endif


#endif // PREPROC_DEFS_H


