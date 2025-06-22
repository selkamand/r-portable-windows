/* Rversion.h.  Generated automatically. */
#ifndef R_VERSION_H
#define R_VERSION_H

#ifdef __cplusplus
extern "C" {
#endif

#define R_VERSION 263425
#define R_NICK "Great Square Root"
#define R_Version(v,p,s) (((v) * 65536) + ((p) * 256) + (s))
#define R_MAJOR  "4"
#define R_MINOR  "5.1"
#define R_STATUS ""
#define R_YEAR   "2025"
#define R_MONTH  "06"
#define R_DAY    "13"
#define R_SVN_REVISION 88306
#ifdef __llvm__
# define R_FILEVERSION    4,51,22770,0
#else
# define R_FILEVERSION    4,51,88306,0
#endif

#ifdef __cplusplus
}
#endif

#endif /* not R_VERSION_H */
