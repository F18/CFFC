#if defined(__MACH__)
# include <sys/malloc.h>
#else
# include <malloc.h>
#endif
#define new(type,size)	(type *) malloc((unsigned) (size)*sizeof(type))
#define renew(type,size,old) (type *)realloc(old,(unsigned)(size)*sizeof(type))
