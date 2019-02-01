#ifndef VERBOSE_H_
#define VERBOSE_H_


//#define VERBOSE

//#define VERBOSE_DYN

#ifdef VERBOSE
	#include<cstdio>
	#ifdef VERBOSE_DYN
		extern bool verbose_mode;
		//#define eprintf(format, ...) if(verbose_mode) printf(format, ##__VA_ARGS__);
		#define eprintf(...) if(verbose_mode) printf(__VA_ARGS__);
	#else
		//#define eprintf(format, ...) printf(format, ##__VA_ARGS__);
		#define eprintf(...) printf(__VA_ARGS__);
	#endif
#else
	//#define eprintf(format, ...);
	#define eprintf(...);
#endif





#endif
