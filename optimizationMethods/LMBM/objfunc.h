
#ifndef OBJFUNC_H

// determine the OS-type - important for the naming conventions
#if defined(_WIN32)
	#define __OS_WIN__
#endif

#ifdef __OS_WIN__
	//Windows specific names
	#define __CURRENT_OBJFUNC_NAME__  OBJFUNC
#else
	//Linux specific names
	#define __CURRENT_OBJFUNC_NAME__  objfunc_
#endif

extern "C" {
    void __CURRENT_OBJFUNC_NAME__(int *n, double *x, double *f, double *g);
}
void set_obj_func(const char *name);

#define OBJFUNC_H
#endif
