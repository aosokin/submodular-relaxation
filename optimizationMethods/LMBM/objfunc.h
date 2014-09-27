
#ifndef OBJFUNC_H

//void objfunc_(int *n, double *x, double *f, double *g);
extern "C" {
    void OBJFUNC(int *n, double *x, double *f, double *g);
}
void set_obj_func(const char *name);

#define OBJFUNC_H
#endif
