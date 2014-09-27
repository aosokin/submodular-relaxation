
#ifndef _GRAPH_CUT_MEMORY_H_
#define _GRAPH_CUT_MEMORY_H_

#include <tmwtypes.h>
#include <limits>
#include <cmath>

#include "graphCutMex.h"
#include "mex.h"

#define INFTY INT_MAX

int round(double a);
int isInteger(double a);

//define types
typedef double EnergyType;
//mxClassID MATLAB_ENERGYTERM_TYPE = mxDOUBLE_CLASS;
#define MATLAB_ENERGYTERM_TYPE (mxDOUBLE_CLASS)

typedef double EnergyTermType;
//mxClassID MATLAB_ENERGY_TYPE = mxDOUBLE_CLASS;
#define MATLAB_ENERGY_TYPE  (mxDOUBLE_CLASS)

typedef double LabelType;
//mxClassID MATLAB_LABEL_TYPE = mxDOUBLE_CLASS;
#define MATLAB_LABEL_TYPE  (mxDOUBLE_CLASS)

typedef Graph<EnergyTermType,EnergyTermType,EnergyType> GraphType; 

typedef void* GraphHandle;

/* pointer types in 64 bits machines */
#ifdef A64BITS
#define MATLAB_POINTER_TYPE mxUINT64_CLASS
#else
#define MATLAB_POINTER_TYPE mxUINT32_CLASS
#endif

#ifdef A64BITS
#define POINTER_CAST    int64_T
#else
#define POINTER_CAST    int
#endif

template<class T>
void GetScalar(const mxArray* x, T& scalar)
{
    if ( mxGetNumberOfElements(x) != 1 )
        mexErrMsgIdAndTxt("graphCutMemory:GetScalar","input is not a scalar!");
    void *p = mxGetData(x);
    switch (mxGetClassID(x)) {
        case mxCHAR_CLASS:
            scalar = *(char*)p;
            break;
        case mxDOUBLE_CLASS:
            scalar = *(double*)p;
            break;
        case mxSINGLE_CLASS:
            scalar = *(float*)p;
            break;
        case mxINT8_CLASS:
            scalar = *(char*)p;
            break;
        case mxUINT8_CLASS:
            scalar = *(unsigned char*)p;
            break;
        case mxINT16_CLASS:
            scalar = *(short*)p;
            break;
        case mxUINT16_CLASS:
            scalar = *(unsigned short*)p;
            break;
        case mxINT32_CLASS:
            scalar = *(int*)p;
            break;
        case mxUINT32_CLASS:
            scalar = *(unsigned int*)p;
            break;
#ifdef A64BITS            
        case mxINT64_CLASS:
            scalar = *(int64_T*)p;
            break;
        case mxUINT64_CLASS:
            scalar = *(uint64_T*)p;
            break;
#endif /* 64 bits machines */            
        default:
            mexErrMsgIdAndTxt("graphCutMemory:GetScalar","unsupported data type");
    }
}

/* memory allocations - redirect to MATLAB memory menager */
void* operator new(size_t size);
void* operator new[](size_t size);
void operator delete(void* ptr);
void operator delete[](void* ptr);

inline int round(double a)
{
	return (int)floor(a + 0.5);
}

inline int isInteger(double a)
{
	return (a - round(a) < 1e-6);
}

GraphType* getGraphHandle(const mxArray *x); // extract handle from mxArray 

#endif /* _GRAPH_CUT_MEMORY_H_ */
