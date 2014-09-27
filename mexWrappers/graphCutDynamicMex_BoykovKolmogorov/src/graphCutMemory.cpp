#include "graphCutMemory.h"

/* memory management */
void* operator new(size_t size)
{
    void *ptr = NULL;
//    mexWarnMsgTxt("Overloaded new operator");
    ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}
void* operator new[](size_t size)
{
    void *ptr = NULL;
//    mexWarnMsgTxt("Overloaded new[] operator");
    ptr = mxMalloc(size);
    mexMakeMemoryPersistent(ptr);
    return ptr;
}
void operator delete(void* ptr)
{
//    mexWarnMsgTxt("Overloaded delete operator");
    mxFree(ptr);
}
void operator delete[](void* ptr)
{
//    mexWarnMsgTxt("Overloaded delete[] operator");
    mxFree(ptr);
}

GraphType* getGraphHandle(const mxArray *x)
{
    GraphHandle gch = 0;
    GraphType* g = 0;
    
    if ( mxGetClassID(x) != MATLAB_POINTER_TYPE ) {
        mexErrMsgIdAndTxt("graphCutMemory:handleWrongType", "Graph handle argument is not of proper type");
    }
	if ( mxGetNumberOfElements(x) != 1 ) {
        mexErrMsgIdAndTxt("graphCutMemory:handleWrongSize", "Too many graph handles");
    }
    
    gch = (GraphHandle*)mxGetData(x);
	g = (GraphType*)(*(POINTER_CAST*)gch);
    if ( g == NULL ) {
        mexErrMsgIdAndTxt("graphCutMemory:badHandle", "Graph handle is not valid");
    }
    return g;
}
