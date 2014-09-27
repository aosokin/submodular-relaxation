%************************************************************************
%                                                                       *
%     LMBM-MEX - Limited Memory Bundle Method for Large-Scale           *
%                Nonsmooth Optimization                                 *
%                                                                       *
%     Original LMBM (fortran sources):                                  *
%       Napsu Karmitsa (maiden name Haarala) 2002 - 2004,               *
%       last modified 2007.                                             *
%                                                                       *
%     MEX-version:                                                      *
%       Seppo Pulkkinen 2009.                                           *
%                                                                       *
%************************************************************************
%     
%
%     * Codes included *
% 
%     MatLab/MEX-sources:
%       lmbm_mex.m    - help file for LMBM-MEX (type "help lmbm_mex" 
%                       in MatLab to read the file).
%       lmbm_driver.m - help file for LMBM_DRIVER (type "help lmbm_diver" 
%                       in MatLab to read the file).
%       lmbm_driver.c - driver for LMBM-MEX.
%       objfunc.c      
%       objfunc.h     
%     
%       testfunc1.m   - an example test problem for LMBM-MEX.
%       testfunc2.m   - an example test problem for LMBM-MEX.
%
%     FORTRAN sources:
%       lmbm_mex.f    - subprogram for MEX -limited memory bundle method.
%       lmbm.f        - limited memory bundle method.
%       lmsub.f       - subprograms for limited memory bundle method.
%       matcal.f      - matrix and vector calculus.
%
%     Other:
%       Makefile      - makefile.
%       examples.txt  - examples for calling LMBM-MEX.
%
%
%     To build the mex file, you need a recent version of MATLAB and also 
%     C and Fortran compilers. A makefile for building the mex file on 
%     Linux is provided in the package. It is known to build with MATLAB 
%     R2007b, gcc 4.2 and gfortran 4.2. Other configurations have not 
%     been yet tested, and we would be pleased to receive feedback from 
%     people with different MATLAB versions and compilers.
%
%     Type "help lmbm_driver" and see examples.txt for more information 
%     of running the program.
%
%
