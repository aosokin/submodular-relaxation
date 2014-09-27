*************************************************************************
*                                                                       *
*     LMBM-MEX - Limited Memory Bundle Method for Large-Scale           *
*                Nonsmooth Optimization                                 *
*                                                                       *
*************************************************************************
*     
*
*     Codes included:
* 
*     MatLab/MEX-sources
*       lmbm_driver.m - help file for LMBM-MEX (type "help lmbm_diver" 
*                       in MatLab to read the file).
*       lmbm_driver.c - driver for LMBM-MEX. 
*       objfunc.c      
*       objfunc.h     
*     
*       testfunc1.m   - an example test problem for LMBM-MEX.
*       testfunc2.m   - an example test problem for LMBM-MEX.
*
*     FORTRAN sources:
*       lmbm_mex.f    - subprogram for MEX -limited memory bundle method.
*       lmbm.f        - limited memory bundle method.
*       lmsub.f       - subprograms for limited memory bundle method.
*       matcal.f      - matrix and vector calculus.
*
*     Other
*       Makefile      - makefile.
*       examples.txt  - examples for calling LMBM-MEX
*
*
*     To build the mex file, you need a recent version of MATLAB and also 
*     C and Fortran compilers. A makefile for building the mex file on 
*     Linux is provided in the package. It is known to build with MATLAB 
*     R2007b, gcc 4.2 and gfortran 4.2. Other configurations have not 
*     been yet tested, and we would be pleased to receive feedback from 
*     people with different MATLAB versions and compilers.
*
*
*     References:
*
*     M.Haarala, K.Miettinen and M.M.Mäkelä: "New Limited Memory Bundle 
*     Method for Large-Scale Nonsmooth Optimization". Optimization 
*     Methods and Software 19(6): 673 - 692, 2004.
*
*     N.Haarala, K.Miettinen and M.M.Mäkelä: "Globally Convergent 
*     Limited Memory Bundle Method for Large-Scale Nonsmooth 
*     Optimization". Mathematical Programming, Vol. 109, No. 1, 
*     pp. 181-205, 2007.
*
*     J.Vlcek, L.Luksan: Globally Convergent Variable Metric
*     Method for Nonconvex Nondifferentiable Unconstrained
*     Minimization. Journal of Optimization Theory and Applications
*     111(2): 407 - 430, 2001.
*
*
*
*     Original LMBM (fortran sources):
*           Napsu Karmitsa (maiden name Haarala) 2002 - 2004
*           last modified 2007      
*
*     MEX-version:
*           Seppo Pulkkinen 2009
*
*
*************************************************************************
*
*     Remark:
*
*     At the beginning of each FORTRAN- file, there is a list of the 
*     subroutines and functions included to that file. Moreover, at the 
*     beginning of each subroutine, there is a description of parameters 
*     used in the routine (at least those needed in calling sequences). 
*     The types of the parameters (or arguments) are introduced with two
*     letters. The first letter is either I for integer arguments or R 
*     for double precision real arguments.  The second letter specifies 
*     whether the argument must have a value defined on the entry to
*     the subroutine (I), whether it is a value which will be returned
*     (O), or both (U), or whether it is an auxiliary value (A). Note 
*     that the arguments of the types II and RI can be changed on output
*     under some circumstances: especially, if improper input values are
*     given or if set zero. In the latter case the default values will 
*     be used (if applicable).
*
*
*************************************************************************
*************************************************************************
*************************************************************************
*
*
*     lmbm_mex.f includes the following subroutines
*
*     S   FUNDER        Computation of the value and the subgradient 
*                         of the objective function.
*      
*
************************************************************************
*      
*     * SUBROUTINE FUNDER *
*
*      
*     * Purpose *
*
*     Computation of the value and the subgradient of the objective
*     function.
*
*      
*     * Calling sequence *
*
*     CALL FUNDER(N,X,F,G,ITERM)
*
*
*     * Parameters *
*      
*     II  N             Number of variables.
*     RI  X(N)          Vector of variables.
*     RO  F             Value of the objective function.
*     RO  G(N)          Subgradient of the objective function.
*     IO  ITERM         Cause of termination:
*                          0  - Everything is ok.
*                         -3  - Failure in function or subgradient
*                               calculations (assigned by the user).
*
*
*     * Variables in COMMON /PROB/ *
*
*     I   NEXT          Number of the test problem.
*
*
*     * Subprograms used *
*      
*     S   FUNC          Computation of the value and the subgradient for
*                       problem next.
*     
*      
*     Napsu Haarala (2002-2004)
*     Mex-version Seppo Pulkkinen 2009
*      
     
      SUBROUTINE FUNDER(N,X,F,G,ITERM)

*     Scalar Arguments
      INTEGER N,ITERM
      DOUBLE PRECISION F

*     Array Arguments
      DOUBLE PRECISION G(*),X(*)

     
*     
*     Function and subgradient evaluation
*    

      CALL objfunc(N, X, F, G)
      
      RETURN
      END

