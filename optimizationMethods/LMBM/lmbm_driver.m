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
%     LMBM_DRIVER 
%
%     * Purpose *
%   
%     Run the LMBM algorithm.
%
%
%     * Calling sequence *
%
%     [X,FVAL,NITER,NFEVAL,TERM,TIME]=LMBM_DRIVER(F,X0,N,PRINT,[MAXTIME,NA,MCU,MC,RPAR,IPAR]) 
%
%     invokes the LMBM algorithm with the user-specified function F, 
%     starting point X0 and problem dimension N. This function returns 
%     the found minimizer X and its corresponding function falue FVAL. 
%     The numbers of used iterations and function evaluations are stored 
%     into NITER, and NFEVAL, respectively. The used computation time is 
%     stored into TIME.
%
%
%     * Required parameters *
%
%     F      A string specifying name of the objective function.
%     X0     starting point for the iteration.
%     N      problem dimension, must be equal to the length of X0.
%     PRINT  verbosity level: 
%              0 = no output.
%              1 = print input parameters and final results.
%
%
%     * Optional parameters *
%
%     MAXTIME  Maximum time until the iteration terminates.
%     NA       Maximum bundle dimension, NA >= 2.
%     MCU      Upper limit for maximum number of stored corrections, MCU >= MC.
%     MC       Maximum number of stored corrections, MC >= 3.
%     RPAR     Array of double scalar parameters.
%     IPAR     Array of integer scalar parameters.
%
%
%     * Detailed parameter description *
%
%     F must be name of the M-file that computes the function value 
%     and (sub)gradient at the given point. The following example 
%     function evaluates f(x)=x^2+y^2 and its gradient at the given 
%     point X:
%
%       function [f, g] = testfunc(x)
%         f = x(1)*x(1) + x(2)*x(2);
%         g(1) = 2*x(1);
%         g(2) = 2*x(2);
%       end
%
%     RPAR:
%       RPAR(1)       Tolerance for change of function values.
%       RPAR(2)       Second Tolerance for change of function values.
%       RPAR(3)       Tolerance for the function value.
%       RPAR(4)       Tolerance for the first termination criterion.
%       RPAR(5)       Tolerance for the second termination criterion.
%       RPAR(6)       Distance measure parameter, 0 <= RPAR(6).
%       RPAR(7)       Line search parameter, 0 < RPAR(7) < 0.25.
%       RPAR(8)       Maximum stepsize, 1 < RPAR(8).
%                       If RPAR(I) <= 0 for I=1,3,4,5,7, and 8 the
%                       default value of the parameter will be used.
%                       If RPAR(2) < 0 the the parameter and the
%                       corresponding termination criterion will be
%                       ignored. If RPAR(2) = 0 default value will
%                       be used. If RPAR(6) < 0 the default value
%                       will be used.
%
%     IPAR:
%       IPAR(1)       Exponent for distance measure.
%       IPAR(2)       Maximum number of iterations.
%       IPAR(3)       Maximum number of function evaluations.
%       IPAR(4)       Maximum number of iterations with changes of
%                       function values smaller than RPAR(1).
%       IPAR(5)       Printout specification:
%                       With mex-version IPAR(5) have to be -1 or MatLab
%                       collapses (code chances this automaticly if a 
%                       different value is given). 
%       IPAR(6)       Selection of the method:
%                       0  - Limited memory bundle method.
%                       1  - L-BFGS bundle method.
%       IPAR(7)       Selection of the scaling:
%                       0  - Scaling at every iteration with STU/UTU.
%                       1  - Scaling at every iteration with STS/STU.
%                       2  - Interval scaling with STU/UTU.
%                       3  - Interval scaling with STS/STU.
%                       4  - Preliminary scaling with STU/UTU.
%                       5  - Preliminary scaling with STS/STU.
%                       6  - No scaling.
%
%     The cause for termination (TERM) is one of the following:
%       1  - The problem has been solved with desired accuracy.
%       2  - Changes in function values < RPAR(1) in IPAR(4) subsequent 
%            iterations.
%       3  - Changes in function values < RPAR(2)*SMALL*MAX(|F_k|,|F_k+1|,1), 
%            where SMALL is the smallest positive number such that 
%            1.0 + SMALL > 1.0.
%       4  - Number of function calls > IPAR(3).
%       5  - Number of iterations > IPAR(2).
%       6  - Time limit exceeded. 
%       7  - F < RPAR(3).
%      -1  - Two consecutive restarts or number of restarts > maximum 
%            number of restarts.
%      -2  - TMAX < TMIN in two subsequent iterations.
%      -4  - Failure in attaining the demanded accuracy.
%      -5  - Invalid input parameters.
%      -6  - Not enough working space.
%
%       
%

