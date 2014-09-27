%build LMBM

% % fortran
% mex lmbm.f -c
% mex lmsub.f -c
% mex lmbmex.f -c
% mex matcal.f -c


% C++: Intel C++
mex lmbm_driver.cpp  objfunc.cpp  lmbm.obj lmbmex.obj lmsub.obj matcal.obj  -largeArrayDims -output lmbm_driver

 