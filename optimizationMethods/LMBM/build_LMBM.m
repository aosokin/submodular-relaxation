function build_LMBM

% fortran
mex lmbm.f -c 
mex lmsub.f -c
mex lmbmex.f -c
mex matcal.f -c

isWindows = strcmp(mexext, 'mexw32') || strcmp(mexext, 'mexw64');

objExt = 'o';
if isWindows
    objExt = 'obj';
end
    
objFiles = {'lmbm', 'lmbmex', 'lmsub', 'matcal'};

mexCmd = 'mex lmbm_driver.cpp  objfunc.cpp -largeArrayDims -output lmbm_driver';

if ~isWindows
mexCmd = [mexCmd, ' -lgfortran'];
end

for iFile = 1 : length(objFiles)
    mexCmd = [mexCmd, ' ', objFiles{iFile}, '.', objExt, ' '];
end

eval(mexCmd);

deleteCmd = 'delete ';
for iFile = 1 : length(objFiles)
    deleteCmd = [deleteCmd, ' ', objFiles{iFile}, '.', objExt, ' '];
end
eval(deleteCmd);
 
