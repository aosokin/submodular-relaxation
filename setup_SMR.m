function setup_SMR
% setup_SMR adds all SMR folders to the function search paths and compiles all the required binaries
%
% Anton Osokin (firstname.lastname@gmail.com),  22.09.2014


curDir = pwd;
smrRootDir = fileparts(mfilename('fullpath'));
cd(smrRootDir);
addpath(smrRootDir);
addpath(genpath(fullfile(smrRootDir,'data')));
addpath(genpath(fullfile(smrRootDir,'mexWrappers')));
addpath(genpath(fullfile(smrRootDir,'optimizationMethods')));
addpath(genpath(fullfile(smrRootDir,'oracles')));
addpath(genpath(fullfile(smrRootDir,'utils')));
addpath(genpath(fullfile(smrRootDir,'experiments')));
cd(curDir);

build_SMR(false);

end
