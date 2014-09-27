function setupSMR
% setupSMR adds all SMR folders to the function search paths
%
% Anton Osokin (firstname.lastname@gmail.com),  22.09.2014

% unzip('http://vision.csd.uwo.ca/code/gco-v3.0.zip', 'mexWrappers/gco-v3.0');

curDir = pwd;
cd(fileparts(mfilename('fullpath')));
addpath(genpath(pwd));
cd(curDir);
