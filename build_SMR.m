function build_SMR(forceBuild)
% build_SMR builds all mex-files required for SMR
%
% Input:
%   forceBuild - if true - rebuilds all the binaries of the package;
%                  if false - rebuils only missing binaries.
%                  (default: false)
%
% Anton Osokin (firstname.lastname@gmail.com),  22.09.2014

if ~exist('forceBuild', 'var')
    forceBuild = false;
end
if ~islogical(forceBuild) || ~isscalar(forceBuild)
    error('build_SMR:badForceBuild', 'forceBuild input should be a single logical value');
end

curDir = pwd;
smrRootDir = fileparts(mfilename('fullpath'));

if exist('graphCutMex', 'file') ~= 3  ||  forceBuild
    % build graphCutMex_BoykovKolmogorov
    fprintf('Building graphCutMex...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'graphCutMex_BoykovKolmogorov'));
    build_graphCutMex;
    cd(curDir);
end

if exist('graphCutDynamicMex', 'file') ~= 3 || ...
   exist('updateUnaryGraphCutDynamicMex', 'file') ~= 3 || ...
   exist('deleteGraphCutDynamicMex', 'file') ~= 3  ||  forceBuild
    % build graphCutDynamicMex_BoykovKolmogorov
    fprintf('Building graphCutDynamicMex...\n')
    cd(fullfile(smrRootDir, 'mexWrappers','graphCutDynamicMex_BoykovKolmogorov'));
    build_graphCutDynamicMex;
    cd(curDir);
end

if exist('icmPottsMex', 'file') ~= 3  ||  forceBuild
    % build icmPottsMex
    fprintf('Building icmPottsMex...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'icmPottsMex'));
    build_icmPottsMex;
    cd(curDir);
end

if exist('qpboMex', 'file') ~= 3  ||  forceBuild
    % build qpboMex
    fprintf('Building qpboMex...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'qpboMex'));
    build_qpboMex;
    cd(curDir);
end

if exist('trwsMex_time', 'file') ~= 3  ||  forceBuild
    % build trwsMex_time
    fprintf('Building trwsMex_time...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'trwsMex_time'));
    build_trwsMex_time;
    cd(curDir);
end

if exist('viterbiPottsMex', 'file') ~= 3  ||  forceBuild
    % build viterbiPottsMex
    fprintf('Building viterbiPottsMex_time...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'viterbiPottsMex'));
    build_viterbiPottsMex;
    cd(curDir);
end

if exist('gco_matlab', 'file') ~= 3  ||  forceBuild
	% build gco_matlab
    if exist('GCO_BuildLib.m', 'file') ~= 2
        fprintf('Downloading gco_matlab...\n')
        unzip('http://vision.csd.uwo.ca/code/gco-v3.0.zip', fullfile(smrRootDir, 'mexWrappers', 'gco-v3.0'));
        addpath(fullfile(smrRootDir, 'mexWrappers', 'gco-v3.0'));
        addpath(fullfile(smrRootDir, 'mexWrappers', 'gco-v3.0', 'matlab'));
    end
    fprintf('Building gco_matlab...\n')
    cd(fullfile(smrRootDir, 'mexWrappers', 'gco-v3.0', 'matlab'));
    GCO_BuildLib;
	cd(curDir);
end

if exist('lmbm_driver', 'file') ~= 3 ||  forceBuild
    fprintf('Building LMBM. Fortran compiler with compatible C++ compaler are required. A valid combination is Intel Visual Fortran Composer XE 2013 and Intel C++ Composer XE 2013.')
    cd(fullfile(smrRootDir, 'optimizationMethods', 'LMBM'));
	build_LMBM;
    cd(curDir);
end

end
