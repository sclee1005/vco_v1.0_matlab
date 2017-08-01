function buildMrfMinimizeMexModifiedForVCO
% mex command to build mrfMinimizeMex
% modified by syshin for VCO 2017.08
trwsPath = 'TRW_S-v1.3/';

srcFiles = { 'mrfMinimizeMexModifiedForVCO.cpp', ...
            fullfile(trwsPath, 'ordering.cpp'), ...
            fullfile(trwsPath, 'MRFEnergy.cpp'), ...
            fullfile(trwsPath, 'treeProbabilities.cpp'), ...
            fullfile(trwsPath, 'minimize.cpp') };
allFiles = '';
for iFile = 1 : length(srcFiles)
    allFiles = [allFiles, ' ', srcFiles{iFile}];
end

cmdLine = ['mex ', allFiles, ' -output mrfMinimizeMexModifiedForVCO -largeArrayDims'];
eval(cmdLine);
