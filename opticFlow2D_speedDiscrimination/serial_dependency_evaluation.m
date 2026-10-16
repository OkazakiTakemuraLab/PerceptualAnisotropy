clear; close all;

%% Set path
% addpath(genpath('../Violinplot-Matlab-master'));
% addpath(genpath('../boxplotGroup'));

%% Load data and calculate CP(0.5) and CP(0.75)

fileInformation = dir('OpticalFlow*.mat');
fileNum = length(fileInformation);
subjectNames = cell(fileNum, 1);
for i = 1:length(fileInformation)
    fileName = fileInformation(i).name;
    subjectNames(i) = extractBetween(fileName, '_', '_');
end
subjectNames = unique(subjectNames);
subjectNum = length(subjectNames);

tic;
for subjectIndex = 1:subjectNum
    
    fprintf('\n***************\nProcessing %s\n***************\n', subjectNames{subjectIndex})
    
    fileList = dir(strcat('*', subjectNames{subjectIndex}, '_*.mat'));

    speedsRefe = [];
    speedsTest = [];
    responses  = [];
    for i = 1:numel(fileList)
        load(fileList(i).name);
        speedsRefe = [speedsRefe, taskData.speedsTestDone]; %#ok<AGROW>
        speedsTest = [speedsTest, round(taskData.speedsCompDone, 4)]; %#ok<AGROW>
        responses  = [responses, taskData.responses]; %#ok<AGROW>
    end

    speedConditions = taskData.speedConditions;
    speedLevels = unique(abs(speedConditions), 'stable');
    speedCondNum = numel(speedLevels);
    
    if subjectIndex == 1
        thresholdsConLow = nan(subjectNum, speedCondNum);
        disRangeConLow = nan(subjectNum, speedCondNum);
        thresholdsExpLow = nan(subjectNum, speedCondNum);
        disRangeExpLow = nan(subjectNum, speedCondNum);

        thresholdsConSam = nan(subjectNum, speedCondNum);
        disRangeConSam = nan(subjectNum, speedCondNum);
        thresholdsExpSam = nan(subjectNum, speedCondNum);
        disRangeExpSam = nan(subjectNum, speedCondNum);

        thresholdsConHig = nan(subjectNum, speedCondNum);
        disRangeConHig = nan(subjectNum, speedCondNum);
        thresholdsExpHig = nan(subjectNum, speedCondNum);
        disRangeExpHig = nan(subjectNum, speedCondNum);

        thresholdsConAll = nan(subjectNum, speedCondNum);
        disRangeConAll = nan(subjectNum, speedCondNum);
        thresholdsExpAll = nan(subjectNum, speedCondNum);
        disRangeExpAll = nan(subjectNum, speedCondNum);
    end

    prevRefe = [NaN, abs(speedsRefe(1:end-1))];
    currRefe = abs(speedsRefe);
    spdDiff = currRefe - prevRefe;

    [thresholdsCon, disRangeCon, thresholdsExp, disRangeExp] = ...
        psychometricCal('Low', spdDiff, speedsRefe, speedsTest, responses, speedLevels);
    thresholdsConLow(subjectIndex, :) = thresholdsCon;
    disRangeConLow(subjectIndex, :)   = disRangeCon;
    thresholdsExpLow(subjectIndex, :) = thresholdsExp;
    disRangeExpLow(subjectIndex, :)   = disRangeExp;

    [thresholdsCon, disRangeCon, thresholdsExp, disRangeExp] = ...
        psychometricCal('Sam', spdDiff, speedsRefe, speedsTest, responses, speedLevels);
    thresholdsConSam(subjectIndex, :) = thresholdsCon;
    disRangeConSam(subjectIndex, :)   = disRangeCon;
    thresholdsExpSam(subjectIndex, :) = thresholdsExp;
    disRangeExpSam(subjectIndex, :)   = disRangeExp;

    [thresholdsCon, disRangeCon, thresholdsExp, disRangeExp] = ...
        psychometricCal('Hig', spdDiff, speedsRefe, speedsTest, responses, speedLevels);
    thresholdsConHig(subjectIndex, :) = thresholdsCon;
    disRangeConHig(subjectIndex, :)   = disRangeCon;
    thresholdsExpHig(subjectIndex, :) = thresholdsExp;
    disRangeExpHig(subjectIndex, :)   = disRangeExp;

    [thresholdsCon, disRangeCon, thresholdsExp, disRangeExp] = ...
        psychometricCal('All', spdDiff, speedsRefe, speedsTest, responses, speedLevels);
    thresholdsConAll(subjectIndex, :) = thresholdsCon;
    disRangeConAll(subjectIndex, :)   = disRangeCon;
    thresholdsExpAll(subjectIndex, :) = thresholdsExp;
    disRangeExpAll(subjectIndex, :)   = disRangeExp;

end

elapsedTime = toc;
disp(['Process took ', num2str(elapsedTime/60), ' min.'])

save('serial_dependency_data.mat');

%% Helper function

function [thresholdsCon, disRangeCon, thresholdsExp, disRangeExp] = ...
    psychometricCal(serialCond, spdDiff, speedsRefe, speedsTest, responses, speedLevels)

serialCond = char(serialCond);

filterFlag = 1;

switch serialCond
    case 'Low' % previous reference faster than current
        spdSelect = spdDiff < 0;
    case 'Sam' % previous reference the same as current
        spdSelect = spdDiff == 0;
    case 'Hig' % previous reference slower than current
        spdSelect = spdDiff > 0;
    case 'All'
        spdSelect = true(1, numel(spdDiff));
        filterFlag = 0;
    otherwise
        error('serialCond must be ''Low'', ''Sam'', ''Hig'', or ''All''.');
end

speedsRefe = speedsRefe(spdSelect);
speedsTest = speedsTest(spdSelect);
responses  = responses(spdSelect);

speedCondNum = numel(speedLevels);

thresholdsCon = nan(1, speedCondNum);
disRangeCon   = nan(1, speedCondNum);
thresholdsExp = nan(1, speedCondNum);
disRangeExp   = nan(1, speedCondNum);

PF = @PAL_CumulativeNormal;
paramsFree = [1 1 0 0];
searchGrid.beta = logspace(-1, 1, 100);
searchGrid.gamma = 0.00;
searchGrid.lambda = 0.00;
B = 500;

for n = 1:speedCondNum
    thisSpeed = speedLevels(n);

    % Contraction trials
    conMask = speedsRefe == -thisSpeed;
    cmpCon = speedsTest(conMask);
    respCon = responses(conMask);

    if ~isempty(cmpCon)
        [stiLevelsCon, numPosCon, numTotalCon] = ...
            PAL_PFML_GroupTrialsbyX(cmpCon, respCon, ones(1, numel(cmpCon)));

        stiLevelsCon = -stiLevelsCon; % make positive before log2

        if ~isempty(stiLevelsCon) && ...
           all(stiLevelsCon > 0) && ...
           numel(unique(stiLevelsCon)) >= 2 && ...
           any(numPosCon > 0) && any(numPosCon < numTotalCon)

            stiLevelsCon = log2(stiLevelsCon);
            searchGrid.alpha = min(stiLevelsCon):0.01:max(stiLevelsCon);
            if isempty(searchGrid.alpha)
                searchGrid.alpha = min(stiLevelsCon);
            end

            [paramsValues, ~, exitFlag] = ...
                PAL_PFML_Fit(stiLevelsCon, numPosCon, numTotalCon, ...
                searchGrid, paramsFree, PF);
            
            if filterFlag
                if exitFlag > 0    
                    thresholdsCon(n) = 2^PF(paramsValues, 0.5, 'Inverse');
                    disRangeCon(n)   = 2^PF(paramsValues, 0.75, 'Inverse');
                end
            else
                thresholdsCon(n) = 2^PF(paramsValues, 0.5, 'Inverse');
                disRangeCon(n)   = 2^PF(paramsValues, 0.75, 'Inverse');
            end
        end
    end

    % Expansion trials
    expMask = speedsRefe == thisSpeed;
    cmpExp = speedsTest(expMask);
    respExp = responses(expMask);

    if ~isempty(cmpExp)
        [stiLevelsExp, numPosExp, numTotalExp] = ...
            PAL_PFML_GroupTrialsbyX(cmpExp, respExp, ones(1, numel(cmpExp)));

        if ~isempty(stiLevelsExp) && ...
           all(stiLevelsExp > 0) && ...
           numel(unique(stiLevelsExp)) >= 2 && ...
           any(numPosExp > 0) && any(numPosExp < numTotalExp)

            stiLevelsExp = log2(stiLevelsExp);
            searchGrid.alpha = min(stiLevelsExp):0.01:max(stiLevelsExp);
            if isempty(searchGrid.alpha)
                searchGrid.alpha = min(stiLevelsExp);
            end

            [paramsValues, ~, exitFlag] = ...
                PAL_PFML_Fit(stiLevelsExp, numPosExp, numTotalExp, ...
                searchGrid, paramsFree, PF);

            if filterFlag

                if exitFlag > 0
                    thresholdsExp(n) = 2^PF(paramsValues, 0.5, 'Inverse');
                    disRangeExp(n)   = 2^PF(paramsValues, 0.75, 'Inverse');
                end
            else
                thresholdsExp(n) = 2^PF(paramsValues, 0.5, 'Inverse');
                disRangeExp(n)   = 2^PF(paramsValues, 0.75, 'Inverse');
            end
        end
    end
    
end

end
