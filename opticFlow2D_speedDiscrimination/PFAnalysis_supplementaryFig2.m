function PFAnalysis_supplementaryFig2

% This code aims to replicate Supplementary Figure 2 in the following article:
% 
% Luo, J., Ashida, H., Yokoi, I. & Takemura, H. (2026) Anisotropic speed discrimination sensitivity for expanding and contracting optic flow. Journal of Vision, 26(10):18, 1–17.
%
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)
%             Violinplot-Matlab
%             (https://github.com/bastibe/Violinplot-Matlab)

clear; clc; close all;

%% Optional paths
addpath(genpath('../Violinplot-Matlab-master'));
addpath(genpath('../Palamedes'));

%% Load data and estimate PSE
fileInformation = dir('OpticalFlow*.mat');
if isempty(fileInformation)
    error('No files matching OpticalFlow*.mat were found in the current folder.');
end

subjectNames = cell(numel(fileInformation), 1);
for iFile = 1:numel(fileInformation)
    fileName = fileInformation(iFile).name;
    subjectNames(iFile) = extractBetween(fileName, '_', '_');
end
subjectNames = unique(subjectNames);
subjectNum = numel(subjectNames);

tmp = load(fileInformation(1).name, 'taskData');
speedConds = sort(unique(abs(tmp.taskData.speedConditions(:))))';
nConds = numel(speedConds);

pseExp = nan(subjectNum, nConds);
pseCon = nan(subjectNum, nConds);
fitOKExp = false(subjectNum, nConds);
fitOKCon = false(subjectNum, nConds);

for subjectIndex = 1:subjectNum
    fileList = dir(strcat('*', subjectNames{subjectIndex}, '_*.mat'));
    speedsTest = [];
    speedsComp = [];
    responses = [];

    for iFile = 1:numel(fileList)
        load(fileList(iFile).name, 'taskData');
        speedsTest = [speedsTest, taskData.speedsTestDone]; %#ok<AGROW>
        speedsComp = [speedsComp, round(taskData.speedsCompDone, 4)]; %#ok<AGROW>
        responses = [responses, taskData.responses]; %#ok<AGROW>
    end

    speedsCompAbs = abs(speedsComp);

    for cond = 1:nConds
        thisSpeed = speedConds(cond);

        idxExp = speedsTest == thisSpeed;
        if any(idxExp)
            [xExp, nPosExp, nTotalExp] = PAL_PFML_GroupTrialsbyX( ...
                speedsCompAbs(idxExp), responses(idxExp), ones(1, sum(idxExp)));
            [pseExp(subjectIndex, cond), fitOKExp(subjectIndex, cond)] = ...
                fitCumulativeNormalPSE(log2(xExp), nPosExp, nTotalExp);
        end

        idxCon = speedsTest == -thisSpeed;
        if any(idxCon)
            [xCon, nPosCon, nTotalCon] = PAL_PFML_GroupTrialsbyX( ...
                speedsCompAbs(idxCon), responses(idxCon), ones(1, sum(idxCon)));
            [pseCon(subjectIndex, cond), fitOKCon(subjectIndex, cond)] = ...
                fitCumulativeNormalPSE(log2(xCon), nPosCon, nTotalCon);
        end
    end
end

%% Central tendency: PSE bias in octaves
% Positive values indicate PSE > reference speed; negative values indicate
% PSE < reference speed.
pseBiasExp = log2(pseExp ./ repmat(speedConds, subjectNum, 1));
pseBiasCon = log2(pseCon ./ repmat(speedConds, subjectNum, 1));

medianBiasExp = median(pseBiasExp, 1, 'omitnan');
medianBiasCon = median(pseBiasCon, 1, 'omitnan');

CIBias_exp = nan(2, nConds);
CIBias_con = nan(2, nConds);
rng(1);
for cond = 1:nConds
    CIBias_exp(:, cond) = bootstrapMedianCI(pseBiasExp(:, cond));
    CIBias_con(:, cond) = bootstrapMedianCI(pseBiasCon(:, cond));
end

%% Figure: central tendency / PSE bias
figure('Color', [1, 1, 1], 'Position', [300, 200, 850, 330], ...
    'Name', sprintf('Central tendency, N = %d', subjectNum));
hold on

for cond = 1:nConds
    xExp = pseBiasExp(:, cond);
    xCon = pseBiasCon(:, cond);

    Violin({xExp(~isnan(xExp))}, cond - 0.05, ...
        'ShowData', true, 'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, ...
        'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[1, 0, 0]}, ...
        'Width', 0.22, 'ViolinAlpha', {0.5}, 'HalfViolin', 'left');

    Violin({xCon(~isnan(xCon))}, cond + 0.05, ...
        'ShowData', true, 'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, ...
        'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[0, 0, 1]}, ...
        'Width', 0.22, 'ViolinAlpha', {0.5}, 'HalfViolin', 'right');

    errNegExp = medianBiasExp(cond) - CIBias_exp(1, cond);
    errPosExp = CIBias_exp(2, cond) - medianBiasExp(cond);
    errNegCon = medianBiasCon(cond) - CIBias_con(1, cond);
    errPosCon = CIBias_con(2, cond) - medianBiasCon(cond);

    errorbar(cond - 0.3, medianBiasExp(cond), errNegExp, errPosExp, ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], ...
        'MarkerEdgeColor', [1, 0, 0], 'Color', 'r', 'LineWidth', 2);
    errorbar(cond + 0.3, medianBiasCon(cond), errNegCon, errPosCon, ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], ...
        'MarkerEdgeColor', [0, 0, 1], 'Color', 'b', 'LineWidth', 2);
end

yL = yline(0, ':k');
uistack(yL, 'bottom');
set(gca, 'TickDir', 'out', 'Box', 'off', ...
    'XTick', 1:nConds, 'XTickLabel', compose('%g', speedConds), 'FontSize', 12);
title('Central tendency / PSE bias');
xlabel('Speed condition (deg/s)');
ylabel('PSE bias, log_2(PSE/reference speed)');

pExp = plot(nan, nan, 's', 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [1, 0, 0], 'MarkerSize', 10, 'DisplayName', 'Expansion');
pCon = plot(nan, nan, 's', 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0, 0, 1], 'MarkerSize', 10, 'DisplayName', 'Contraction');
legend([pExp pCon], 'Location', 'best', 'Box', 'off', 'FontSize', 10);

end

%% Local helper functions
function [pseLinear, fitOK] = fitCumulativeNormalPSE(xLog2, numPos, numTotal)
% Fit cumulative-normal psychometric function on log2 speed.
% Return PSE in linear speed units.

    PF = @PAL_CumulativeNormal;
    paramsFree = [1 1 0 0];

    xLog2 = xLog2(:)';
    numPos = numPos(:)';
    numTotal = numTotal(:)';

    searchGrid.alpha = min(xLog2):0.01:max(xLog2);
    searchGrid.beta = logspace(-1, 1, 100);
    searchGrid.gamma = 0.00;
    searchGrid.lambda = 0.00;

    [paramsValues, ~, exitFlag] = PAL_PFML_Fit( ...
        xLog2, numPos, numTotal, searchGrid, paramsFree, PF);

    fitOK = exitFlag > 0;
    pseLinear = 2.^PF(paramsValues(1, :), 0.50, 'Inverse');
end

function ci = bootstrapMedianCI(x)
% 95% BCa bootstrap CI for the median, with a percentile fallback.

    x = x(:);
    x = x(~isnan(x));

    if numel(x) < 2
        ci = [NaN; NaN];
        return;
    end

    try
        ci = bootci(10000, {@median, x}, 'type', 'bca');
    catch
        bootStat = bootstrp(10000, @median, x);
        ci = prctile(bootStat, [2.5 97.5])';
    end
end
