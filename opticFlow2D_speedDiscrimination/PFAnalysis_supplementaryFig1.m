function PFAnalysis_supplementaryFig1

% This code aims to replicate Supplementary Figure 1 in the following article:
% 
% Luo, J., Ashida, H., Yokoi, I. & Takemura, H. (2026) Anisotropic speed discrimination sensitivity for expanding and contracting optic flow. Journal of Vision, 26(10):18, 1–17.
%
% Note: The N reported in the final paired statistics may be smaller than the
% total number of participants because trial splitting by previous-speed
% condition can leave too few usable trials for some subject x speed x direction
% psychometric fits, and estimates were retained only when the fit was valid.
% Paired tests therefore include only subjects with valid values for both
% expansion and contraction.
%
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)
%             Violinplot-Matlab
%             (https://github.com/bastibe/Violinplot-Matlab)

clear; close all; clc;

%% Set path
addpath(genpath('../Violinplot-Matlab-master'));
addpath(genpath('../boxplotGroup'));
addpath(genpath('../Palamedes'));

%% Load data

load('serial_dependency_data.mat');

speedLevels = speedLevels(:)';

if ~exist('subjectNum', 'var') || isempty(subjectNum)
    subjectNum = size(thresholdsExpHig, 1);
end

selectedSpeeds = [4, 16];
[isFound, selectedIdx] = ismember(selectedSpeeds, speedLevels);

if any(~isFound)
    missingSpeeds = selectedSpeeds(~isFound);
    error('The following selected speed levels were not found in speedLevels: %s', ...
        num2str(missingSpeeds));
end

rng(1);  % reproducible bootstrap confidence intervals

%% Alias variables

pseExpPrevSlow = thresholdsExpHig;
pseConPrevSlow = thresholdsConHig;
pseExpPrevFast = thresholdsExpLow;
pseConPrevFast = thresholdsConLow;

p75ExpPrevSlow = disRangeExpHig;
p75ConPrevSlow = disRangeConHig;
p75ExpPrevFast = disRangeExpLow;
p75ConPrevFast = disRangeConLow;

outlierCriterionPSE = 5;

[pseExpPrevSlow, outlierNumPSE_ExpPrevSlow] = removeOutliersMAD(pseExpPrevSlow, outlierCriterionPSE);
[pseConPrevSlow, outlierNumPSE_ConPrevSlow] = removeOutliersMAD(pseConPrevSlow, outlierCriterionPSE);
[pseExpPrevFast, outlierNumPSE_ExpPrevFast] = removeOutliersMAD(pseExpPrevFast, outlierCriterionPSE);
[pseConPrevFast, outlierNumPSE_ConPrevFast] = removeOutliersMAD(pseConPrevFast, outlierCriterionPSE);

%% Calculate Weber fractions and sequential-shift terms

weberFractionExpPrevSlow = (p75ExpPrevSlow - pseExpPrevSlow) ./ pseExpPrevSlow;
weberFractionConPrevSlow = (p75ConPrevSlow - pseConPrevSlow) ./ pseConPrevSlow;
weberFractionExpPrevFast = (p75ExpPrevFast - pseExpPrevFast) ./ pseExpPrevFast;
weberFractionConPrevFast = (p75ConPrevFast - pseConPrevFast) ./ pseConPrevFast;

outlierCriterionWF = 3;

[weberFractionExpPrevSlow, outlierNumWF_ExpPrevSlow] = removeOutliersMAD(weberFractionExpPrevSlow, outlierCriterionWF);
[weberFractionConPrevSlow, outlierNumWF_ConPrevSlow] = removeOutliersMAD(weberFractionConPrevSlow, outlierCriterionWF);
[weberFractionExpPrevFast, outlierNumWF_ExpPrevFast] = removeOutliersMAD(weberFractionExpPrevFast, outlierCriterionWF);
[weberFractionConPrevFast, outlierNumWF_ConPrevFast] = removeOutliersMAD(weberFractionConPrevFast, outlierCriterionWF);

seqShiftExpWF = weberFractionExpPrevSlow - weberFractionExpPrevFast;
seqShiftConWF = weberFractionConPrevSlow - weberFractionConPrevFast;

%% Plot results

figure('Color', [1, 1, 1], ...
       'Position', [300, 250, 620, 300], ...
       'Name', 'Weber fraction sequential shift, selected speeds');

tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

titleStrings = {'Speed = 4 ˚/s', 'Speed = 16 ˚/s'};

for iPanel = 1:numel(selectedIdx)
    nexttile;
    hold on

    iCond = selectedIdx(iPanel);
    x = seqShiftExpWF(:, iCond);
    y = seqShiftConWF(:, iCond);

    validPair = ~isnan(x) & ~isnan(y);
    xValid = x(~isnan(x));
    yValid = y(~isnan(y));

    % Paired lines.
    if any(validPair)
        plot([ones(1, sum(validPair)) - 0.2; ones(1, sum(validPair)) + 0.2], ...
             [x(validPair)'; y(validPair)'], '-', ...
             'Color', [0.8, 0.8, 0.8]);
    end

    % Opaque half-violins.
    if ~isempty(xValid)
        drawViolinOrScatter(xValid, 1 - 0.3, [1, 0, 0], 'left');
        plot(repmat(1 - 0.2, numel(xValid), 1), xValid, ...
            'o', 'MarkerSize', 6, 'MarkerFaceColor', [1, 0, 0], ...
            'MarkerEdgeColor', [1, 1, 1]);
        plotMedianCI(xValid, 1 - 0.56, [1, 0, 0]);
    end

    if ~isempty(yValid)
        drawViolinOrScatter(yValid, 1 + 0.3, [0, 0, 1], 'right');
        plot(repmat(1 + 0.2, numel(yValid), 1), yValid, ...
            'o', 'MarkerSize', 6, 'MarkerFaceColor', [0, 0, 1], ...
            'MarkerEdgeColor', [1, 1, 1]);
        plotMedianCI(yValid, 1 + 0.56, [0, 0, 1]);
    end

    yline(0, ':k');

    set(gca, 'TickDir', 'out', 'Box', 'off', ...
        'XTick', [0.7, 1.3], ...
        'XTickLabel', {'Expansion', 'Contraction'}, ...
        'FontSize', 12);

    title(titleStrings{iPanel});
    ylabel('\Delta Weber fraction');
    xlim([0.3, 1.7]);
    ylim([-0.6, 0.6]);
end

%% Calculate only the p-values reported in the manuscript
% These p-values test whether the sequential-shift term in Weber fraction
% differs between expansion and contraction at 4 and 16 deg/s.
% Sequential shift is defined as previous-slower minus previous-faster.

pWF_seqShift = nan(1, numel(selectedIdx));

for iPanel = 1:numel(selectedIdx)
    iCond = selectedIdx(iPanel);
    x = seqShiftExpWF(:, iCond);
    y = seqShiftConWF(:, iCond);

    validPair = ~isnan(x) & ~isnan(y);

    if sum(validPair) >= 3
        pWF_seqShift(iPanel) = signrank(x(validPair), y(validPair));
    end
end

fprintf('Delta Weber fraction, expansion vs contraction, 4 deg/s: p = %.4g\n', pWF_seqShift(1));
fprintf('Delta Weber fraction, expansion vs contraction, 16 deg/s: p = %.4g\n', pWF_seqShift(2));
end

%% Helper functions

function [X, outlierNum] = removeOutliersMAD(X, k)
    if nargin < 2
        k = 3;
    end

    outlierNum = zeros(1, size(X, 2));

    for c = 1:size(X, 2)
        x = X(:, c);
        valid = ~isnan(x);

        if sum(valid) < 3
            X(:, c) = x;
            continue
        end

        medx = median(x(valid), 'omitnan');
        madx = mad(x(valid), 1) * 1.4826;

        if madx == 0 || isnan(madx)
            X(:, c) = x;
            continue
        end

        robust_z = abs(x - medx) ./ madx;
        outlierMask = robust_z > k;

        x(outlierMask) = NaN;
        outlierNum(c) = sum(outlierMask, 'omitnan');

        X(:, c) = x;
    end
end

function statsTable = pairedStatsByCond(X, Y, speedLevels, labelName)
    nConds = numel(speedLevels);

    nValid      = nan(nConds, 1);
    medianX     = nan(nConds, 1);
    medianY     = nan(nConds, 1);
    medianDiff  = nan(nConds, 1);
    pValue      = nan(nConds, 1);
    rRankBis    = nan(nConds, 1);
    cohenDz     = nan(nConds, 1);

    for i = 1:nConds
        x = X(:, i);
        y = Y(:, i);

        valid = ~isnan(x) & ~isnan(y);
        x = x(valid);
        y = y(valid);

        nValid(i) = numel(x);

        if nValid(i) >= 1
            d = x - y;
            medianX(i)    = median(x, 'omitnan');
            medianY(i)    = median(y, 'omitnan');
            medianDiff(i) = median(d, 'omitnan');
        end

        if nValid(i) >= 3
            pValue(i) = signrank(x, y);
            rRankBis(i) = pairedRankBiserial(x, y);

            d = x - y;
            if std(d, 'omitnan') > 0
                cohenDz(i) = mean(d, 'omitnan') ./ std(d, 'omitnan');
            end
        end
    end

    [hFDR, pFDR] = bhFDR(pValue, 0.05);

    statsTable = table( ...
        repmat({labelName}, nConds, 1), ...
        speedLevels(:), nValid, medianX, medianY, medianDiff, ...
        pValue, pFDR(:), hFDR(:), rRankBis, cohenDz, ...
        'VariableNames', {'Analysis', 'Speed', 'N', ...
        'MedianExp', 'MedianCon', 'MedianDiff_ExpMinusCon', ...
        'p', 'p_FDR_BH', 'h_FDR_BH', 'rankBiserial', 'cohenDz'});
end

function r = pairedRankBiserial(x, y)
    d = x - y;
    d = d(~isnan(d) & d ~= 0);

    if isempty(d)
        r = NaN;
        return
    end

    ranks = tiedrank(abs(d));
    Wpos = sum(ranks(d > 0));
    Wneg = sum(ranks(d < 0));

    r = (Wpos - Wneg) / (Wpos + Wneg);
end

function [hFDR, pAdj] = bhFDR(p, q)
    if nargin < 2
        q = 0.05;
    end

    originalSize = size(p);
    pVec = p(:);

    pAdjVec = nan(size(pVec));
    hVec = false(size(pVec));

    valid = ~isnan(pVec);
    pValid = pVec(valid);

    if isempty(pValid)
        hFDR = reshape(hVec, originalSize);
        pAdj = reshape(pAdjVec, originalSize);
        return
    end

    [pSorted, sortIdx] = sort(pValid);
    m = numel(pSorted);

    adjSorted = pSorted .* m ./ (1:m)';
    adjSorted = flipud(cummin(flipud(adjSorted)));
    adjSorted = min(adjSorted, 1);

    validIdx = find(valid);
    pAdjVec(validIdx(sortIdx)) = adjSorted;

    hVec = pAdjVec <= q;

    hFDR = reshape(hVec, originalSize);
    pAdj = reshape(pAdjVec, originalSize);
end

function drawViolinOrScatter(data, xPosition, colorValue, halfSide)
    if exist('Violin', 'file') == 2
        Violin({data}, xPosition, ...
            'ShowData', false, ...
            'EdgeColor', [1, 1, 1], ...
            'BoxWidth', 0.02, ...
            'BoxColor', [0.3, 0.3, 0.3], ...
            'ViolinColor', {colorValue}, ...
            'Width', 0.2, ...
            'ViolinAlpha', {1}, ...
            'HalfViolin', halfSide);
    else
        scatter(repmat(xPosition, size(data)), data, 25, ...
            'filled', 'MarkerFaceColor', colorValue, ...
            'MarkerEdgeColor', [1, 1, 1]);
    end
end

function plotMedianCI(data, xPosition, colorValue)
    data = data(~isnan(data));

    if isempty(data)
        return
    end

    medianValue = median(data, 'omitnan');

    if numel(data) >= 2
        ci = bootstrapMedianCI(data, 10000);
        yNeg = medianValue - ci(1);
        yPos = ci(2) - medianValue;
        errorbar(xPosition, medianValue, yNeg, yPos, ...
            'o', ...
            'MarkerSize', 4, ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerEdgeColor', colorValue, ...
            'Color', colorValue, ...
            'LineWidth', 2, ...
            'CapSize', 6);
    else
        plot(xPosition, medianValue, ...
            'o', ...
            'MarkerSize', 4, ...
            'MarkerFaceColor', [1, 1, 1], ...
            'MarkerEdgeColor', colorValue, ...
            'Color', colorValue, ...
            'LineWidth', 2);
    end
end

function ci = bootstrapMedianCI(data, nBoot)
    data = data(~isnan(data));
    n = numel(data);

    if n < 2
        ci = [NaN, NaN];
        return
    end

    bootMedian = nan(nBoot, 1);
    for b = 1:nBoot
        idx = randi(n, n, 1);
        bootMedian(b) = median(data(idx), 'omitnan');
    end

    ci = prctile(bootMedian, [2.5, 97.5]);
end
