function PFAnalysis_2D_populational
% This code aims to replicate the population analysis of Experiment 1 (Figure 3) in
% the following article:
% 
% Luo, J., Ashida, H., Yokoi, I., & Takemura, H. 
% The anisotropic sensitivities of perceptual speed between expanding and
% contracting optic flows.
% 
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)

clear; close all;

%% Set path
addpath(genpath('../Violinplot-Matlab-master'));
addpath(genpath('../boxplotGroup'));

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
for subjectIndex = 1:subjectNum
    fileList = dir(strcat('*', subjectNames{subjectIndex}, '_*.mat'));
    speedsTest = [];
    speedsComp = [];
    responses = [];
    for i = 1:length(fileList)
        load(fileList(i).name);
        speedsTest = [speedsTest, taskData.speedsTestDone]; %#ok<AGROW>
        speedsComp = [speedsComp, round(taskData.speedsCompDone, 4)]; %#ok<AGROW>
        responses = [responses, taskData.responses]; %#ok<AGROW>
    end

    speedConditions = taskData.speedConditions;
    speedCondNum = length(speedConditions)/2;

    speedsTestConCond = [];
    speedsCompConCond = [];
    responsesConCond = [];
    speedsTestExpCond = [];
    speedsCompExpCond = [];
    responsesExpCond = [];
    for i = 1:length(speedConditions)
        if speedConditions(i)<0
            speedsTestConCond = [speedsTestConCond; speedsTest(speedsTest==speedConditions(i))]; %#ok<AGROW>
            speedsCompConCond = [speedsCompConCond; speedsComp(speedsTest==speedConditions(i))]; %#ok<AGROW>
            responsesConCond = [responsesConCond; responses(speedsTest==speedConditions(i))]; %#ok<AGROW>
        else
            speedsTestExpCond = [speedsTestExpCond; speedsTest(speedsTest==speedConditions(i))]; %#ok<AGROW>
            speedsCompExpCond = [speedsCompExpCond; speedsComp(speedsTest==speedConditions(i))]; %#ok<AGROW>
            responsesExpCond = [responsesExpCond; responses(speedsTest==speedConditions(i))]; %#ok<AGROW>
        end
    end

    for i = 1:size(speedsTestExpCond, 1)
        [stiLevelsConG(i, :), numPosConG(i, :), numTotalConG(i, :)] = PAL_PFML_GroupTrialsbyX(speedsCompConCond(i, :), responsesConCond(i, :), ones(1, size(speedsCompConCond, 2))); %#ok<SAGROW>
        stiLevelsConG(i, :) = -stiLevelsConG(i, :); %#ok<SAGROW>
        [stiLevelsExpG(i, :), numPosExpG(i, :), numTotalExpG(i, :)] = PAL_PFML_GroupTrialsbyX(speedsCompExpCond(i, :), responsesExpCond(i, :), ones(1, size(speedsCompExpCond, 2))); %#ok<SAGROW>
    end

    stiLevelsConG = log2(stiLevelsConG);
    stiLevelsExpG = log2(stiLevelsExpG);

    for n = 1:size(stiLevelsExpG, 1)
        PF = @PAL_CumulativeNormal;
        paramsFree = [1 1 0 0];
        searchGrid.alpha = min(stiLevelsExpG(n, :)):0.01:max(stiLevelsExpG(n, :));
        searchGrid.beta = logspace(-1, 1, 100);
        searchGrid.gamma = 0.00;
        searchGrid.lambda = 0.00;
        [paramsValues, LL, exitFlag] = PAL_PFML_Fit(stiLevelsExpG(n, :), numPosExpG(n, :), numTotalExpG(n, :), searchGrid, paramsFree, PF);
        thresholdsExp(subjectIndex, n) = 2^PF(paramsValues(1, :), 0.5, 'Inverse'); %#ok<SAGROW>
        disRangeExp(subjectIndex, n) = 2^PF(paramsValues(1, :), 0.75, 'Inverse'); %#ok<SAGROW>
    end

    for n = 1:size(stiLevelsConG, 1)
        PF = @PAL_CumulativeNormal;
        paramsFree = [1 1 0 0];
        searchGrid.alpha = min(stiLevelsConG(n, :)):0.01:max(stiLevelsConG(n, :));
        searchGrid.beta = logspace(-1, 1, 100);
        searchGrid.gamma = 0.00;
        searchGrid.lambda = 0.00;
        [paramsValues, LL, exitFlag] = PAL_PFML_Fit(stiLevelsConG(n, :), numPosConG(n, :), numTotalConG(n, :), searchGrid, paramsFree, PF);
        thresholdsCon(subjectIndex, n) = 2^PF(paramsValues(1, :), 0.5, 'Inverse'); %#ok<SAGROW>
        disRangeCon(subjectIndex, n) = 2^PF(paramsValues(1, :), 0.75, 'Inverse'); %#ok<SAGROW>
    end
end

speedConds = speedConditions(1:length(speedConditions)/2);

%% Weber fraction results

weberFractionExp = (disRangeExp - thresholdsExp)./thresholdsExp;
weberFractionCon = (disRangeCon - thresholdsCon)./thresholdsCon;

% Summary Statistics on Weber Fraction
medianWeberFractionExp = median(weberFractionExp, 1);
medianWeberFractionCon = median(weberFractionCon, 1);

% Calculate confidence interval by bootstrapping
CI_exp = nan(2, length(speedConds));
CI_con = nan(2, length(speedConds));
for n = 1:length(speedConds)
    CI_exp(:,n) = bootci(10000, {@median, weberFractionExp(:,n)}, 'type', 'bca');
    CI_con(:,n) = bootci(10000, {@median, weberFractionCon(:,n)}, 'type', 'bca');
end

% Difference from median
yneg_exp = medianWeberFractionExp - CI_exp(1,:);
ypos_exp = CI_exp(2,:) - medianWeberFractionExp;
yneg_con = medianWeberFractionCon - CI_con(1,:);
ypos_con = CI_con(2,:) - medianWeberFractionCon;

% Plot results
figure('Color', [1, 1, 1], 'Position', [300, 200, 1000, 600], 'Name', sprintf('Number of Participants: %d', subjectNum))
yUpperLim = [1.4, 1.2, 1, 0.8, 0.8];
tiledlayout(2, round(length(medianWeberFractionExp)/2), 'Padding','compact', 'TileSpacing','loose');
for i = 1:length(medianWeberFractionExp)
    nexttile;
    hold on
    plot([ones(1, subjectNum) - 0.2; ones(1, subjectNum) + 0.2], [weberFractionExp(:, i)'; weberFractionCon(:, i)'], '-', 'Color', [0.8, 0.8, 0.8]);
    Violin({weberFractionExp(:, i)}, 1 - 0.3,...
        'ShowData', false, 'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, 'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[1, 0, 0]}, 'Width', 0.2, 'ViolinAlpha', {1}, 'HalfViolin','left');
    Violin({weberFractionCon(:, i)}, 1 + 0.3,...
        'ShowData', false, 'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, 'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[0, 0, 1]}, 'Width', 0.2, 'ViolinAlpha', {1}, 'HalfViolin','right');
    plot(repmat(1 - 0.2, subjectNum, 1), weberFractionExp(:, i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [1, 1, 1]);
    plot(repmat(1 + 0.2, subjectNum, 1), weberFractionCon(:, i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', [1, 1, 1]);
    errorbar(1 - 0.56, medianWeberFractionExp(i), yneg_exp(i), ypos_exp(i), ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', [1, 0, 0], 'Color', 'r', 'LineWidth', 2);
    errorbar(1 + 0.56, medianWeberFractionCon(i), yneg_con(i), ypos_con(i), ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', [0, 0, 1], 'Color', 'b', 'LineWidth', 2);

    set(gca, 'TickDir', 'out', 'box', 'off', 'XTick', [0.7, 1.3], 'XTickLabel', {'Expansion', 'Contraction'})
    title(strcat({'Speed ='}, {' '}, {num2str(speedConditions(i))}, {' '},  {'°/s'}))
    ylabel('Weber Fraction')
    xlim([0.3, 1.7])
    ylim([0, yUpperLim(i)])

    if i == 1
        text(0.35, 1.34, ['N = ', num2str(subjectNum)]);
    end
end

nexttile;
hold on
boxplotGroup({weberFractionExp, weberFractionCon}, 'Colors', 'rb', 'Symbol', '.','GroupType', 'betweenGroups', 'notch', 'on');
hw = findobj(gca, 'Type', 'Line', 'LineStyle', '--');
set(hw, 'LineStyle', ':');
set(gca, 'TickDir', 'out', 'box', 'off', 'XTick', [1.5, 4.5, 7.5, 10.5, 13.5], 'XTickLabel', compose('%g', speedConditions(1:i)));
xlabel('Speed conditions (deg/s)')
ylabel('Weber Fraction');
ylim([0, max(yUpperLim)])

pExp = plot(nan, nan, 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [1, 0, 0], ...
            'MarkerSize', 10, 'DisplayName', 'Expansion');
pCon = plot(nan, nan, 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0, 0, 1], ...
            'MarkerSize', 10, 'DisplayName', 'Contraction');

lgd = legend([pExp pCon], 'Location', 'northeast', 'Box', 'off');
lgd.ItemTokenSize = [20 10];

ax = findall(gcf,'type','axes');
ax = ax(arrayfun(@(h) isprop(h,'YLabel'), ax));
pos = cat(1, ax.Position);
lefts = pos(:,1);
tol = 1e-3;
[uniqueLefts, ~, colIdx] = uniquetol(lefts, tol);
margin = 0.035;
targetFigX = uniqueLefts - margin;
for k = 1:numel(ax)
    a = ax(k);
    a.YLabel.Units = 'normalized';
    x_norm = (targetFigX(colIdx(k)) - a.Position(1)) / a.Position(3);
    p = a.YLabel.Position;
    p(1) = x_norm;
    p(2) = 0.5;
    a.YLabel.Position = p;
    a.YLabel.Clipping = 'off';
end

%% Statistics

deltaWeb = weberFractionExp - weberFractionCon;

nConds = length(speedConds);
pNPWeb = nan(1,nConds);
hNPWeb = nan(1,nConds);
r_rb   = nan(1,nConds);
effectSizeWeb = nan(1,nConds);
for cond = 1:nConds

    % Nonparametric test (paired)
    [pNPWeb(cond), hNPWeb(cond)] = signrank(weberFractionExp(:, cond), weberFractionCon(:, cond));

    % Rank-biserial correlation (paired)
    d = deltaWeb(:, cond);
    d = d(d ~= 0);

    if isempty(d)
        r_rb(cond) = NaN;
    else
        ranks = tiedrank(abs(d));
        Wpos = sum(ranks(d > 0));
        Wneg = sum(ranks(d < 0));
        r_rb(cond) = (Wpos - Wneg) / (Wpos + Wneg);
    end

    % Robust Cohen's d
    effectSizeWebTable = meanEffectSize(weberFractionExp(:, cond), weberFractionCon(:, cond), Paired = true, Effect = "robustcohen"); 
    effectSizeWeb(cond) = effectSizeWebTable.Effect;
end

% FDR correction (Benjamini-Hochberg)
p = pNPWeb(:);
m = numel(p);
q = 0.05;

% Sort p-values
[p_sorted, idx] = sort(p);
rank = (1:m)';

% BH adjustment
thresholds = (1:m) / m * q;
h_sorted = p_sorted <= transpose(thresholds);

% Restore original order
hFDR(idx) = h_sorted;
