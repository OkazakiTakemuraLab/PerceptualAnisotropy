function PFAnalysis_3D_populational
% This code aims to replicate the population analysis of Experiment 2 (Figure 6) in
% the following article:
% 
% Luo, J., Ashida, H., Yokoi, I. & Takemura, H. (2026) Anisotropic speed discrimination sensitivity for expanding and contracting optic flow. Journal of Vision, 26(10):18, 1–17.
% 
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)

clear; close all;

%% Set path
addpath(genpath('../Violinplot-Matlab-master'));
addpath(genpath('../boxplotGroup'));
addpath(genpath('../Palamedes'));

%% Load data and calculate CP(0.5) and CP(0.75)

fileInformation = dir('OpticalFlow*.mat');
subjectNum = length(fileInformation);

for n = 1:subjectNum
    load(fullfile(fileInformation(n).folder, fileInformation(n).name));

    speedsTest = taskData.speedConditions/1000; % m/s
    speedsComp = cat(1, taskData.staircaseResults.x);
    speedsComp = speedsComp(:, 1:end-1); % the final values are not used in the task
    
    responses = cat(1, taskData.staircaseResults.response);

    PF = @PAL_CumulativeNormal;
    paramsFree = [1 1 0 0];
    searchGrid.beta = logspace(-1, 1, 100);
    searchGrid.gamma = 0.00;
    searchGrid.lambda = 0.00;
    fittingOptions = PAL_minimize('options');
    for i = 1:length(speedsTest)
        [stiLevelsG, numPosG, numTotalG] =...
            PAL_PFML_GroupTrialsbyX(speedsComp(i, :), responses(i, :), ones(1, size(speedsComp, 2)));

        searchGrid.alpha = min(stiLevelsG):0.01:max(stiLevelsG);
        [paramsValues, LL, exitFlag] =...
            PAL_PFML_Fit(stiLevelsG, numPosG, numTotalG, searchGrid, paramsFree, PF, 'searchOptions', fittingOptions);

        if i <= length(speedsTest)/2
            thresholdsExp(n, i) = PF(paramsValues, 0.5, 'Inverse'); %#ok<SAGROW>
            disRangeExp(n, i) = PF(paramsValues, 0.75, 'Inverse'); %#ok<SAGROW>
        else
            thresholdsCon(n, i - length(speedsTest)/2) = PF(paramsValues, 0.5, 'Inverse'); %#ok<SAGROW>
            disRangeCon(n, i - length(speedsTest)/2) = PF(paramsValues, 0.75, 'Inverse'); %#ok<SAGROW>
        end
    end
end

disRangeExp = 2.^disRangeExp;
thresholdsExp = 2.^thresholdsExp;
disRangeCon = 2.^disRangeCon;
thresholdsCon = 2.^thresholdsCon;

speedConds = speedsTest(1:length(speedsTest)/2);

%% Weber fraction results

weberFractionExp = (disRangeExp - thresholdsExp)./thresholdsExp; 
weberFractionCon = (disRangeCon - thresholdsCon)./thresholdsCon;

% Summary Statistics on Weber Fraction
medianWeberFractionExp = median(weberFractionExp, 1);
medianWeberFractionCon = median(weberFractionCon, 1);

% Calculate confidence interval by bootstrapping
CI_exp = nan(2, length(speedConds));
CI_con = nan(2, length(speedConds));
rng(1); % Fix the random seed to ensure reproducibility of bootstrap results
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
figure('Color', [1, 1, 1], 'Position', [300, 200, 1000, 600], 'Name', sprintf('Number of Participants: %d', subjectNum));
yUpperLim = [4, 2, 3, 0.8, 0.8];
tiledlayout(2, round(length(speedsTest)/4) + 1, 'Padding','compact', 'TileSpacing','loose');
for i = 1:length(speedsTest)/2
    nexttile;
    hold on
    plot([ones(1, subjectNum) - 0.2; ones(1, subjectNum) + 0.2], [weberFractionExp(:, i)'; weberFractionCon(:, i)'], '-', 'Color', [0.8, 0.8, 0.8]);
    Violin({weberFractionExp(:, i)}, 1 - 0.3, 'ShowData', false, 'ShowMean', false, 'ShowMedian', true, ...
        'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, 'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[1, 0, 0]}, 'Width', 0.2, 'ViolinAlpha', {1}, 'HalfViolin','left');
    Violin({weberFractionCon(:, i)}, 1 + 0.3, 'ShowData', false, 'ShowMean', false, 'ShowMedian', true,...
        'EdgeColor', [1, 1, 1], 'BoxWidth', 0.02, 'BoxColor', [0.3, 0.3, 0.3], 'ViolinColor', {[0, 0, 1]}, 'Width', 0.2, 'ViolinAlpha', {1}, 'HalfViolin','right');
    plot(repmat(1 - 0.2, subjectNum, 1), weberFractionExp(:, i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [1, 1, 1]);
    plot(repmat(1 + 0.2, subjectNum, 1), weberFractionCon(:, i), 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', [1, 1, 1]);
    errorbar(1 - 0.56, medianWeberFractionExp(i), yneg_exp(i), ypos_exp(i), ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', [1, 0, 0], 'Color', 'r', 'LineWidth', 2);
    errorbar(1 + 0.56, medianWeberFractionCon(i), yneg_con(i), ypos_con(i), ...
        'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', [1, 1, 1], 'MarkerEdgeColor', [0, 0, 1], 'Color', 'b', 'LineWidth', 2);

    set(gca, 'TickDir', 'out', 'box', 'off', 'XTick', [0.7, 1.3], 'XTickLabel', {'Exp', 'Con'})
    title(strcat({'Speed ='}, {' '}, {num2str(speedsTest(i))}, {' '},  {'m/s'}))
    if i == 1
        xlabel('Opic Flow Patterns')
    end
    ylabel('Weber Fraction')
    xlim([0.3, 1.7])
    ylim([0, yUpperLim(i)])

    if i == 1
        text(0.35, 3.85, ['N = ', num2str(subjectNum)]);
    end
end

nexttile;
hold on
boxplotGroup({weberFractionExp, weberFractionCon}, 'Colors', 'rb', 'Symbol', '.','GroupType', 'betweenGroups', 'notch', 'on');
hw = findobj(gca, 'Type', 'Line', 'LineStyle', '--');
set(hw, 'LineStyle', ':');
set(gca, 'TickDir', 'out', 'box', 'off', 'XTick', [1.5, 4.5, 7.5, 10.5], 'XTickLabel', compose('%g', speedsTest(1:i)));
xlabel('Speed conditions (m/s)')
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