function PFAnalyisis_3D_single

% Replicate analysis of single subject data by fitting psychometric
% functions `to data acquired in Experiment 2.
% To replicate Figure 5 of Luo et al., please select .mat file including
% the name of "P13".
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)

clear;

%% Load data

[fileName, filePath] = uigetfile('*.mat');
load(fullfile(filePath, fileName));

speedsRefe = taskData.speedConditions/1000; % m/s
speedCondNum = length(speedsRefe)/2;
speedsComp = cat(1, taskData.staircaseResults.x);
speedsComp = speedsComp(:, 1:end-1); % the final values are not used in the task

responses = cat(1, taskData.staircaseResults.response);

%% Fit psychometric functions

PF = @PAL_CumulativeNormal;
paramsFree = [1 1 0 0];
searchGrid.beta = logspace(-1, 1, 100);
searchGrid.gamma = 0.00;
searchGrid.lambda = 0.00;
fittingOptions = PAL_minimize('options');
for i = 1:length(speedsRefe)
    [stiLevelsG{i}, numPosG{i}, numTotalG{i}] =...
        PAL_PFML_GroupTrialsbyX(speedsComp(i, :), responses(i, :), ones(1, size(speedsComp, 2))); %#ok<SAGROW>

    searchGrid.alpha = min(stiLevelsG{i}):0.01:max(stiLevelsG{i});
    [paramsValues(i, :), LL(i, :), exitFlag(i, :)] =...
        PAL_PFML_Fit(stiLevelsG{i}, numPosG{i}, numTotalG{i}, searchGrid, paramsFree, PF, 'searchOptions', fittingOptions); %#ok<SAGROW>
end

%% Plot psychometric functions and Weber fractions

figure('Color', [1, 1, 1], 'Position', [300, 300, 900, 400]);
t = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

markerShape = {'v', 'square', 'o', 'diamond'};

alphaRange = [-2, 8]; % same for all conditions
stiLevelFine = min(alphaRange):0.01:max(alphaRange);

% --- Top plot (expansion) ---
nexttile(1, [1, 2]);
hold on
for i = 1:speedCondNum
    plot(nan, nan, markerShape{i}, 'MarkerFaceColor','r', 'MarkerEdgeColor','w', 'MarkerSize', 20, 'LineWidth', 0.5);
end
for i = 1:speedCondNum
    thresholdsExp(i) = PF(paramsValues(i, :), 0.5, 'Inverse'); %#ok<SAGROW>
    disRangeExp(i) = PF(paramsValues(i, :), 0.75, 'Inverse'); %#ok<SAGROW>
    fitPsyFunc = PF(paramsValues(i, :), stiLevelFine);
    plot(stiLevelFine, fitPsyFunc, '-', 'Color', [1, 0, 0], 'LineWidth', 1);

    for k = 1:length(stiLevelsG{i})
        plot(stiLevelsG{i}(k), numPosG{i}(k)/numTotalG{i}(k),...
            markerShape{i}, 'MarkerFaceColor','r', 'MarkerEdgeColor','w', 'MarkerSize', 20*sqrt(numTotalG{i}(k)/sum(numTotalG{i})), 'LineWidth', 0.5);
    end
end
yline(0.5, ':k'); yline(0.75, ':k');
lgdExp = legend(compose('%g', speedsRefe(1:speedCondNum)), 'Box', 'off');
lgdExp.Position = lgdExp.Position + [0.04, -0.14, -0.02, 0];
text(-1.8, 0.95, 'Expansion', 'FontSize', 12, 'FontWeight', 'bold');
text(-2.65, 0.75, '0.75');
set(gca, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 0.5, ...
    'FontName', 'Arial', 'FontSize', 10);

% --- Bottom plot (contraction) ---
nexttile(4, [1, 2]);
hold on
for i = 1:speedCondNum
    plot(nan, nan, markerShape{i}, 'MarkerFaceColor','b', 'MarkerEdgeColor','w', 'MarkerSize', 20, 'LineWidth', 0.5);
end
for i = (1:speedCondNum) + speedCondNum
    thresholdsCon(i - speedCondNum) = PF(paramsValues(i, :), 0.5, 'Inverse'); %#ok<SAGROW>
    disRangeCon(i - speedCondNum) = PF(paramsValues(i, :), 0.75, 'Inverse'); %#ok<SAGROW>
    fitPsyFunc = PF(paramsValues(i, :), stiLevelFine);
    plot(stiLevelFine, fitPsyFunc, '-', 'Color', [0, 0, 1], 'LineWidth', 1);

    for k = 1:length(stiLevelsG{i})
        plot(stiLevelsG{i}(k), numPosG{i}(k)/numTotalG{i}(k),...
            markerShape{i - speedCondNum}, 'MarkerFaceColor','b', 'MarkerEdgeColor','w', 'MarkerSize', 20*sqrt(numTotalG{i}(k)/sum(numTotalG{i})), 'LineWidth', 0.5);
    end
end
yline(0.5, ':k'); yline(0.75, ':k');
lgdCon = legend(compose('%g', speedsRefe(1:speedCondNum)), 'Box', 'off');
lgdCon.Position = lgdCon.Position + [0.04, -0.1, -0.02, 0];
xlabel('Speed of Comparison Stimulus [log_2(m/s)]');
text(-1.8, 0.95, 'Contraction', 'FontSize', 12, 'FontWeight', 'bold');
text(-2.65, 0.75, '0.75');
set(gca, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 0.5, ...
    'FontName', 'Arial', 'FontSize', 10);

ylabel(t, sprintf('Choice Proportion\n(Comparison > Reference)'));

% --- Left plot (Weber fraction) ---
disRangeExp = 2.^disRangeExp;
thresholdsExp = 2.^thresholdsExp;
disRangeCon = 2.^disRangeCon;
thresholdsCon = 2.^thresholdsCon;
weberFractionExp = (disRangeExp - thresholdsExp)./thresholdsExp;
weberFractionCon = (disRangeCon - thresholdsCon)./thresholdsCon;
nexttile(3, [2, 1]);
hold on
bar((1:speedCondNum) - 0.18, weberFractionExp, 'EdgeColor', 'none', 'FaceColor', [1, 0, 0], 'BarWidth', 0.3);
bar((1:speedCondNum) + 0.18, weberFractionCon, 'EdgeColor', 'none', 'FaceColor', [0, 0, 1], 'BarWidth', 0.3);

legend({'Expansion', 'Contraction'}, 'Location', 'northeast', 'Box', 'off');
set(gca, 'TickDir', 'out', 'XTick', 0:1:6, 'XTickLabel', {'', '0.5', '1.5', '20', '100', ''}, 'FontName', 'Arial', 'FontSize', 10);
xlim([0, 5]);
xlabel('Speed conditions (m/s)');
ylabel('Weber Fraction');
