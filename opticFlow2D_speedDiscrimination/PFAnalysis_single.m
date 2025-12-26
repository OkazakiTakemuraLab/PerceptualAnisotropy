function PFAnalysis_single

% Replicate analysis of single subject data by fitting psychometric
% functions to data acquired in Experiment 1 of the following paper:
% 
% Luo, J., Ashida, H., Yokoi, I., & Takemura, H. 
% The anisotropic sensitivities of perceptual speed between expanding and
% contracting optic flows.
%
% To replicate Figure 2 of Luo et al., please select .mat file including
% the name of "P12".
% Dependency: Palamedes Toolbox (https://www.palamedestoolbox.org/)

clear; close all;

%% Load data

% Select a data file (.mat) with participant ID first using user interface.
% Then the code automatically select all data files from the same participant.  
fileName = uigetfile('*.mat');
subjectName = extractBetween(fileName, '_', '_');
fileList = dir(strcat('*', subjectName{1}, '*.mat'));

speedsRefe = [];
speedsComp = [];
responses = [];
for i = 1:size(fileList, 1)
    load(fileList(i).name);
    speedsRefe = [speedsRefe, taskData.speedsTestDone]; %#ok<AGROW>
    speedsComp = [speedsComp, round(taskData.speedsCompDone, 4)]; %#ok<AGROW>
    responses = [responses, taskData.responses]; %#ok<AGROW>
end
speedConditions = taskData.speedConditions;
speedCondNum = length(speedConditions)/2;

speedsRefeConCond = [];
speedsCompConCond = [];
responsesConCond = [];
speedsRefeExpCond = [];
speedsCompExpCond = [];
responsesExpCond = [];
for i = 1:length(speedConditions)
    if speedConditions(i)<0
        speedsRefeConCond = [speedsRefeConCond; speedsRefe(speedsRefe==speedConditions(i))]; %#ok<AGROW>
        speedsCompConCond = [speedsCompConCond; speedsComp(speedsRefe==speedConditions(i))]; %#ok<AGROW>
        responsesConCond = [responsesConCond; responses(speedsRefe==speedConditions(i))]; %#ok<AGROW>
    else
        speedsRefeExpCond = [speedsRefeExpCond; speedsRefe(speedsRefe==speedConditions(i))]; %#ok<AGROW>
        speedsCompExpCond = [speedsCompExpCond; speedsComp(speedsRefe==speedConditions(i))]; %#ok<AGROW>
        responsesExpCond = [responsesExpCond; responses(speedsRefe==speedConditions(i))]; %#ok<AGROW>
    end
end

for i = 1:size(speedsRefeExpCond, 1)
    [stiLevelsConG(i, :), numPosConG(i, :), numTotalConG(i, :)] =...
        PAL_PFML_GroupTrialsbyX(speedsCompConCond(i, :), responsesConCond(i, :), ones(1, size(speedsCompConCond, 2))); %#ok<SAGROW>
    stiLevelsConG(i, :) = -stiLevelsConG(i, :); %#ok<SAGROW>
    [stiLevelsExpG(i, :), numPosExpG(i, :), numTotalExpG(i, :)] =...
        PAL_PFML_GroupTrialsbyX(speedsCompExpCond(i, :), responsesExpCond(i, :), ones(1, size(speedsCompExpCond, 2))); %#ok<SAGROW>
end

stiLevelsConG = log2(stiLevelsConG);
stiLevelsExpG = log2(stiLevelsExpG);

%% Fit psychometric functions

PF = @PAL_CumulativeNormal;
paramsFree = [1 1 0 0];
searchGrid.beta = logspace(-1, 1, 100);
searchGrid.gamma = 0.00;
searchGrid.lambda = 0.00;
fittingOptions = PAL_minimize('options');

for n = 1:speedCondNum
    searchGrid.alpha = min(stiLevelsExpG(n, :)):0.01:max(stiLevelsExpG(n, :));
    [paramsValuesExp(n, :), LLExp(n, :), exitFlagExp(n, :)] =...
        PAL_PFML_Fit(stiLevelsExpG(n, :), numPosExpG(n, :), numTotalExpG(n, :), searchGrid, paramsFree, PF, 'searchOptions', fittingOptions); %#ok<SAGROW>
    thresholdsExp(n) = 2^PF(paramsValuesExp(n, :), 0.5, 'Inverse'); %#ok<SAGROW>
    disRangeExp(n) = 2^PF(paramsValuesExp(n, :), 0.75, 'Inverse'); %#ok<SAGROW>


    searchGrid.alpha = min(stiLevelsConG(n, :)):0.01:max(stiLevelsConG(n, :));
    [paramsValuesCon(n, :), LLCon(n, :), exitFlagCon(n, :)] =...
        PAL_PFML_Fit(stiLevelsConG(n, :), numPosConG(n, :), numTotalConG(n, :), searchGrid, paramsFree, PF, 'searchOptions', fittingOptions); %#ok<SAGROW>
    thresholdsCon(n) = 2^PF(paramsValuesCon(n, :), 0.5, 'Inverse'); %#ok<SAGROW>
    disRangeCon(n) = 2^PF(paramsValuesCon(n, :), 0.75, 'Inverse'); %#ok<SAGROW>
end

%% Plot psychometric functions and Weber fractions

figure('Color', [1, 1, 1], 'Position', [300, 300, 900, 400]);
t = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'compact'); 

markerShape = {'v', 'square', 'o', 'diamond', '^'};
legendSpeeds = compose('%g', speedConditions(1:speedCondNum));
legends = cell(1, 2*speedCondNum);
legends(2:2:end) = legendSpeeds;
legends(1:2:end) = {''};

stiLevelFineExp = min(stiLevelsExpG, [], 'all'):0.01:max(stiLevelsExpG, [], 'all');
stiLevelFineCon = min(stiLevelsConG, [], 'all'):0.01:max(stiLevelsConG, [], 'all');

% --- Top plot (expansion) ---
nexttile(1, [1, 2]);
hold on;
for n = 1:speedCondNum
    fitPsyFuncExp = PF(paramsValuesExp(n, :), stiLevelFineExp);
    plot(stiLevelFineExp, fitPsyFuncExp, '-', 'Color', [1, 0, 0], 'LineWidth', 1);
    plot(stiLevelsExpG(n, :), numPosExpG(n, :)./numTotalExpG(n, :), markerShape{n}, ...
        'MarkerSize', 10, 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 0.5);
end
yline(0.5, ':k'); yline(0.75, ':k');
text(-0.8, 0.95, 'Expansion', 'FontSize', 12, 'FontWeight', 'bold');
text(-1.65, 0.75, '0.75');
lgdExp = legend(legends, 'Box', 'off');
lgdExp.Position = lgdExp.Position + [0.04, -0.11, -0.02, 0];
set(gca, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 0.5, ...
    'FontName', 'Arial', 'FontSize', 10);

% --- Bottom plot (contraction) ---
nexttile(4, [1, 2]);
hold on;
for n = 1:speedCondNum
    fitPsyFuncCon = PF(paramsValuesCon(n, :), stiLevelFineCon);
    plot(stiLevelFineCon, fitPsyFuncCon, '-', 'Color', [0, 0, 1], 'LineWidth', 1);
    plot(stiLevelsConG(n, :), numPosConG(n, :)./numTotalConG(n, :), markerShape{n}, ...
        'MarkerSize', 10, 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', [1, 1, 1], 'LineWidth', 0.5);
end
yline(0.5, ':k'); yline(0.75, ':k');
xlabel('Speed of Comparison Stimulus [log_2(°/s)]')
text(-0.8, 0.95, 'Contraction', 'FontSize', 12, 'FontWeight', 'bold');
text(-1.65, 0.75, '0.75');
lgdCon = legend(legends, 'Box', 'off');
lgdCon.Position = lgdCon.Position + [0.04, -0.1, -0.02, 0];
set(gca, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on', 'LineWidth', 0.5, ...
    'FontName', 'Arial', 'FontSize', 10);

ylabel(t, sprintf('Choice Proportion\n(Comparison > Reference)'));

% --- Left plot (Weber fraction) ---
weberFractionExp = (disRangeExp - thresholdsExp)./thresholdsExp;
weberFractionCon = (disRangeCon - thresholdsCon)./thresholdsCon;
nexttile(3, [2 1]);
hold on;
bar((1:speedCondNum) - 0.18, weberFractionExp, 'EdgeColor', 'none', 'FaceColor', [1, 0, 0], 'BarWidth', 0.3);
bar((1:speedCondNum) + 0.18, weberFractionCon, 'EdgeColor', 'none', 'FaceColor', [0, 0, 1], 'BarWidth', 0.3);

legend({'Expansion', 'Contraction'}, 'Location', 'northeast', 'Box', 'off');
set(gca, 'TickDir', 'out', 'XTick', 0:1:6, 'XTickLabel', {'', '1', '4', '16', '64', '128', ''}, 'FontName', 'Arial', 'FontSize', 10);
xlim([0, 6]);
xlabel('Speed conditions (°/s)');
ylabel('Weber Fraction');