%% CCA with MINDy Parameters
% 
% 
% Ruiqi Chen, 2023/07/13
% 
% The CCA analysis was performed using https://github.com/rq-Chen/HCP_CCA_1200_OSF 

clear; close all

cca_file = fullfile('data', 'HCP_CCA', ...
    'HCP_Rest_FIX_Simple_Mdl200_sess_W_average', 'hcp_1200_cca_MINDy_2023-11-05_T12-31-27.mat');
BaseFont = 14;
flSign = -1;  % Flip the sign for visualization

%% 
% 
%% CCA correlation
% 
% 
% Replicate Figure 1B&C in (Smith et al., 2015). We copied some codes from the 
% (Goyal et al., 2020) replication.

% Load data
load(cca_file, 'grotRpval', 'grotR', 'grot_NET', 'grot_VARS');
I = 1:20;  % CCA modes to show

f = figure('Units', 'inches', 'Position', [0 0.5 16 9]);
t = tiledlayout(6, 3);

% CCA correlation
nexttile(3, [2 1]);
hold on;
% plot(I, grotR(I), 'k', 'Marker', '.');
stem(I, grotR(I), 'k', 'filled')
for i = 1:length(I)
    if grotRpval(i) < 0.05
        text(i, grotR(i) * 1.03, '*', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
end
text(max(I) / 4, grotR(1), '*: p < .05 (permutation test)')
ylim([0.9 * min(grotR(I)), 1.1 * max(grotR(I))])
ylabel('Correlation')
fontsize(gca(), BaseFont, 'points');
title('Correlation between subject measures and MINDy CCA modes', ...
    'FontSize', BaseFont + 2, 'FontWeight', 'bold')

% Variance of MINDy parameters explained by CCA modes
nexttile(9, [2 1]);
hold on;
% Draw the rectangles for null distributions per mode
for i=1:length(I)
  rectangle('Position',[i-0.5 grot_NET(i,2) 1 grot_NET(i,4)-grot_NET(i,2)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
plot(grot_NET(I,3),'k'); plot(grot_NET(I,1),'b'); plot(grot_NET(I,1),'b.'); % plot(grot_NET(I,5),'g');  % turned off showing the PCA equivalent plots
yl = ylim();
text(max(I) / 6, yl(1) * 1/5 + yl(2) * 4/5, 'Shade: 5%-95% of null distribution')
ylabel('% of variance')
xlim([1 20])
% ylim([0.3 1.2]), ...
fontsize(gca(), BaseFont, 'points');
title('MINDy parameters variance explained', ...
    'FontSize', BaseFont + 2, 'FontWeight', 'bold')

% Variance of subject measures explained by CCA modes
nexttile(15, [2 1]);
hold on;
for i=1:length(I)
  rectangle('Position',[i-0.5 grot_VARS(i,2) 1 grot_VARS(i,4)-grot_VARS(i,2)],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
end
plot(grot_VARS(I,3),'k'); plot(grot_VARS(I,1),'b'); plot(grot_VARS(I,1),'b.'); % plot(grot_VARS(I,5),'g');
yl = ylim();
text(max(I) / 6, yl(1) * 1/5 + yl(2) * 4/5, 'Shade: 5%-95% of null distribution')
ylabel('% of variance')
xlim([1 20])
% ylim([0 2])
% yticks([0 0.5 1 1.5 2])
fontsize(gca(), BaseFont, 'points');
title('Subject measures variance explained', ...
    'FontSize', BaseFont + 2, 'FontWeight', 'bold')


%% 
% 
%% Most correlated behavioral measures with the CCA mode
% 
% 
% Load data:

load(cca_file, 'subs', 'SMs', 'grotBBd', 'varskeep');

% Flip sign
grotBBd = grotBBd * flSign;


%% 
% Rename several subject measures to match the dictionary:

SMs(strcmp(SMs, 'age')) = {'Age'};
SMs(strcmp(SMs, 'race')) = {'Race'};
SMs(strcmp(SMs, 'handedness')) = {'Handedness'};
SMs(strcmp(SMs, 'subject')) = {'Subject'};
%% 
% Exclude nominal variables (some were recoded in |sm_file_update.py| in the 
% HCP-CCA replication scripts):

exclVars = ["subject", "sex", "quarter/release", "race", "color_vision", "handedness"];
exclVarIdx = arrayfun(@(x) find(strcmpi(SMs, x)), exclVars);
grotBB(exclVarIdx, :) = nan;
grotBBd(exclVarIdx, :) = nan;
%% 
% 
%% 
% Load detailed explations of the variables:

HCP_DataDict_File = fullfile('data', 'HCP_S1200_DataDictionary_April_20_2018.csv');
datDict = readtable(HCP_DataDict_File);
datDict = datDict(:, ["columnHeader", "fullDisplayName", "category", "assessment", "description"]);
%% 
% Add missing variables to the dictionary (some can actually be found on the 
% HCP <https://wiki.humanconnectome.org/display/PublicData/HCP-YA+Data+Dictionary-+Updated+for+the+1200+Subject+Release 
% website> but were somehow not included in the CSV):

tmp = SMs(~ismember(SMs, datDict.columnHeader));
tmp = struct('columnHeader', tmp, 'fullDisplayName', tmp, 'category', [], 'assessment', [], 'description', tmp);
datDict = [datDict; struct2table(tmp)];
%% 
% Add row names:

% Remove duplication and set column header as row names
datDict(strcmp(datDict.description, 'Age group of Participant.'), :) = [];
datDict.Properties.RowNames = datDict.columnHeader;

% Remove the trailing parentheses in the full display name
% e.g, in 'Variable Short Penn Line Orientation: Total Number Correct (VSPLOT_TC)'
datDict.DisplayName = regexprep(datDict.fullDisplayName, ' \([^\s]+\)$', '');

% head(datDict(:, ["DisplayName", "description"]))
%% 
% 
%% 
% Plot the correlation between the first two CCA modes and behavioral measurements 
% (after deconfounding - results similar without deconfounding):

nSM_plot = 18;
redundant_thres = 0.8;
maxFontSize = 36;

% After deconfounding

smcorr = corrcoef(grotBBd');  % Not the true correlation, but can still remove some redundant measures
[c, cidx] = sort(grotBBd(:, 1), 'descend', 'MissingPlacement', 'last');
cidx = cidx(~isnan(c));

% Remove redundant measures
pltIdx = nan(nSM_plot, 2);
pltIdx(1, :) = [cidx(1) cidx(end)];
i1 = 2; i2 = length(cidx) - 1;
for i = 2:nSM_plot
    while max(smcorr(cidx(i1), pltIdx(1:i - 1, 1))) >= redundant_thres
        i1 = i1 + 1;
    end
    pltIdx(i, 1) = cidx(i1);
    i1 = i1 + 1;
    while max(smcorr(cidx(i2), pltIdx(1:i - 1, 2))) >= redundant_thres
        i2 = i2 - 1;
    end
    pltIdx(i, 2) = cidx(i2);
    i2 = i2 - 1;
end

% Colored by included or not
tmp_color = [0.5 0.5 0.5; 0 0 1];
included = zeros(numel(SMs), 1); included(varskeep) = 1;

for k = 1:2
    nexttile(9 * k - 8, [3 2])
    hold on;
    rg = [min(grotBBd(pltIdx(:, k), 1)) max(grotBBd(pltIdx(:, k), 1))];
    tmpy = linspace(rg(1), rg(2), nSM_plot);
    if k == 1
        tmpy = tmpy(end:-1:1);
    end
    rg = rg + [-0.05 0.05];
    for i = 1:nSM_plot
        text(1, tmpy(i), datDict{SMs{pltIdx(i, k)}, 'DisplayName'}, ...
            'Color', tmp_color(included(pltIdx(i, k)) + 1, :), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Interpreter', 'none', 'FontSize', abs(grotBBd(pltIdx(i, k), 1)) * maxFontSize);
    end
    ylim(rg);
    xlim([0.8 1.2]);
    text(0.8, tmpy(1), sprintf('\\rho = %.2f', tmpy(1)), "FontSize", abs(tmpy(1)) * 30, ...
        'Color', 'k', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    tmpI = round(length(tmpy) / 2); tmpP = grotBBd(pltIdx(tmpI, k), 1);
    text(0.8, tmpy(tmpI), sprintf('\\rho = %.2f', tmpP), "FontSize", abs(tmpP) * maxFontSize, ...
        'Color', 'k', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    text(0.8, tmpy(end), sprintf('\\rho = %.2f', tmpy(end)), "FontSize", abs(tmpy(end)) * 30, ...
        'Color', 'k', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    axis off
    if k == 1
        title(sprintf('Posthoc correlation with mode %d (de-confound)', 1), ...
            'FontSize', BaseFont + 2, 'FontWeight', 'bold');
    end
end

nexttile('east'); hold on; axis off 
text(0.5, 0.52, 'Included in CCA', 'Color', [0 0 1], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', BaseFont);
text(0.5, 0.48, 'Excluded', 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', BaseFont);
xlim([0.4 0.6]); ylim([0 1])


%%

title(t, 'Canonical Correlation Analysis between MINDy connectivity W and HCP phenotypic measures', ...
    'FontSize', 20, 'FontWeight', 'bold');
PrintAsSeen(fullfile('figures', 'HCP_Rest_FIX_Simple_Mdl200_sess', 'CCA_MINDy200_W'), '-dsvg', '-vector');
PrintAsSeen(fullfile('figures', 'HCP_Rest_FIX_Simple_Mdl200_sess', 'CCA_MINDy200_W'), '-dpdf', '-vector');