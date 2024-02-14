%% Validation_Figure

clear; clc; close all;
init_prj('MINDy_Base_v1.0');


%% Parameters

% Model and data
mdlName = 'HCP_Rest_FIX_Simple_Mdl200_sess';
mdlFile = fullfile('data', [mdlName, '.mat']);
datFile = replace(mdlFile, 'Mdl', 'Deconv');

% Figure outputs
figdir = fullfile('figures', mdlName);
if ~exist(figdir, 'dir')
    mkdir(figdir);
end


%% Load data and model

load(mdlFile, 'allMdl', 'runIdx');
load(datFile, 'allDat');
[nSubs, nSess] = size(allMdl);
nParcels = size(allMdl{1, 1}.Param{5}, 1);

% Function to get data and derivatives
switch derivType
    case 'two-point'
        GetXDX = @(dat) deal(cell2mat(cellfun(@(x) x(:, 1:end-2), dat, "UniformOutput", false)), ...
            cell2mat(cellfun(@(x) (x(:, 3:end) - x(:, 1:end-2)) / 2, dat, "UniformOutput", false)));
    % other cases
end


%% Cross validation between sessions

cvR2 = nan(nSubs, nSess, nParcels, 3);  % (training, testing on other session, testing on other subject)
cvMSE = nan(nSubs, nSess, nParcels, 3);  % (training, testing on other session, testing on other subject)
for i = 1:nSubs
    disp("Subject " + i + "/" + nSubs);
    [X, DX] = cellfun(@(idx) GetXDX(allDat(i, idx)), runIdx(i, :), "UniformOutput", false);
    R2 = nan(nSess, nSess, nParcels);
    MSE = nan(nSess, nSess, nParcels);
    for j = 1:nSess
        for k = 1:nSess
            sims = MINDyInt_00(allMdl{i, j}, X{k}, 1, 1, 0, 1);
            sims = squeeze(sims(:, 2, :) - sims(:, 1, :));  % (nParcels, nTRs)
            MSE(j, k, :) = mean((sims - DX{k}).^2, 2);
            R2(j, k, :) = 1 - MSE(j, k, :) ./ reshape(var(DX{k}, 0, 2), 1, 1, []);
        end
    end
    for l = 1:nParcels
        cvR2(i, :, l, 1) = diag(R2(:, :, l));
        cvR2(i, :, l, 2) = offdiagnanmean(R2(:, :, l)');
        cvMSE(i, :, l, 1) = diag(MSE(:, :, l));
        cvMSE(i, :, l, 2) = offdiagnanmean(MSE(:, :, l)');
    end

    iSub = randi(nSubs);
    [X, DX] = cellfun(@(idx) GetXDX(allDat(iSub, idx)), runIdx(iSub, :), "UniformOutput", false);
    for j = 1:nSess
        sims = MINDyInt_00(allMdl{i, j}, X{j}, 1, 1, 0, 1);
        sims = squeeze(sims(:, 2, :) - sims(:, 1, :));  % (nParcels, nTRs)
        cvMSE(i, j, :, 3) = mean((sims - DX{j}).^2, 2);
        cvR2(i, j, :, 3) = 1 - cvMSE(i, j, :, 3) ./ reshape(var(DX{k}, 0, 2), 1, 1, []);
    end
end


%% Visualize cross validation results

tmpStr = {'training', 'testing on the other session', 'testing on another participant'};

f = figure('Units', 'inches', 'Position', [1 3 8 6]);
t = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');

nexttile();
boxchart(reshape(mean(cvR2, 3), [], 3));
xticklabels(tmpStr);
ylabel('Mean R^2');
title('Mean R^2 over parcels for all models')

cl = [min(mean(cvR2, [1 2]), [], "all") max(mean(cvR2, [1 2]), [], "all")];
for i = 1:3
    nexttile();
    PlotYeoSurface(squeeze(mean(cvR2(:, :, :, i), [1 2])));
    title(['Mean R^2 for ' tmpStr{i}]);
    clim(cl);
    colorbar;
end

set(findall(gcf,'type','text'), 'FontSize', 14, 'FontName','helvetica');
title(t, 'Training and cross-validated R squared', 'FontSize', 16);
PrintAsSeen(f, fullfile(figdir, ['cross_validation_R2_' mdlName]), '-dpng', '-r300');


%% Same visualization for MSE

f = figure('Units', 'inches', 'Position', [1 3 8 6]);
t = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');

nexttile();
boxchart(reshape(mean(cvMSE, 3), [], 3));
xticklabels(tmpStr);
ylabel('Mean MSE');
title('Mean MSE over parcels for all models')

cl = [min(mean(cvMSE, [1 2]), [], "all") max(mean(cvMSE, [1 2]), [], "all")];
for i = 1:3
    nexttile();
    PlotYeoSurface(squeeze(mean(cvMSE(:, :, :, i), [1 2])));
    title(['Mean MSE for ' tmpStr{i}]);
    clim(cl);
    colorbar;
end

set(findall(gcf,'type','text'), 'FontSize', 14, 'FontName','helvetica');
title(t, 'Training and cross-validated MSE', 'FontSize', 16);
PrintAsSeen(f, fullfile(figdir, ['cross_validation_MSE_' mdlName]), '-dpng', '-r300');


%% Utilities
function m = offdiagnanmean(M, varargin)
%OFFDIAGNANMEAN Mean of off-diagonal, non-NaN elements (if any; otherwise nan)
for i = 1:min(size(M))
    M(i, i) = nan;
end
m = mean(M, varargin{:}, "omitnan");
end