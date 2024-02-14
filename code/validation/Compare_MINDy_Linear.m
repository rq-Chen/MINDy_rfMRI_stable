%% Compare_MINDy_Linear - Compare the goodness of fit of MINDy and linear models

clear; clc; close all;
init_prj('MINDy_Base_v1.0');


%% Parameters

datFile = fullfile('data', 'allDat100.mat');
runIdx = {1:2, 3:4};
nSubs = 3;  % Use a subset of data
nSess = size(runIdx, 2);

% Training
derivType = 'two-point';
nBatch = 5000;  % Actual number of batches used for training
batchPerCP = 50;  % Save checkpoint after this number of batches
nStopCP = nBatch / batchPerCP;

% Whether to exclude the diagonal elements from reliability estimation or not
trr_no_diag = true;

load(datFile, 'allDat');


%% Train models

mindys = cell(nSubs, nSess);
linears = cell(nSubs, nSess);
denses = cell(nSubs, nSess);  % Dense linear model
densesE = cell(nSubs, nSess);  % Dense linear model error
for i = 1:nSubs
    tmpX = cell(1, nSess);
    tmpDX = cell(1, nSess);
    for j = 1:nSess
        [mindys{i, j}, Dat, dDat, tmpX{j}, tmpDX{j}] = MINDy_from_Deconv(allDat(i, runIdx{j}), 'Base', derivType);
        linears{i, j} = MINDy_from_Deconv(allDat(i, runIdx{j}), 'Linear', derivType);
        Dat = [Dat{:}];
        dDat = [dDat{:}];
        denses{i, j} = (Dat' \ dDat')';
        densesE{i, j} = mean((denses{i, j} * Dat - dDat) .^ 2, 2);
    end
    tmpX = cellfun(@(x) x{1}, tmpX, 'UniformOutput', false);
    tmpDX = cellfun(@(x) x{1}, tmpDX, 'UniformOutput', false);
    % Cross validation (pending)
end


%% Check reliability and error

% Number of recorded checkpoints
nCheck = size(mindys{1, 1}.RecW{1}, 3);

% Compute the Pearson correlation for W between sessions
diagTRR = zeros(nCheck, nSubs, 3);
TRR = zeros(nCheck, nSubs, 3);  % [MINDy, Linear, denses]
for i = 1:nSubs
    for j = 1:nCheck
        W1 = mindys{i, 1}.RecW{1}(:, :, j);
        W2 = mindys{i, 2}.RecW{1}(:, :, j);
        diagTRR(j, i, 1) = corr(W1(:), W2(:), "rows", "complete");
        if trr_no_diag
            for k = 1:size(W1, 1)
                W1(k, k) = nan;
                W2(k, k) = nan;
            end
        end
        TRR(j, i, 1) = corr(W1(:), W2(:), "rows", "complete");
        W1 = linears{i, 1}.RecW{1}(:, :, j);
        W2 = linears{i, 2}.RecW{1}(:, :, j);
        diagTRR(j, i, 2) = corr(W1(:), W2(:), "rows", "complete");
        if trr_no_diag
            for k = 1:size(W1, 1)
                W1(k, k) = nan;
                W2(k, k) = nan;
            end
        end
        TRR(j, i, 2) = corr(W1(:), W2(:), "rows", "complete");
    end
    W1 = denses{i, 1};
    W2 = denses{i, 2};
    diagTRR(:, i, 3) = corr(W1(:), W2(:), "rows", "complete");
    if trr_no_diag
        for k = 1:size(W1, 1)
            W1(k, k) = nan;
            W2(k, k) = nan;
        end
    end
    TRR(:, i, 3) = corr(W1(:), W2(:), "rows", "complete");
end

% Compute the magnitude of error
err = zeros(nCheck, nSubs, nSess, 3);  % [MINDy, Linear, denses]
for i = 1:nSubs
    for j = 1:nSess
        err(:, i, j, 1) = mean(mindys{i, j}.E);
        err(:, i, j, 2) = mean(linears{i, j}.E);
        err(:, i, j, 3) = mean(densesE{i, j});
    end
end
% Collapse across sessions
err = squeeze(mean(err, 3));


%% Visualize

figure;
tiledlayout(1, 3);
nexttile(1, [1 2]);

yyaxis left;
plot(squeeze(mean(err, 2)), 'LineWidth', 2);
ylabel('Mean square error (log scale)');
yticks([0.1:0.1:0.5 1 5 10]);
ylim([0.05 20])
set(gca, 'YScale', 'log')
hold on
tmpI = nStopCP:nStopCP:nCheck;
tmpC = gca().Children(1).Color;
for i = 1:length(tmpI)
    tmp = mean(err(tmpI(i), :, 1), "all");
    text(tmpI(i), tmp * 1.1, sprintf("(%d, %.2f)", tmpI(i), tmp), ...
        'HorizontalAlignment','center', 'VerticalAlignment', 'bottom', 'Color', tmpC)
end
tmp = mean(tmpDX{1} .^ 2, 'all');
yline(tmp, '--k');
text(nCheck / 2, tmp * 1.2, 'Mean squares of the target')

yyaxis right;
plot(squeeze(mean(TRR, 2)), 'LineWidth', 2);
if trr_no_diag
    ylabel('Mean \rho(W_{sess1}, W_{sess2}), off-diagonal only');
else
    ylabel('Mean \rho(W_{sess1}, W_{sess2})');
end
lgd = legend({'MINDy', 'Sparse Linear', 'Dense Linear'});
lgd.AutoUpdate = "off";
xlabel('Checkpoint');
hold on;
tmpI = nStopCP:nStopCP:nCheck;
tmpC = gca().Children(1).Color;
for i = 1:length(tmpI)
    tmp = mean(TRR(tmpI(i), :, 1), "all");
    text(tmpI(i), tmp - 0.02, sprintf("(%d, %.2f)", tmpI(i), tmp), ...
        'HorizontalAlignment','center', 'VerticalAlignment', 'top', 'Color', tmpC)
end
xline(nStopCP, '--k');
text(105, mean(ylim()), 'Default stopping point');

xlim(xlim() .* [1 1.1])
title(sprintf('Error and Reliability of MINDy, sparse linear and dense linear models (%d subjects)', nSubs));
subtitle(sprintf('%d parcel MINDy, %s derivative', size(mindys{1, 1}.RecW{1}, 1), derivType));


%% Check whether training fails (parameters become all zero)

Wnorm = cellfun(@(x) squeeze(pagenorm(x.RecW{1}, 'fro')), mindys, 'UniformOutput', false);  % (nCheck, 1) each cell
Wnorm = cell2mat(reshape(Wnorm, 1, []));

nexttile(3);
plot(mean(Wnorm, 2))
xlabel('Checkpoint');
ylabel('Mean Frobenius norm of W');
title('Mean Frobenius norm of MINDy W across subjects and sessions')


%% Test MINDy on data generated by the dense linear system

% Size
nTR = 2000;
nParcels = size(denses{1, 1}, 1);

% Model
tmpDat = [allDat{1, runIdx{1}}];
tmpFC = corr(tmpDat');
tmpX = tmpDat(:, 1:end - 1);
tmpY = tmpDat(:, 2:end);
A = (tmpX' \ tmpY')';
epstd = std(A * tmpX - tmpY, 0, 2);

% Generate data
tmpDat = nan(nParcels, nTR);
tmpDat(:, 1) = randn(nParcels, 1);
for i = 2:nTR
    tmpDat(:, i) = A * tmpDat(:, i - 1) + randn(nParcels, 1) .* epstd;
end

% Run MINDy
mdl = MINDy_from_Deconv(tmpDat, 'Base', derivType, nBatch);  % Default stopping point

% Get the fitted MINDy model at nBatch batches (checkpoint nStopCP)
mdl0 = struct();
mdl0.Param = mindys{1, 1}.Param;
mdl0.Param{2} = mindys{1, 1}.RecA{1}(:, nStopCP);
mdl0.Param{5} = mindys{1, 1}.RecW{1}(:, :, nStopCP);
mdl0.Param{6} = mindys{1, 1}.RecD{1}(:, nStopCP);


%% Visualization

figure('Units', 'inches', 'Position', [0 0.5 20 10.0208]);
t = tiledlayout(2, 3);

nexttile;
MyNetMatYeo(tmpFC);
cl = clim();
title('Data FC');

nexttile;
MyNetMatYeo(mdl0.Param{5});
cl2 = clim();
title('MINDy connectivity fitted on data')

nexttile;
[eq, lc, sims] = GetEqLC(mdl0);
lcmtf = LC2Motif(lc);
PlotAttractorLandscape(sims, eq, lcmtf, 3, 20);
title('MINDy landscape fitted on data');

nexttile;
MyNetMatYeo(corr(tmpDat'));
clim(cl);
title('Linear system simulation FC');

nexttile;
MyNetMatYeo(mdl.Param{5});
clim(cl2);
title('MINDy connectivity fitted on simulation');

nexttile;
[eq, lc, sims] = GetEqLC(mdl);
lcmtf = LC2Motif(lc);
PlotAttractorLandscape(sims, eq, lcmtf, 3, 20);
title('MINDy landscape fitted on simulation');

title(t, sprintf('%d parcel MINDy fitted on real and surrogate data', nParcels))
subtitle(t, sprintf('%s derivative', derivType))