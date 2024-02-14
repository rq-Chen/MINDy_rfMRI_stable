function lgd = PlotAttractorLandscape(sims, eq, LCmotifs, dim, nSkipTR, varargin)
%PLOTATTRACTORLANDSCAPE Plot the simulated trajectories and identified motifs
%
%   Plot the trajectories in PCA space along with identified stable fixed
%   points and limit cycles PC1 & PC2.
%
%   Inputs:
%     sims: (nParcels, T, nSims) array of simulated trajectories
%     eq: optional, (nParcels, nEq) matrix of stable fixed points.
%     LCmotifs: optional, (1, nLC) cell. Output of |LC2Motif()|.
%     dim: visualization dimension, either 2 or 3 (default 3).
%     nSkipTR: number of initial TRs to skip in the plot (default 25). Skipping
%       initial TRs makes the plot less cluttered.
%     varargin: additional arguments for |myPCA|, e.g., projecting to specific
%       PCs.
%
%   Outputs:
%     lgd: legend handle. The legend is not displayed by default. Set lgd.Visible
%       to 'on' to display it.
%
%   Note that the default plots might contain too many data points for
%   vector format. Try to use |sims(:, 1:5:end, 1:60)| to down sample the
%   data if the figure needs to be saved in vector format.


%% Handle inputs

if nargin < 5 || isempty(nSkipTR)
    nSkipTR = 25;
end
if nargin < 4 || isempty(dim)
    dim = 3;
end
assert(dim == 2 || dim == 3, "dim must be 2 or 3.");
assert(size(sims, 2) > nSkipTR, "|nSkipTR| must be smaller than the number of TRs in |sim|.");
if nargin < 3 || isempty(LCmotifs)  % (1, nLC) of (nParcels, 5) matrices
    LCmotifs = {};
end
if nargin < 2 || isempty(eq)  % (nParcels, nEq) matrix
    eq = [];
end


%% Plotting parameters (all using |plot| or |plot3|)

LCGhostOnly = true;  % If true, only plot LCmotif 1 & 2 but no lines

% Scale the sizes by the size of current axis
ax = gca;
un = get(ax, 'Units');
set(ax, 'Units', 'points');
pos = get(ax, 'Position');
scaling = pos(3) / 259.6;  % Parameters were tuned for a 259.6 * 244.0 axis
set(ax, 'Units', un);
clear ax un pos

% For (trajectory, start, end, fixed point, LCax1, LCax2) respectively:
% Color scheme
CLS = [70, 70, 70; ...
    0, 0, 0; ...
    204, 187, 68; ...
    238, 102, 119; ...
    68, 119, 170; ...
    34, 136, 51] / 255;  % https://personal.sron.nl/~pault/#sec:qualitative
CLS = mat2cell(CLS, ones(size(CLS, 1), 1));
% Line style
if LCGhostOnly
    LSTY = {'-'; 'none'; 'none'; 'none'; 'none'; 'none'};
else
    LSTY = {'-'; 'none'; 'none'; 'none'; '-.'; ':'};
end
% Line width
LW = [0.3; 0.5; 0.5; 0.5; 3; 3] * scaling;
LW = mat2cell(LW, ones(size(LW, 1), 1));
% Marker shape
MK = {'none'; 'o'; 'o'; 'pentagram'; '^'; 'v'};
% Marker size
MKSZ = [2; 2; 4; 14; 8; 8] * scaling;
MKSZ = mat2cell(MKSZ, ones(size(MKSZ, 1), 1));
% Legend strings
if LCGhostOnly
    ldgstr = {'Trajectories', 'Start', 'End', 'Stable fixed point', ...
        'Limit cycle slowest point'};
else
    ldgstr = {'Trajectories', 'Start', 'End', 'Stable fixed point', ...
        'Limit cycle major axis', 'Limit cycle minor axis'};
end


%% Compute PC projections

[score, coeff, ~, mu] = myPCA(sims, dim, varargin{:});  % score: (T, nSims, dim); coeff: (nParcels, dim); mu: (1, nParcels)
if ~isempty(eq)
    eq = (eq' - mu) * coeff;  % (nEq, dim)
end
if ~isempty(LCmotifs)
    LCproj = cellfun(@(x) (x' - mu) * coeff, LCmotifs, 'UniformOutput', false);  % (1, nLC) of (5, dim) matrices
    LCproj = cat(3, LCproj{:});  % (4 or 5, dim, nLC)
    LCproj = permute(LCproj, [1 3 2]);  % (4 or 5, nLC, dim)
    if LCGhostOnly && size(LCproj, 2) > 1  % More than one LC (thus LC might be asymmetric)
        LCax1 = LCproj(1, :, :);  % Only the slowest point
    else
        LCax1 = LCproj([3 1], :, :);  % neg major - pos major
    end
    LCax2 = LCproj([4 2], :, :);  % neg minor - pos minor
else
    LCax1 = [];
    LCax2 = [];
end


%% Plot

hold on
plt_dat = {score(nSkipTR + 1:end, :, :); score(nSkipTR + 1, :, :); score(end, :, :); eq; LCax1; LCax2};
if LCGhostOnly
    plt_res = cellfun(...
        @(x, cls, lsty, lw, mk, mksz) myplot(x, dim, 'LineStyle', lsty, ...
            'LineWidth', lw, 'Marker', mk, 'MarkerSize', mksz, 'Color', cls, ...
            'MarkerFaceColor', cls, 'MarkerEdgeColor', cls), ...
        plt_dat(1:5), CLS(1:5), LSTY(1:5), LW(1:5), MK(1:5), MKSZ(1:5), 'UniformOutput', false);    
else
    plt_res = cellfun(...
        @(x, cls, lsty, lw, mk, mksz) myplot(x, dim, 'LineStyle', lsty, ...
            'LineWidth', lw, 'Marker', mk, 'MarkerSize', mksz, 'Color', cls, ...
            'MarkerFaceColor', cls, 'MarkerEdgeColor', cls), ...
        plt_dat, CLS, LSTY, LW, MK, MKSZ, 'UniformOutput', false);
end

% View
if dim == 3
    view([-37.5, 30])
end
grid on

% Axes
axis tight
xlabel("PC1")
ylabel("PC2")
if dim == 3
    zlabel("PC3")
end

% Legend
lgd = legend(cellfun(@(x) x(1), plt_res), ldgstr);
set(lgd, 'Location', 'best', 'Visible', 'off');

end


%% Utilities
function lines = myplot(dat, dim, varargin)
%MYPLOT Plot with potentially empty input
%
%   Usage:
%     lines = myplot(dat, dim, varargin); % Plot first |dim| slices
%     lines = myplot([], dim, varargin);  % Empty plot
%     lines = myplot(dat, [], varargin);  % Plot all slices
if ndims(dat) == 2
    dat = reshape(dat, size(dat, 1), 1, []);
end
if nargin < 2 || isempty(dim)
    dim = size(dat, 3);
end
assert(dim == 2 || dim == 3, "dim must be 2 or 3.");
if isempty(dat)
    if dim == 2
        lines = plot(nan, nan, varargin{:});
    elseif dim == 3
        lines = plot3(nan, nan, nan, varargin{:});
    end
else
    if dim == 2
        lines = plot(dat(:, :, 1), dat(:, :, 2), varargin{:});
    elseif dim == 3
        lines = plot3(dat(:, :, 1), dat(:, :, 2), dat(:, :, 3), varargin{:});
    end
end

end


%% A similar one for scatter plots
function sc = myscatter(dat, dim, varargin)
%MYSCATTER Scatter plot with potentially empty input
%
%   Usage:
%     sc = myscatter(dat, dim, varargin); % Plot first |dim| slices
%     sc = myscatter([], dim, varargin);  % Empty plot
%     sc = myscatter(dat, [], varargin);  % Plot all slices

if ndims(dat) == 2
    dat = reshape(dat, size(dat, 1), 1, []);
end
if nargin < 2 || isempty(dim)
    dim = size(dat, 3);
end
assert(dim == 2 || dim == 3, "dim must be 2 or 3.");
if isempty(dat)
    if dim == 2
        sc = scatter(nan, nan, varargin{:});
    elseif dim == 3
        sc = scatter3(nan, nan, nan, varargin{:});
    end
else
    if dim == 2
        sc = scatter(dat(:, :, 1), dat(:, :, 2), varargin{:});
    elseif dim == 3
        sc = scatter3(dat(:, :, 1), dat(:, :, 2), dat(:, :, 3), varargin{:});
    end
end

end