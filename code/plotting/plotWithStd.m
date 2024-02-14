%% Plot lines with shades representing SEM or Std on current figure
%
% Author: Ruiqi Chen
% Version: 12/14/2021
%
% Usage: ax = plotWithStd(x, y, std, varargin)
%
% `x` can be [], in which case the data will be plotted along the longest
% dimension of `y` with `x = 1:length(y)`. `varargin` is passed to `plot()`
% function to specify the line style.
%
% `std` should have the same size with `y`.
%
% The return value is the same as normal `plot`.
function ax = plotWithStd(x, y, sem, varargin)

if isempty(x)
    x = 1:length(y);
end

if iscolumn(x)
    x = x';
end
if size(y, 2) ~= size(x, 2)
    y = y';
end
if size(sem, 2) ~= size(x, 2)
    sem = sem';
end
ax = plot(x, y, varargin{:});

currHold = ishold;
if ~currHold
    hold on;
end
for i = 1:length(ax)
    fill([x, x(end:-1:1)], [y(i, :) - sem(i, :), ...
        y(i, end:-1:1) + sem(i, end:-1:1)], ...
        ax(i).Color, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
end
if ~currHold
    hold off;
end

end