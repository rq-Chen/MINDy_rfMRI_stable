function f = PrintAsSeen(varargin)
%PRINTASSEEN Print a figure as it appears on the screen
%
%   Syntax: (varargin will be passed to print)
%      f = PrintAsSeen(varargin);
%      f = PrintAsSeen(f, varargin);
%
%   The physical size of the output will be the same as
%   the figure's size on the screen and there will be no
%   white space around the figure.

if ishandle(varargin{1})
    f = varargin{1};
    varargin(1) = [];
else
    f = gcf;
end

% Get the physical size of the figure
set(f, 'Units', 'inches');
pos = get(f, 'Position');
width = pos(3);
height = pos(4);

% Set the paper size to match the figure's size
set(f, 'PaperUnits', 'inches', ...
    'PaperSize', [width height], ...
    'PaperPosition', [0 0 width height]);

% Print the figure
print(f, varargin{:});

end