function varargout = plotGifti(varargin)
% Modified by Ruiqi Chen, 2023/02/07
% 
% plot method for GIfTI objects
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: plot.m 5888 2014-02-19 19:54:12Z guillaume $

% if ishandle(varargin{1})
%     h = figure(varargin{1});
% else
%     h = figure;
%     %axis equal;
%     %axis off;
%     %camlight;
%     %camlight(-80,-10);
%     %lighting phong;
% end
% cameramenu;

if nargin == 1 || ~ishandle(varargin{1})
    h = gcf;
    ax = axes('Parent', h);
    these = varargin;
else
    ax = varargin{1};
    % In order to handle layout organizer
    h = get(ax,'parent');
    while ~strcmp(get(h,'type'),'figure')
        h = get(h, 'parent');
    end
    these = varargin(2:end);
end
this = these{1};
if length(these) > 1
    cdata = subsref(these{2},struct('type','.','subs','cdata'));
else
    cdata = [];
end
if length(these) > 2
    indc = these{3};
else
    indc = 1;
end

axis(ax,'equal');
axis(ax,'off');
hp = patch(struct(...
    'vertices',  subsref(this,struct('type','.','subs','vertices')),...
    'faces',     subsref(this,struct('type','.','subs','faces'))),...
    'FaceColor', 'b',...
    'EdgeColor', 'none',...
    'Parent',ax);

if ~isempty(cdata)
    set(hp,'FaceVertexCData',cdata(:,indc), 'FaceColor','interp')
end

axes(ax);
camlight;
camlight(-80,-10);
lighting phong;
axes(ax);
% cameramenu;

if nargout
    varargout{1} = hp;
end
