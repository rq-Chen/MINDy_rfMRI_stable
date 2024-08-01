function[]=MyNetMatYeo(MAT,nNet,ShowXText,ShowYText,netLines)
%MYNETMATYEO Plots a nX by nX matrix in Yeo network order
%
%  Inputs:
%  - MAT: nX by nX matrix to plot
%  - nNet: Number of Yeo networks to plot (7 or 17, default 17)
%  - ShowXText: 'bottom', 'top' or 'none' (default 'bottom')
%  - ShowYText: 'left', 'right' or 'none' (default 'left')
%  - netLines: Thick lines separating networs, true or false (default true)

nX=size(MAT,1);
if nargin<2 || isempty(nNet)
    nNet=17;
end
if nargin<3 || isempty(ShowXText)
    ShowXText='bottom';
end
if nargin<4 || isempty(ShowYText)
    ShowYText='left';
end
if nargin<5 || isempty(netLines)
    netLines=true;
end
assert(nNet==7 || nNet==17,'nNet must be 7 or 17')
nSubCort=mod(nX,100);
nCort=nX-nSubCort;
gg=GetYeoNetworks(nCort,nNet);
MAT(:,1:nCort)=MAT(:,gg.SortInd);
MAT(1:nCort,:)=MAT(gg.SortInd,:);

matrix=MAT;

if nSubCort==0
    networks = gg.Names;
else
    networks=[gg.Names {'SubCort'}];
end

colors_new = distinguishable_colors(numel(networks));
tickpos=zeros(1,numel(networks));
for i=1:numel(gg.Names)
    tickpos(i)=median(find(gg.SortNet==i));
    atlas_params.transitions(i)=find(gg.SortNet==i,1,'last');
end
if nSubCort~=0
    tickpos(numel(networks))=nCort+(nSubCort)/2;
else
    atlas_params.transitions(end)=[];
end

imagesc(matrix);
colormap(gca,'turbo');

ax = axis;

set(gca,'XTick',[],'YTick',[],'Xlim',[ax(1) ax(2)]);

set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

if ~strcmpi(ShowXText,'none')
    if strcmpi(ShowXText,'bottom')
        tx=text(tickpos,ones(1,length(tickpos))*(nX+1),networks);
        set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',30);
    else
        tx=text(tickpos,ones(1,length(tickpos))*(-1),networks);
        set(tx,'HorizontalAlignment','right','VerticalAlignment','bottom','Rotation',-30);
    end
    for i = 1:length(tx)
        set(tx(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');  
    end
end

if ~strcmpi(ShowYText,'none')
    if strcmpi(ShowYText,'left') 
        ty=text(ones(1,length(tickpos))*(-1),tickpos,networks);
        set(ty,'HorizontalAlignment','right','VerticalAlignment','middle')
    else
        ty=text(ones(1,length(tickpos))*(nX+2),tickpos,networks);
        set(ty,'HorizontalAlignment','left','VerticalAlignment','middle')
    end
    for i = 1:length(ty)
        set(ty(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
    end
end

colorbar;
set(gca,'FontWeight','bold','FontSize',10,'LineWidth',2.5);


%Specify color and width of horizontal and verticel axis lines that
%seperate networks
if netLines
    for p = 1:numel(atlas_params.transitions)
        x = line([atlas_params.transitions(p) atlas_params.transitions(p)]+0.5,[0 nX]+0.5);
        set(x,'Color', 'Black','LineWidth',1.5);
        y = line([0 nX]+0.5,[atlas_params.transitions(p) atlas_params.transitions(p)]+0.5);
        set(y,'Color', 'Black','LineWidth',1.5);
    end
end
