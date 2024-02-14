function[]=MyNetMatYeo(MAT,doSort,makeShort,varargin)
%% Do sort indicates whether to resort into network blocks ('n'/'y')
%% Makeshort indicates whether to use the 17 networks ('n') or the 8 network ('y') divisions.
nX=size(MAT,1);
if nargin==2
    makeShort='n';
end
if nargin==1
%     disp('Assuming do Sort')
    doSort='y';
    makeShort='n';
end
nSubCort=mod(nX,100);
nCort=nX-nSubCort;
gg=GetYeoNetworks(nCort,makeShort);
if strcmpi(doSort(1),'y')
    MAT(:,1:nCort)=MAT(:,gg.SortInd);
    MAT(1:nCort,:)=MAT(gg.SortInd,:);
end

matrix=MAT;

if nSubCort==0
networks = gg.Names;
else
    networks=[gg.Names {'SubCort'}];
end

colors_new = distinguishable_colors(numel(networks));
tickpos=zeros(1,numel(networks));
for i=1:numel(gg.Names)
    tickpos(i)=median(find(ceil(gg.SortNet/2)==i));
    atlas_params.transitions(i)=find(ceil(gg.SortNet/2)==i,1,'last');
end
if nSubCort~=0
    tickpos(numel(networks))=nCort+(nSubCort)/2;
else
    atlas_params.transitions(end)=[];
end


%h = figure('Color',[0.9 0.9 0.9],'Position',[56 143 1295 807]); %[56 143 1095 807]

if nargin>7
   if ~isempty(varargin{1}) && ~isempty(varargin{2})
        climlow = varargin{1};
        climhigh = varargin{2};
       imagesc(matrix,[climlow climhigh]);
   else
       imagesc(matrix);
   end
else
  imagesc(matrix);
end
%load /data/cn5/caterina/PDgrant/scripts/better_jet_colormap.mat;
%colormap(better_jet_colormap_diff);    

colormap(gca,'jet');

%vline_new(atlas_params.transitions,'k',3);
%hline_new(atlas_params.transitions,'k',3);

ax = axis;

set(gca,'XTick',[],'YTick',[],'Xlim',[ax(1) ax(2)]);

set(gca,'XTicklabel','');
set(gca,'YTicklabel','');    

tx= text(tickpos,ones(1,length(tickpos))*(nX+1),networks);
set(tx,'HorizontalAlignment','right','VerticalAlignment','top','Rotation',45);

%specify color, font and size of network labels (lines 42 - 45)
for i = 1:length(tx)
     set(tx(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');  
end
set(gca,'FontWeight','bold','FontSize',10);


ty= text(-1*ones(1,length(tickpos)),tickpos,networks);
set(ty,'HorizontalAlignment','right','VerticalAlignment','middle')

for i = 1:length(ty)
    set(ty(i),'Color',[colors_new(i,1) colors_new(i,2) colors_new(i,3)],'FontName','Helvetica','FontSize',10,'FontWeight','bold');   
end
colorbar;
set(gca,'FontWeight','bold','FontSize',10, 'LineWidth'  ,2.5);


%Specify color and width of horizontal and verticel axis lines that
%seperate networks (lines 67 - 73)
for p = 1:numel(atlas_params.transitions)
       x = line([atlas_params.transitions(p) atlas_params.transitions(p)]+0.5,[0 nX]+0.5);
       set(x,'Color', 'Black','LineWidth',1.5);
       
       y = line([0 nX]+0.5,[atlas_params.transitions(p) atlas_params.transitions(p)]+0.5);
       set(y,'Color', 'Black','LineWidth',1.5);
end
