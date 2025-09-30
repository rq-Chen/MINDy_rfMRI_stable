%% PLOTYEOSURFACE - Plot Schaefer parcellated data on brain surface
%
%  Importantly, make sure the Schaefer parcellation `dlabel` file under
%  `data/atlas/` folder is the one you used to parcellate your data.
%  The version provided here is not the latest one, but it is the one
%  we used to process our data in the `Singh2020Preproc` repository.
%  Inconsistency between parcellation version will lead to wrong parcel
%  ordering.
%
%  To use different brain geometry file, replace `nL` and `nR` with
%  desired filenames in `data/atlas/` folder.
function[Out,ommpL,ommpR,osL,osR]=PlotYeoSurface(X,mmpL,mmpR,sL,sR)

% fprintf("Remark: PlotYeoSurface is for pre-sorting data (original space).");

bDir=fullfile('data', 'atlas');

HZspace=5;
VZspace=5;
szX=size(X);
if szX(1)==1
    X=X';
end
if mod(numel(X),100)~=0
    doSub='y';
    SubX=X((100*floor(numel(X)/100))+1:end);
    SubMask=load(fullfile(bDir, 'SubMask.mat'));
    SubMask=SubMask.SubMask;
else
    doSub='n';
end

X=X(1:(100*floor(numel(X)/100)));
nY=numel(X);

if nargin == 1
    nL='S900.L.surf.gii';
    nR='S900.R.surf.gii';
    sL=gifti(fullfile(bDir,nL));
    sR=gifti(fullfile(bDir,nR)); 
    mmpL=gifti(fullfile(bDir,'mmpL.func.gii'));
    mmpR=gifti(fullfile(bDir,'mmpR.func.gii'));
    dlabel=ft_read_cifti(fullfile(bDir,['Schaefer2018_' num2str(nY) 'Parcels_17Networks_order.dlabel.nii']));
    yy=dlabel.parcels;
    mmpL.cdata=yy(1:end/2);
    mmpR.cdata=yy((1+end/2):end);
end
ommpL = mmpL;
ommpR = mmpR;
osL = sL;
osR = sR;

xL=mmpL;
xR=mmpR;

for i=setdiff(unique(mmpL.cdata),0)'
    xL.cdata(mmpL.cdata==i)=X(i);
end
for i=setdiff(unique(mmpR.cdata),0)'
    xR.cdata(mmpR.cdata==i)=X(i);
end

BadVertL=(mmpL.cdata==0);
BadVertR=(mmpR.cdata==0);


OrigIndL=1:numel(mmpL.cdata);
OrigIndR=1:numel(mmpR.cdata);


xL.cdata(BadVertL)=[];
xR.cdata(BadVertR)=[];

sL.vertices(BadVertL,:)=[];
sR.vertices(BadVertR,:)=[];

[fL,~]=find(ismember(sL.faces,find(BadVertL)));
[fR,~]=find(ismember(sR.faces,find(BadVertR)));
sL.faces(unique(fL),:)=[];
sR.faces(unique(fR),:)=[];

newL=1:numel(xL.cdata);
newR=1:numel(xR.cdata);

OrigIndL(~BadVertL)=newL;
OrigIndL(BadVertL)=nan;
OrigIndR(~BadVertR)=newR;
OrigIndR(BadVertR)=nan;
sL.faces=OrigIndL(sL.faces);
sR.faces=OrigIndR(sR.faces);

sR.vertices(:,1)=-(sR.vertices(:,1));
sL.vertices=sL.vertices-min(sL.vertices,[],1);
sR.vertices=sR.vertices-min(sR.vertices,[],1);

sR.vertices(:,2)=HZspace+max(sR.vertices(:,2))-sR.vertices(:,2)+max(sL.vertices(:,2));
sC=sL;
sC.vertices=[sL.vertices;sR.vertices];
sC.faces=[sL.faces;sR.faces+numel(xL.cdata)];
xC=xL;xC.cdata=[xC.cdata;xR.cdata];

sC.faces=[sC.faces;sC.faces+size(sC.vertices,1)];
tmp=(sC.vertices.*[-1 -1 1])+[0 0 VZspace+max(sC.vertices(:,3))];
tmp(:,2)=tmp(:,2)-min(tmp(:,2));

ppL=tmp(1:size(sL.vertices,1),2);
ppR=tmp((1+numel(ppL)):size(sC.vertices,1),2);
ppL2=(ppL-min(ppL))+min(ppR);
ppR2=(ppR-max(ppR))+max(ppL);
tmp(:,2)=[ppL2;ppR2];

sC.vertices=[sC.vertices;tmp];
xC.cdata=repmat(xC.cdata,2,1);

if strcmpi(doSub(1),'y')
    aa{1}=sC;aa{2}=xC;
    load(fullfile(bDir, 'SubMask.mat')); %#ok<LOAD>
    VZspace=0; %#ok<NASGU>
    HZspace=5;
    numCort=size(aa{1}.vertices,1);
    numR=size(SubMask.Right.vertices,1);
    numL=size(SubMask.Left.vertices,1);
    sR=SubMask.Right;
    sL=SubMask.Left;
    sR.vertices=sR.vertices(:,[2 1 3])*2;
    sL.vertices=sL.vertices(:,[2 1 3])*2;
    CortMax=max(aa{1}.vertices,[],1);
    CortMin=min(aa{1}.vertices,[],1);
    sR.vertices=sR.vertices-min(sR.vertices,[],1);
    sL.vertices=sL.vertices-min(sL.vertices,[],1);
    sLturn=sL;sLturn.vertices=sLturn.vertices.*[-1 -1 1];
    sRturn=sR;sRturn.vertices=sRturn.vertices.*[-1 -1 1];
    sRturn.vertices=sRturn.vertices-min(sRturn.vertices,[],1);
    sLturn.vertices=sLturn.vertices-min(sLturn.vertices,[],1);
    sRmax=max(sR.vertices,[],1); %#ok<NASGU>
    sLmax=max(sL.vertices,[],1); %#ok<NASGU>

    sRbottom=sRturn;sRbottom.vertices=[0 HZspace+CortMax(2) 0]+sRturn.vertices;
    sRbottom.vertices(:,3)=CortMax(3)+(sRbottom.vertices(:,3)-max(sRbottom.vertices(:,3)));
    sRbottom.vertices(:,1)=CortMin(1)+(sRbottom.vertices(:,1)-min(sRbottom.vertices(:,1)));

    sLtop=sLturn;sLtop.vertices=[0 HZspace+max(sRbottom.vertices(:,2)) 0]+sLturn.vertices;
    sLtop.vertices(:,3)=CortMin(3)+(sLtop.vertices(:,3)-min(sLtop.vertices(:,3)));
    sLtop.vertices(:,1)=CortMax(1)+(sLtop.vertices(:,1)-max(sLtop.vertices(:,1)));



    sLbottom=sL;sLbottom.vertices=[0 HZspace+max(sRbottom.vertices(:,2)) 0]+sL.vertices;
    sLbottom.vertices(:,3)=CortMax(3)+(sLbottom.vertices(:,3)-max(sLbottom.vertices(:,3)));
    sLbottom.vertices(:,1)=CortMin(1)+(sLbottom.vertices(:,1)-min(sLbottom.vertices(:,1)));

    sRtop=sR;sRtop.vertices=[0 HZspace+CortMax(2) 0]+sR.vertices;
    sRtop.vertices(:,3)=CortMin(3)+(sRtop.vertices(:,3)-min(sRtop.vertices(:,3)));
    sRtop.vertices(:,1)=CortMax(1)+(sRtop.vertices(:,1)-max(sRtop.vertices(:,1)));


    sLbottom.faces=sLbottom.faces+numCort;
    sRbottom.faces=numL+numCort+sRbottom.faces;
    sLtop.faces=numL+numR+numCort+sLtop.faces;
    sRtop.faces=(2*numL)+numR+numCort+sRtop.faces;
    SubCdata=repmat([sLbottom.facevertexcdata; sRbottom.facevertexcdata],2,1);
    unique(SubCdata);
    SubCdata=SubX(SubCdata);
    aa{2}.cdata=[aa{2}.cdata; SubCdata];
    aa{1}.vertices=[aa{1}.vertices;sLbottom.vertices;sRbottom.vertices;sLtop.vertices;sRtop.vertices];
    aa{1}.faces=[aa{1}.faces;sLbottom.faces;sRbottom.faces;sLtop.faces;sRtop.faces];
    sC=aa{1};xC=aa{2};
end


% figure;
plotGifti(gca, sC, xC);
if nargout>0
    Out={sC,xC};
else
    Out=[];
end
%figure;
%plot(sL,xL);
%hold on;
%plot(sR,xR);

view([90 0])

end