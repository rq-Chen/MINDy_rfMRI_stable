%% Does the recording of variables in a separate module to clean code
%% Doesn't record error--need to do that separate
if strcmpi(ParStr.RecCorr(1),'y')
    xPt=(W0*(PolyZ)-((Decay).*rX0)+C);
    yPt=dY(:,bInd);
    Rec.RecCorr{iRep}(:,ceil(iBatch/RecWrate))=DiagCorr(yPt',xPt');
    Rec.RecTotCorr{iRep}(:,ceil(iBatch/RecWrate))=corr(xPt(:),yPt(:));
end
if strcmpi(doDisp(1),'y')
    disp([iBatch iRep])
end
if strcmpi(RecW(1),'y')
    %% Records W0 (the full connectivity matrix) instead of W1 (the sparse part)!
    Rec.RecW{iRep}(:,:,ceil(iBatch/RecWrate))=W0;
end
% Change to record the actual curvature alpha in equations
if iBatch==1&&(size(A,2)~=1)&&(ismatrix(ndims(Rec.RecA{iRep})))
    ss=size(Rec.RecA{iRep});
    Rec.RecA{iRep}=zeros(ss(1),size(A,2),ss(2));
end
if strcmpi(RecA(1),'y')
    if size(A,2)==1
        Rec.RecA{iRep}(:,ceil(iBatch/RecWrate))=sqrt(A./(B1(:,1).^2)-.25);
    else            
        Rec.RecA{iRep}(:,:,ceil(iBatch/RecWrate))=sqrt(A./(B1(:,1).^2)-.25);
    end      
end                                             
    
%% I changed to do B1 instead!!!
if strcmpi(RecB2(1),'y')
    Rec.RecB2{iRep}(:,ceil(iBatch/RecWrate))=B1(:,1);
end
if strcmpi(RecDdiff(1),'y')
    Rec.RecD{iRep}(:,ceil(iBatch/RecWrate))=Decay;  % Change this from D to Decay!
end