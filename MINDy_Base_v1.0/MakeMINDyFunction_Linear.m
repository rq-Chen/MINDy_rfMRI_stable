function[ooP]=MakeMINDyFunction_Linear(ooP)

% 'Tran': activation given input
% 'dTran': derivative of activation over input
% 'FastFun': time derivative of the dynamics given state
if ~isempty(ooP.Param{5})
    W=ooP.Param{5};
else
    W=ooP.Param{1};
end
slope=ooP.Param{2};
D=ooP.Param{6};
Wfast = (W .* slope - diag(D));

ooP.Tran=@(xx) slope .* xx;
ooP.dTran=@(xx) slope;
if (~isempty(ooP.Param{4}))&&(mean((ooP.Param{4})==0)~=1)
    c=ooP.Param{4};
    ooP.FastFun=@(yy)Wfast*yy+c;
else
    ooP.FastFun=@(yy)Wfast*yy;
end
end