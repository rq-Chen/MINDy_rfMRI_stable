%% Flip the majority of sign to negative or positive
function dat = FlipSign(dat, dim, sgn)
    if nargin < 3
        sgn = -1;
    end
    if nargin < 2
        dim = 1;
    end
    tmp = sign(mean(sign(dat), dim));
    tmp(tmp == 0) = 1;
    dat = dat .* tmp * sgn;
end


