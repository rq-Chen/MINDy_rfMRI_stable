classdef EqFinder
    %EQFINDER Find equilirbiums of MINDy model
    %   Find the equilibriums of a MINDy model using the method from
    %   (Sussillo & Barak, 2013), searching for the local minimums of 
    %   $q(X) = 1/2 * |F(X)|^2$. Also compute the derivative and Jacobian
    %   of the model.
    
    properties
        W
        D
        b
        aob2
        p5ob
        nParcels
    end
    
    methods
        function obj = EqFinder(mdl)
            %EQFINDER Constructor
            %   Set the properties
            obj.W = mdl.Param{5};  % (nParcels, nParcels)
            obj.D = mdl.Param{6};  % (nParcels, 1)
            obj.b = mdl.Param{3}(1);  % double
            A = mdl.Param{2};  % (nParcels, 1)
            obj.aob2 = (A ./ obj.b) .^ 2;
            obj.p5ob = 0.5 ./ obj.b;
            obj.nParcels = size(obj.W, 1);
        end
        
        function out = dxdt(obj, X)
            %DXDT Get the derivative of the model at X
            %   Output size is the same as X (nParcels, ...)
            %   (X will be transposed if it's row vector)
            if isrow(X); X = X'; end          
            xp5b = X + obj.p5ob;
            xm5b = X - obj.p5ob;
            psix = obj.b .* (sqrt(obj.aob2 + xp5b .^ 2) - sqrt(obj.aob2 + xm5b .^ 2));
            out = pagemtimes(obj.W, psix) - obj.D .* X;
        end

        function out = J(obj, X)
            %J Get Jacobian of the model at X
            %   If X is of size (nParcels, ...) the output will be of size
            %   (nParcels, nParcels, ...).
            if isrow(X); X = X'; end          
            xp5b = X + obj.p5ob;
            xm5b = X - obj.p5ob;
            ppsipx = obj.b .* ((xp5b ./ sqrt(obj.aob2 + xp5b .^ 2)) - ...
                (xm5b ./ sqrt(obj.aob2 + xm5b .^ 2)));  % (nParcels, 1)
            out = obj.W .* reshape(ppsipx, [1 size(ppsipx)]) - diag(obj.D);
        end

        function [q, G, H] = f(obj, X)
            %F Get the value, gradient and Hessian approximation of q(X)
            if isrow(X); X = X'; end
            assert(isvector(X), 'Input must be a vector, not matrix!');
            xdot = obj.dxdt(X);
            q = 0.5 * sum(xdot .* xdot);
            if nargout > 1
                pfpx = obj.J(X);
                G = pfpx' * xdot;
            end
            if nargout > 2
                H = pfpx' * pfpx;
            end
        end
    end
end

