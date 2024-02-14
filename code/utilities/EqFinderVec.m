classdef EqFinderVec
    %EQFINDERVEC linear combination of EqFinders
    properties
        N  % Number of EqFinders
        M  % (N, 1) Array of EqFinders
        W  % (N, 1) Array of weights
        nParcels  % Number of parcels in the model
    end

    methods

        function obj = EqFinderVec(mdls, weights)
            if iscell(mdls)
                mdls = [mdls{:}];
            end
            if isa(mdls, 'struct')
                mdls = arrayfun(@(x) EqFinder(x), mdls);
            end
            assert(numel(mdls) == numel(weights))
            assert(isnumeric(weights), 'Weights must be numerical array!');
            obj.N = numel(mdls);
            obj.M = mdls(:);
            obj.W = weights(:);
            nParcels = arrayfun(@(x) x.nParcels, obj.M);
            assert(all(nParcels == nParcels(1)), 'All models must have the same number of parcels!');
            obj.nParcels = nParcels(1);
        end

        function out = dxdt(obj, X)
            out = arrayfun(@(m, w) w * m.dxdt(X), obj.M, obj.W, 'UniformOutput', false);
            tmp = ndims(out{1}) + 1;
            out = sum(cat(tmp, out{:}), tmp);
        end

        function out = J(obj, X)
            out = arrayfun(@(m, w) w * m.J(X), obj.M, obj.W, 'UniformOutput', false);
            tmp = ndims(out{1}) + 1;
            out = sum(cat(tmp, out{:}), tmp);
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

        function out = run(obj, X0, T, dt)
            %RUN Simulate the model
            %   Inputs:
            %       - X0: (nParcels, nSims), initial condition
            %       - T: integer, simulation length in units TR
            %       - dt: in units TR (default 1), must be 1/N (N is
            %       integer)
            %
            %   Output:
            %       - out: (nParcels, T + 1, nSims), value at 0~T TR
            %
            %   Note that dt only controls the simulation timestep but the
            %   output is always the state at integer TRs.
            if nargin < 4
                dt = 1;
            end
            Samp = 1 / dt;
            assert(T == round(T) && Samp == round(Samp));
            assert(size(X0, 1) == obj.nParcels);
            nSims = size(X0, 2);
            X0 = reshape(X0, obj.nParcels, 1, nSims);
            out = nan(obj.nParcels, T + 1, nSims);
            out(:, 1, :) = X0;
            for t = 1:T
                for i = 1:Samp
                    X0 = X0 + dt * obj.dxdt(X0);
                end
                out(:, t + 1, :) = X0;
            end
        end

    end
end