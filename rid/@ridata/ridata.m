classdef ridata
    properties
        ny = 0; nu = 0; nx = 0;
        M = 0; D = 0;
        T  = []; U = []; Y = [];
        N  = []; 
        dT = []; X = []; V = [];
        trialidx = [];
    end
    
    %% TODO
    % Determine how to select a fixed number of data-points.
    methods

        function d = ridata(Y,U,T,X,V,dT)
            if nargin  == 0; return; end
            if nargin < 3
                error('Y,U,T must be specified.');
            end
            
            if nargin < 6 && nargin > 3
                error('X,V,dT must be specified.');
            end
            
            
            if ~all(size(Y,2) == [size(U,2) size(T,2)]) ...
                    ~all(1 == [size(Y,1) size(U,1) size(T,1)])
                error('Y,U,T must all be 1xM cell arrays.');
            end
            
            d.M = prod(size(Y));
            if d.M < 1, error('No trials given.'); end
            
            d.ny = size(Y{1},1);
            d.nu = size(U{1},1);
            d.N  = zeros(1,d.M);
            d.trialidx = zeros(2,d.M);
            
            for i = 1:d.M,
                d.N(i) = length(T{i});
                d.trialidx(1,i) = sum([1 d.N(1:i-1)]);
                d.trialidx(2,i) = d.trialidx(1,i) + d.N(i) - 1;
                if size(T{i},1) ~= 1
                    error(['T{' int2str(i) '} must be 1xNi']);
                end
                if ~all(size(Y{i}) == [d.ny d.N(i)]);
                    error(['Y{' int2str(i) '} is incorrect size.']);
                end
                if ~all(size(U{i}) == [d.nu d.N(i)]);
                    error(['U{' int2str(i) '} is incorrect size.']);
                end
            end
            
            d.D = sum(d.N);
            d.T = [T{:}]; d.Y = [Y{:}]; d.U = [U{:}];
            
            if nargin > 3
                d.nx = size(X{1},1);
                for i = 1:d.M,
                    if ~all(size(X{i}) == [d.nx d.N(i)]);
                        error(['X{' int2str(i) '} is incorrect size.']);
                    end
                    if ~all(size(V{i}) == [d.nx d.N(i)]);
                        error(['V{' int2str(i) '} is incorrect size.']);
                    end
                end
                d.X = [X{:}]; d.V = [V{:}];
                if ~all(size(dT) == [1,1]), error('dT must be scalar.'); end
                d.dT = dT;
            else
                d.nx = 0;
                d.X = []; d.V = []; d.dT = [];
            end
        end
        
        function dmerged = merge(d1,d2)
            if d1.ny ~= d2.ny || d1.nu ~= d2.nu
                error('d1 and d2 must have same dimensions.');
            end
            if d1.dT ~= d2.dT
                error('d1 and d2 must have same dT.');
            end
            
            dmerged = ridata;
            dmerged.nx = d1.nx;
            dmerged.ny = d1.ny;
            dmerged.nu = d1.nu;
            dmerged.dT = d1.dT;
            dmerged.N = [d1.N d2.N];
            dmerged.D = sum(dmerged.N);
            dmerged.M = length(dmerged.N);
            dmerged.Y = [d1.Y d2.Y];
            dmerged.U = [d1.U d2.U];
            dmerged.T = [d1.T d2.T];
            dmerged.X = [d1.X d2.X];
            dmerged.V = [d1.V d2.V];
            dmerged.trialidx = [d1.trialidx (d1.D + d2.trialidx)];
        end
        
        function dsplit = split(d,n,idx)
        % dsplit = split(d,trialidx)
        % d        -- ridata structure.
        % n        -- trial number
        % idx      -- indicies, increasing
        %
        % Split trial n into length(idx)+1 separate trials.
        % Requires 1 < idx(i) < d.N(n)
        % The first new trial starts on 1 and ends at idx(i)-1.
        % The remaining trials start at idx(i) and end at
        %     idx(i+1)-1.
            
            
            if n < 1 || n > d.M || n ~= floor(n),
                error('n must be an integer in {1,...,d.M}'); 
            end
            
            if any(diff(idx) <= 0) || ...
                    any(idx ~= floor(idx)) || ...
                    any(idx == 1) || any(idx > d.N(n))
                error('idx must be increasing integers and 1 < idx(i) <= d.N(n)');
            end
            
            dsplit = d;
            dsplit.M = d.M + length(idx);
            dsplit.N = [d.N(1:n-1) diff([1 idx d.N(n)+1]) d.N(n+1: ...
                                                              end)];

            dsplit.trialidx = [d.trialidx(:,1:n-1) ...
                               [d.trialidx(1,n)+[1 idx]-1 ; ...
                                [d.trialidx(1,n)+idx-2 d.trialidx(2,n)]]];
            if n < size(d.trialidx,2),
                dsplit.trialidx = [dsplit.trialidx  d.trialidx(:,n+1:end)];
            end
            %                               d.trialidx(:,n:end)];
        end
        
        

        function dsel = select(d,sel)
        % dsel = select(d,sel)
        %
        % d -- an ridata object.
        % sel -- a set of indicies between 1 and d.D
        %
        % dsel is a new ridata with length(sel) trials, each
        % containing a single sample.

            if isempty(d.dT)
                dsel = ridata({d.Y(:,sel)}, ...
                              {d.U(:,sel)}, ...
                              {d.T(:,sel)});
            else
                dsel = ridata({d.Y(:,sel)}, ...
                              {d.U(:,sel)}, ...
                              {d.T(:,sel)}, ...
                              {d.X(:,sel)}, ...
                              {d.V(:,sel)}, ...
                              d.dT);
            end
            dsel = split(dsel,1,2:dsel.N(1));
        end
        
        function [dsel,sel] = unf_select(d,M,dist)
        % dsel = unf_select(d,M)
        %
        % d -- an ridata object.
        % M -- a number of data-points to select.
        %
        % Implements a uniform data space coverage algorithm to
        % select M data-points from d.  Each data-point becomes a
        % one sample trial in dsel.
        %
        % The number of dat points returned will be within +/- 10%
        % of M.
        % 
            
            if ~isempty(d.dT)
                D = [ d.Y ; d.U ; d.X ; d.V ];
            else
                D = [ d.Y ; d.U];
            end
            
%             if nargin < 3
%                 S = diag(1./std(D'));
%                 dist = @(x,A) sum(S*(A - repmat(x,1,size(A,2))).^2,1);
%             end

            f = @(gamma) sim_select(D',(1:size(D,2))',...
                                   gamma*ones(1,size(D,1)));
            

            l = 0; u = nthroot(M,size(D,1));
            [SEL,sel] = f(u);
            while length(sel) < M
                l = u;
                u = 2*l;
                [SEL,sel] = f(u);
            end
            

            
            while abs(length(sel)-M)./M > 0.1
                gamma = (u+l)/2;
                [SEL,sel] = f(gamma);
                if length(sel) > M
                    u = gamma;
                else
                    l = gamma;
                end
            end
            

          %   K = 50;
            
%             sel = randsample(size(D,2),1)*ones(M,1);
            
%             Ds = D(:,sel);
%             Ds = repmat(reshape(Ds,[size(D,1) 1 M]),[1 K 1]);
            
%             for i = 1:M
%                 %rsel = randsample(size(D,2),K);
%                 rsel = 1+floor(size(D,2)*rand(K,1));
%                 Dr = repmat(D(:,rsel),[1 1 M]);
%                 Di = min(sum((Ds - Dr).^2,1),[],3);
%                 [di,i0] = max(Di,[],2);
%                 sel(i) = rsel(i0);
%                 Ds(:,:,i) = repmat(D(:,rsel(i0)),[1 K]);
%             end

%             dists = dist(mean(D,2),D);
            
%             sel = zeros(1,M);

%             for i = 1:M
%                 [di,sel(i)] = max(dists);

%                 dists = min([ dists 
%                               dist(D(:,sel(i)),D)]);
%             end
            
            dsel = select(d,sel);
            
        end
        
        
        function [dscl,affine] = scale(d,affine)
            dscl = d;

            dscl.U = inv(affine.U.T)*(dscl.U - repmat(affine.U.mu,1,size(dscl.U,2)));
            dscl.Y = inv(affine.Y.T)*(dscl.Y - repmat(affine.Y.mu,1,size(dscl.Y,2)));
            if ~isempty(d.dT)
                dscl.X = inv(affine.X.T)*(dscl.X - repmat(affine.X.mu,1,size(dscl.V,2)));
                dscl.V = inv(affine.V.T)*(dscl.V - repmat(affine.V.mu,1,size(dscl.V,2)));
            end
        end
        
        
        function dsub = trials(dorig,trialnos)
            if ~all(trialnos == floor(trialnos)) || ...
                any(trialnos < 1)
                error('trialnos must be positive integers');
            end
            if any(trialnos > dorig.M)
                error('trialnos must not exceed dorig.M');
            end
            
            Y = cell(1,length(trialnos));
            U = Y; T = Y;
            
            if dorig.nx ~= 0, X = Y; V = Y; end
            for i = 1:length(trialnos)
                sel = dorig.trialidx(1,trialnos(i)):dorig.trialidx(2,trialnos(i));
                Y{i} = dorig.Y(:,sel);
                U{i} = dorig.U(:,sel);
                T{i} = dorig.T(:,sel);
                if dorig.nx ~= 0
                    X{i} = dorig.X(:,sel);
                    V{i} = dorig.V(:,sel);
                end
            end
            
            if dorig.nx == 0
                dsub = ridata(Y,U,T);
            else
                dsub = ridata(Y,U,T,X,V,dorig.dT);
            end
        end
    end
end