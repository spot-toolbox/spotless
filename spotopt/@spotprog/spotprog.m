classdef spotprog
%  TODO:
%
%  -- Every spot_decomp_linear should come with a potential error stmt.
%
%  -- Move /all/ variables with the same name so that cstr issued
%  variables can be renamed.
%
%
    properties 
        %
        % Basic Representation:
        %
        % Primal:
        %
        %    min (c1,f)+(c2,x), 
        %    x in K1,  
        %    A1*f+A2*x = b, 
        %    g-F1*f-F2*x = r in K2.
        %
        % Dual:
        %  
        %    max. (y,b) - (q,g),
        %    q in K2*,
        %    c-A2'y-(-F2'*q) = z in K1*,
        %    c1 = A1'*y-F1'*q
        %
        % variables = x(p).
        %
        K1 = struct('l',0,'q',[],'r',[],'s',[]);
        K2 = struct('l',0,'q',[],'r',[],'s',[]);
        
        % w = [x;z],  w(p(i)) corresponds to variable i.
        % Really, I want v(p(i)) corresponds to variable x(i).
        %
        freeVar = [];
        coneVar = [];   % variables for [K1] in the order they are created
        coneToVar = []; % for each coneVar, which elt. of K1 corresponds
        dualVar = [];
        
        coneCstrVar = [];
        coneToCstrVar = [];
        
        A = [];
        b = [];
        
        g = [];
        F = [];
        
        name = '@';
    end
    
    methods (Static)
        function x = varToDec(v,ordering)
            [~,I] = sort(ordering);
            x = v(I);
        end
        function opt = defaultOptions()
            opt = spot_sdp_default_options();
        end
        function n = psdDimToNo(d)
            n = spot_psdDimToNo(d);
        end
        function [d,err] = psdNoToDim(n)
            [d,err] = spot_psdNoToDim(n);
        end
        
        function flag = isScalarDimension(n)
            flag = spot_hasSize(n,[1 1]) && ...
                   spot_isIntGE(n,0);
        end
        function flag = isVectorDimension(n)
            flag = size(n,1) == 1 && spot_isIntGE(n,0);
        end
        function dim = checkCstrDimension(e,dim)
            if ~spotprog.isVectorDimension(dim)
                error('dim must be a vector of non-negative scalars.');
            end
            if size(e,2) ~= 1
                dim = repmat(dim,1,size(e,2));
                e = e(:);
            end
            if sum(dim) ~= length(e)
                error('Dimension and vector size disagree.');
            end
        end
        
        function [nl,nq,nr,ns] = coneDim(K)
            nl = K.l;
            nq = sum(K.q);
            nr = sum(K.r);
            ns = sum(spotprog.psdDimToNo(K.s));
        end
        
        function Pdiag = coneInnerProduct(K)
            [nl,nq,nr,ns] = spotprog.coneDim(K);
            
            Pdiag = ones(nl+nq+nr,1);
            
            ipMap = containers.Map('KeyType','uint32','ValueType','any');
            ssize = unique(K.s);
            for i = 1:length(ssize)
                n = ssize(i);
                I = repmat((1:n)',1,n);
                J = I';
                ipMap(n)=ones(nchoosek(n+1,2),1)+mss_s2v(I~=J);
            end
            
            ips = cell(1,length(K.s));
            for i = 1:length(K.s)
                ips{i} = ipMap(K.s(i));
            end
            if isempty(ips), d = [];
            else, d = vertcat(ips{:});
            end
            Pdiag = [ Pdiag ; d];
        end
        
        function [ol,oq,or,os] = coneOffset(K)
            ol = K.l;
            oq = ol+sum(K.q);
            or = oq+sum(K.r);
            os = or+sum(spotprog.psdDimToNo(K.s));
        end
    end
    
    methods ( Access = private )        
        function n = numVar(pr)
            n = length(pr.coneVar) + length(pr.coneCstrVar) + ...
                length(pr.dualVar) + length(pr.freeVar);
                
        end


        function n = numEquations(pr)
            n = size(pr.A,1);
        end
        
        

        
        function nm = varName(prog)
            nm = [prog.name 'v'];
        end
        
        function off = posOffset(pr)
            [off] = spotprog.coneOffset(pr.K1);
        end
        
        function off = lorOffset(pr)
            [~,off] = spotprog.coneOffset(pr.K1);
        end
        
        function off = rlorOffset(pr)
            [~,~,off] = spotprog.coneOffset(pr.K1);
        end
        
        function off = psdOffset(pr)
            [~,~,~,off] = spotprog.coneOffset(pr.K1);
        end
        
        function off = posCstrOffset(pr)
            [off] = spotprog.coneOffset(pr.K2);
        end
        
        function off = lorCstrOffset(pr)
            [~,off] = spotprog.coneOffset(pr.K2);
        end
        
        function off = rlorCstrOffset(pr)
            [~,~,off] = spotprog.coneOffset(pr.K2);
        end
        
        function off = psdCstrOffset(pr)
            [~,~,~,off] = spotprog.coneOffset(pr.K2);
        end

        function [pr,v] = insertVariables(pr,offset,number)
            pr.A = [ pr.A(:,1:offset) ...
                     sparse(size(pr.A,1),number) ...
                     pr.A(:,offset+1:end) ];
            
            pr.F = [ pr.F(:,1:offset) ...
                     sparse(size(pr.F,1),number) ...
                     pr.F(:,offset+1:end) ];
            
            v = msspoly(pr.varName,[number pr.numVar]);            
        end
        
        function [pr,v] = insertConeVariables(pr,offset,number)
            [pr,v] = insertVariables(pr,offset+pr.numFree,number);
            pr.coneVar = [ pr.coneVar ; v];
            shift = pr.coneToVar > offset;
            pr.coneToVar(shift) = pr.coneToVar(shift) + number;
            pr.coneToVar = [ pr.coneToVar ; offset+(1:number)'];
        end
        
        function pr = insertConstraints(pr,offset,Fnew,gnew)
            off = offset;
            pr.F = [ pr.F(1:off,:)
                     Fnew
                     pr.F(off+1:end,:)];
            pr.g = [ pr.g(1:off,:)
                     gnew
                     pr.g(off+1:end,:)];
        end
        
        function [pr,z] = addConstraintExpr(pr,off,e)
            [Fnew,gnew] = spot_decomp_linear(-e(:),[pr.freeVar;pr.coneVar]);
            [~,I] = sort(pr.coneToVar);
            nf = pr.numFree;
            pr = pr.insertConstraints(off,Fnew(:,[ (1:nf)' ; nf+I]),gnew);
            
            n = size(Fnew,1);

            z = msspoly(pr.varName,[n pr.numVar]);
            
            pr.coneCstrVar = [pr.coneCstrVar ; z];
            shift = pr.coneToCstrVar > off;
            pr.coneToCstrVar(shift) = pr.coneToCstrVar(shift) + n;
            pr.coneToCstrVar = [ pr.coneToCstrVar ; off+(1:n)' ];
        end
    end

    methods
        function pr = spotprog(name)
            if nargin < 1,
                name = '@';
            end
            
            if ~spot_hasSize(name,[1 1]) | ~msspoly.isName(name)
                error(['Argument must be a single character valid ' ...
                       'msspoly name.']);
            end
            pr.name = name;
        end
        
        function v = variables(pr)
            v = [pr.freeVar ; pr.coneVar];
        end
        
        function v = cstrVariables(pr)
            v = pr.coneCstrVar;
        end
        
        function nf = numFree(pr)
            nf = length(pr.freeVar);
        end
        
        
        function v = decToVar(prog,x)
            v = [ x(1:prog.numFree) 
                  x(prog.numFree+prog.coneToVar)];
        end
        
        function y = dualEqVariables(pr)
            y = pr.dualVar;
        end
        
        function [pr,y] = withEqs(pr,e)
            e = e(:);
            y = msspoly(pr.varName,[length(e) pr.numVar]);
            pr.dualVar = [ pr.dualVar ; y];
            [Anew,bnew] = spot_decomp_linear(e,[pr.freeVar;pr.coneVar]);
            nf = pr.numFree;
            [~,I] = sort(pr.coneToVar);
            pr.b = [ pr.b ; bnew];
            pr.A = [ pr.A ; Anew(:,[(1:nf)' ; nf+I]) ];

        end
        
        function [pr,z] = withCone(pr,e,K)
            [nl,nq,nr,ns]=spotprog.coneDim(K);
            [ol,oq,or,os]=spotprog.coneOffset(K);
            if nl+nq+nr+ns ~= size(e,1)
                error('Cone size does not match dim.');
            end

            if size(e,2) == 0, return; end
            
            if nl > 0, [pr,zl,yl] = pr.withPos(e(1:ol));
            else, zl = []; yl = []; end

            if nq > 0, [pr,zq,yq] = pr.withLor(e(ol+(1:oq)),K.q);
            else, zq = []; yq = []; end

            if nr > 0, [pr,zr,yr] = pr.withRLor(e(oq+(1:or)),K.r);
            else, zr = []; yr = []; end

            if ns > 0, [pr,zs,ys] = pr.withBlkPSD(e(or+(1:os)),K.s);
            else, zs = []; ys = []; end
            
            y = [ yl ; yq ; yr ; ys ];
            z = [ zl ; zq ; zr ; zs ];
            
        end
        
        
        function [pr,z] = withPos(pr,e)
            n = prod(size(e));
            off = pr.posCstrOffset;
            pr.K2.l = pr.K2.l + n;
            
            [pr,z] = addConstraintExpr(pr,off,e);
        end
        
        function [pr,z] = withLor(pr,e,dim)
            if nargin < 3,
                dim = size(e,1);
            end
            
            dim = spotprog.checkCstrDimension(e,dim);

            off = pr.lorCstrOffset;
            pr.K2.q = [pr.K2.q dim];
            
            [pr,z] = addConstraintExpr(pr,off,e);
        end
       
        function [pr,z] = withRLor(pr,e,dim)
            if nargin < 3,
                dim = size(e,1);
            end
            
            dim = spotprog.checkCstrDimension(e,dim);
            
            off = pr.rlorCstrOffset;
            pr.K2.r = [pr.K2.r dim];
            
            [pr,z] = addConstraintExpr(pr,off,e);
        end
        
        function [pr,z] = withBlkPSD(pr,e,dim)
            if nargin < 3,
                [dim,v] = spotprog.psdNoToDim(size(e,1));
                if v,
                    error(['expression size incorrect (not upper ' ...
                           'triangular of symmetric matrix).']);
                end
            end
            
            if ~spot_isIntGE(dim,0),
                error('dim must be non-negative integer.');
            end

            no = spotprog.psdDimToNo(dim);
            
            no = spotprog.checkCstrDimension(e,no);

            off = pr.psdCstrOffset;
            pr.K2.s = [pr.K2.s spotprog.psdNoToDim(no)];
            
            [pr,z] = pr.addConstraintExpr(off,e);
        end
        
        function [pr,z] = withPSD(pr,e)
            if size(e,1) ~= size(e,2)
                error('Arugment must be square.');
            end
            [pr,z] = pr.withBlkPSD(mss_s2v(e));
            z = mss_v2s(z);
        end
        
        function pr = withDD(pr,Q)
            if size(Q,1) ~= size(Q,2)
                error('Arugment must be square.');
            end
            
           ndim = spotprog.psdDimToNo(size(Q,1));
           % Constrain Q to be DD (vectorized version)
           n = length(Q);
           [pr,tauvec] = pr.newFree(ndim);
           tau = mss_v2s(tauvec);
                  
           T = diag(ones(1,n)); T = double(~T); 
           tau0d = T.*tau;
           T = triu(T);
           tauT = T.*tau;
           QT = T.*Q;
           tauT_vec = mss_s2v(tauT);
           QT_vec = mss_s2v(QT);
           nzs_tauT_vec = find(tauT_vec);
           nzs_QT_vec = find(QT_vec);
           
           % This should never happen, but I'm keeping the check in just
           % for debugging
           if ~all(nzs_tauT_vec == nzs_QT_vec)
               disp('Something strange happened');
               keyboard;
           end
           
           tauT_vec = tauT_vec(nzs_tauT_vec);
           QT_vec = QT_vec(nzs_QT_vec);
           pr = pr.withPos(tauT_vec - QT_vec); 
           pr = pr.withPos(tauT_vec + QT_vec);
           pr = pr.withPos(diag(Q) - tau0d*ones(n,1));
            
            
        end
        
        
        
        function pr = withSDD(pr,e)
            if size(e,1) ~= size(e,2)
                error('Arugment must be square.');
            end
            [pr,Q] = pr.newSDD(size(e,1));
            pr = pr.withEqs(Q - e);
            
        end
        
        function [pr,v] = newFree(pr,n,m)
            if nargin < 3, m = 1; end
            if ~spotprog.isScalarDimension(n)
                error('n must be a non-negative integer scalar.');
            end
            if ~spotprog.isScalarDimension(m)
                error('m must be a non-negative integer scalar.');
            end
            [pr,v] = pr.insertVariables(pr.numFree,n*m);
            pr.freeVar = [ pr.freeVar ; v];
            v = reshape(v,n,m);
        end
        
        function [pr,v] = newPos(pr,n,m)
            if nargin < 3, m = 1; end
            if ~spotprog.isScalarDimension(n)
                error('n must be a non-negative integer scalar.');
            end
            if ~spotprog.isScalarDimension(m)
                error('m must be a non-negative integer scalar.');
            end
            nl = pr.posOffset;
            pr.K1.l = pr.K1.l + n*m;
            [pr,v] = pr.insertConeVariables(nl,n*m);
            v = reshape(v,n,m);
        end

        function [pr,v] = newLor(pr,dim,m)
            if nargin < 3, m = 1; end
            if ~spotprog.isVectorDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            if ~spotprog.isScalarDimension(m)
                error('m must be a non-negative integer scalar.');
            end
            nq = pr.lorOffset;
            n = sum(dim);
            [pr,v] = pr.insertConeVariables(nq,n*m);
            pr.K1.q = [ pr.K1.q repmat(dim,1,m) ];
            v = reshape(v,n,m);
        end
        
        function [pr,v] = newRLor(pr,dim,m)
            if nargin < 3, m = 1; end
            if ~spotprog.isVectorDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            if ~spotprog.isScalarDimension(m)
                error('m must be a non-negative integer scalar.');
            end
            nr = pr.rlorOffset;
            n = sum(dim);
            [pr,v] = pr.insertConeVariables(nr,n*m);
            pr.K1.r = [ pr.K1.r repmat(dim,1,m) ];
            v = reshape(v,n,m);
        end
        
        function [pr,v] = newBlkPSD(pr,dim,m)
            if nargin < 3, m = 1; end
            if ~spotprog.isVectorDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            if ~spotprog.isScalarDimension(m)
                error('m must be a non-negative integer scalar.');
            end
            ns = pr.psdOffset;
            d = spotprog.psdDimToNo(dim);
            n = sum(d);
            [pr,v] = pr.insertConeVariables(ns,n*m);
            pr.K1.s = [ pr.K1.s repmat(dim,1,m) ];
            v = reshape(v,n,m);
        end
        
        function [pr,V] = newPSD(pr,dim)
            if ~spotprog.isScalarDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            [pr,v] = pr.newBlkPSD(dim);
            V = mss_v2s(v);
        end
        
        function [pr,Q] = newSym(pr, dim)
            if ~spotprog.isScalarDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            
            ndim = spotprog.psdDimToNo(dim);
            
            [pr,Qvec] = pr.newFree(ndim);
            Q = mss_v2s(Qvec);
        end
        
        function [pr,Q] = newDD(pr,dim)
            if ~spotprog.isScalarDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            ndim = spotprog.psdDimToNo(dim);
            
            [pr,Qvec] = pr.newFree(ndim);
            Q = mss_v2s(Qvec);
                        
                      
           % Constrain Q to be DD (vectorized version)
           n = length(Q);
           [pr,tauvec] = pr.newFree(ndim);
           tau = mss_v2s(tauvec);
                  
           T = diag(ones(1,n)); T = double(~T); 
           tau0d = T.*tau;
           T = triu(T);
           tauT = T.*tau;
           QT = T.*Q;
           tauT_vec = mss_s2v(tauT);
           QT_vec = mss_s2v(QT);
           nzs_tauT_vec = find(tauT_vec);
           nzs_QT_vec = find(QT_vec);
           
           % This should never happen, but I'm keeping the check in just
           % for debugging
           if ~all(nzs_tauT_vec == nzs_QT_vec)
               disp('Something strange happened');
               keyboard;
           end
           
           tauT_vec = tauT_vec(nzs_tauT_vec);
           QT_vec = QT_vec(nzs_QT_vec);
           pr = pr.withPos(tauT_vec - QT_vec); 
           pr = pr.withPos(tauT_vec + QT_vec);
           pr = pr.withPos(diag(Q) - tau0d*ones(n,1));
           
           
%            % Constrain Q to be DD (old vectorized version)
%            n = length(Q);
%            [pr,tau] = pr.newFree(n,n);
%            T = diag(ones(1,n)); T = double(~T);
%            tau = T.*tau;
%            Q0d = T.*Q;
%            pr = pr.withPos(tau - Q0d); 
%            pr = pr.withPos(tau + Q0d);
%            pr = pr.withPos(diag(Q) - tau*ones(n,1));
            
%            % Constrain Q to be DD (non-vectorized version)
%            n = length(Q);
%            for k = 1:n
%                 inds = [1:k-1,k+1:n];
%                 [pr,tau] = pr.newFree(n-1); % Slack variables
%                 pr = pr.withPos(tau - Q(k,inds)');
%                 pr = pr.withPos(tau + Q(k,inds)');
%                 pr = pr.withPos(Q(k,k) - sum(tau));
%            end
        end
        
        function [prog,A,pvar] = newSDD(prog,n)
            if ~spotprog.isScalarDimension(n)
                error('dim must be a row of non-negative scalars.');
            end
            [i,j]=ind2sub([n n],mss_s2v(reshape(1:n^2,n,n)));

            offdiag = i ~= j;
            i = i(offdiag);
            j = j(offdiag);

            M = nchoosek(n,2);

            [prog,L] = prog.newLor(3,M); % keyboard;
            prog = prog.withPos(L(1,:)+L(2,:));


            P = [ 1 1 0 ; 0 0 1 ; 1 -1 0]*L;
            I = [ i' ; i' ; j' ];
            J = [ i' ; j' ; j' ];

            Ind = sub2ind([n n],I(:),J(:));
            S=sparse(Ind,(1:3*M)',ones(3*M,1),n*n,3*M);

            %---- option a: these 30 lines.
            [pvar,p,Coeff] = decomp(P(:)); 
            M = S*Coeff;
            N = size(M,1);
            [I,J,C] = find(M);
            dim     = [ N 1 ];
            sub     = [ I ones(size(I)) ];
            coeff   = C;

            % Assumption above is that variables enter linearly,
            % and no affine term appears. As a result, each row corresponds
            % to a single variable, indicated by the column.
            [i,j,e] = find(p);

            if max(sum(p ~= 0,2)) ~= 1 || ...
                    min(sum(p~=0,2)) ~= 1 || ...
                    any(e ~= 1)
                error('Assumption violation.');
            end

            % Fetch variable ID numbers.
            [~,xn] = isfree(pvar);
            xn = xn(j);  % For each (i,j) find the variables for that row.
            xn(i) = xn;

            var = xn(J);
            pow = ones(size(J));

            SP = msspoly(dim,sub,var,pow,coeff);
            A = mss_v2s(mss_s2v(reshape(SP,n,n)));

            %--- option B, this one line
            % A = mss_v2s(mss_s2v(reshape(S*P(:),n,n)));
        end
        
        function [pr,Q] = newDDdual(pr,dim)
            if ~spotprog.isScalarDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end
            ndim = spotprog.psdDimToNo(dim);
            
            [pr,Qvec] = pr.newFree(ndim);
            Q = mss_v2s(Qvec);
                        
                      
           % Constrain Q to be in dual of DD (non-vectorized version)
           n = length(Q);
           pr = pr.withPos(diag(Q));
           for i = 1:n
               for j = [1:i-1, i+1:n] % j \neq i
                   [pr,tau(i,j)] = pr.newFree(1);
                   pr = pr.withPos(tau(i,j) - Q(i,j));
                   pr = pr.withPos(tau(i,j) + Q(i,j));
                   pr = pr.withPos(Q(i,i) + Q(j,j) - tau(i,j));
               end
           end
           
           
        end
        
        function [pr,Q] = newSDDdual(pr,n)
            if ~spotprog.isScalarDimension(n)
                error('dim must be a row of non-negative scalars.');
            end
            
           ndim = spotprog.psdDimToNo(n);
            
           [pr,Qvec] = pr.newFree(ndim);
           Q = mss_v2s(Qvec);
            
            % Constrain Q to be in dual of SDD (non-vectorized version)
           n = length(Q);
           M = [0.5000         0    0.5000
                0.5000         0   -0.5000
                0    1.0000         0];
           % pr = pr.withPos(diag(Q));
           for i = 1:n
               for j = [1:i-1, i+1:n] % j \neq i
                   % Xij = [Q(i,i) Q(i,j); Q(j,i) Q(j,j)];
                   % pr = pr.withPSD(Xij);
                   pr = pr.withLor(M*[Q(i,i);Q(i,j);Q(j,j)]);
               end
           end
            

           
        end
      
        
        function [pr,Q] = newDiag(pr,dim)
            if ~spotprog.isScalarDimension(dim)
                error('dim must be a row of non-negative scalars.');
            end

	    [pr,diagvec] = pr.newPos(dim);

	    Q = diag(diagvec);
        end
        
        
        function pred = isStandardDualForm(pr)
            pred = pr.numEquations == 0 && spotprog.coneDim(pr.K1) == 0;
        end
        
        function pred = isPrimalWithFree(pr)
            pred = 0 == pr.K2.l && ...
                   isempty([ pr.K2.q pr.K2.r pr.K2.s]);
        end
        
        
        
        
        function [sol] = minimize(prog,pobj,solver,options)
            if nargin < 2,
                pobj = 0;
            end
            if nargin < 3,
                solver = @spot_sedumi;
            end
            if nargin < 4,
                options = spotprog.defaultOptions;
            end
            
            pobj = msspoly(pobj);

            
            % Enable removal of redundant equations.
            %[feas,E,F,g,U,V,w,Ad,bd,cd,Kd] = spot_sdp_remove_redundant_eqs(A,b,c,K);
            
            if options.dualize
                if ~prog.isStandardDualForm
                    error(['Sorry dualization is not supported.  ' ...
                           'Please write your program in standard ' ...
                           'dual form.']);
                end
                [dl,dobj] = prog.toDual(pobj);
                [P,A,b,c,K,d] = dl.toSedumi(-dobj);
                
                [x,y,z,info] = solver(A,b,c,K,options);
                nf = dl.numFree;
                xsol = P*x;
                zsol = P(nf+1:end,nf+1:end)*z(nf+1:end);
                sol = spotprogsol(dl,-dobj,xsol,y,zsol,info,1);
            else
                pr = prog.primalize();
                nf = pr.numFree;
                [P,A,b,c,K,d] = pr.toSedumi(pobj);

                % Enable basic facial reduction.
                % save Abck_1sdsos_30.mat A b c K options P pr pobj 
                [x,y,z,info] = solver(A,b,c,K,options);

                if ~isempty(x)
                   xsol = P*x;
                else
                    xsol = [];
                end
                if isempty(z)
                    zsol = [];
                else
                    zsol = P(nf+1:end,nf+1:end)*z(nf+1:end);
                end
                
                sol = spotprogsol(pr,pobj,xsol,y,zsol,info);
            end
        end
        
        function [sol] = minimizeDSOS(prog,pobj,solver,options)
            if nargin < 2,
                pobj = 0;
            end
            if nargin < 3,
                solver = @spot_gurobi;
            end
            if nargin < 4,
                options = spotprog.defaultOptions;
            end
            
            pobj = msspoly(pobj);

            
            % Enable removal of redundant equations.
            %[feas,E,F,g,U,V,w,Ad,bd,cd,Kd] = spot_sdp_remove_redundant_eqs(A,b,c,K);
            
            if ~isfield(options,'dualize')
                options.dualize = false;
            end
            
            if options.dualize
                error('dualization not supported');
                if ~prog.isStandardDualForm
                    error(['Sorry dualization is not supported.  ' ...
                           'Please write your program in standard ' ...
                           'dual form.']);
                end
                [dl,dobj] = prog.toDual(pobj);
                [P,A,b,c,K,d] = dl.toSedumi(-dobj);
                
                [x,y,z,info] = solver(A,b,c,K,options);
                nf = dl.numFree;
                xsol = P*x;
                zsol = P(nf+1:end,nf+1:end)*z(nf+1:end);
                sol = spotprogsol(dl,-dobj,xsol,y,zsol,info,1);
            else
                pr = prog.primalize();
                nf = pr.numFree;
                [P,A,b,c,K,d] = pr.toSedumi(pobj);

                % Enable basic facial reduction.
                [x,y,z,info] = solver(A,b,c,K,options);

                xsol = P*x;
                zsol = []; % FIX THIS P(nf+1:end,nf+1:end)*z(nf+1:end);
                
                sol = spotprogsol(pr,pobj,xsol,y,zsol,info);
            end
        end
        
        function [sol] = minimizeSDSOS(prog,pobj,solver,options)
            if nargin < 2,
                pobj = 0;
            end
            if nargin < 3,
                solver = @spot_sedumi;
            end
            if nargin < 4,
                options = spotprog.defaultOptions;
            end
            
            pobj = msspoly(pobj);

            
            % Enable removal of redundant equations.
            %[feas,E,F,g,U,V,w,Ad,bd,cd,Kd] = spot_sdp_remove_redundant_eqs(A,b,c,K);
            
            if ~isfield(options,'dualize')
                options.dualize = false;
            end
            
            if options.dualize
                if ~prog.isStandardDualForm
                    error(['Sorry dualization is not supported.  ' ...
                           'Please write your program in standard ' ...
                           'dual form.']);
                end
                [dl,dobj] = prog.toDual(pobj);
                [P,A,b,c,K,d] = dl.toSedumi(-dobj);
                
                [x,y,z,info] = solver(A,b,c,K,options);
                nf = dl.numFree;
                xsol = P*x;
                zsol = P(nf+1:end,nf+1:end)*z(nf+1:end);
                sol = spotprogsol(dl,-dobj,xsol,y,zsol,info,1);
            else
                pr = prog.primalize();
                nf = pr.numFree;
                [P,A,b,c,K,d] = pr.toSedumi(pobj);

                % Enable basic facial reduction.
                [x,y,z,info] = solver(A,b,c,K,options);

                xsol = P*x;
                zsol = []; % FIX THIS!!! P(nf+1:end,nf+1:end)*z(nf+1:end);
                
                sol = spotprogsol(pr,pobj,xsol,y,zsol,info);
            end
        end
        
        
    end
    
    methods (Access = public)
        function [P,A,b,c,K,d] = toSedumi(pr,pobj)
            if ~pr.isPrimalWithFree()
                error(['Problem must be primal with free variables, ' ...
                       ' see toPrimalWithFree().']);
            end
            if ~spot_hasSize(pobj,[1 1])
                error('Objective must be scalar.');
            end
            
            [c,nd] = spot_decomp_linear(pobj,[pr.freeVar;pr.coneVar]);
            c1 = c(1:pr.numFree)';
            c2 = c(pr.numFree+1:end)';
            d = -nd;
            
            [~,I] = sort(pr.coneToVar);
            c2 = c2(I);
            
            
            projMap = containers.Map('KeyType','uint32','ValueType','any');
            ssize = unique(pr.K1.s);
            for i = 1:length(ssize)
                n = ssize(i);
                I = mss_s2v(reshape(1:n^2,n,n));
                projMap(ssize(i)) = sparse((1:length(I))',I,ones(size(I)),length(I),n^2);
            end
            
            nf = pr.numFree;
            [nl,nq,nr,ns] = spotprog.coneDim(pr.K1);

            projs = cell(1,length(pr.K1.s));
            for i = 1:length(pr.K1.s)
                projs{i} = projMap(pr.K1.s(i));
            end
            if isempty(projs), P = [];
            else, P = blkdiag(projs{:});
            end
            A = [ pr.A(:,1:nf+nl+nq+nr) pr.A(:,nf+nl+nq+nr+(1:ns))*P];
            
            c = [ c1 
                  c2(1:nl+nq+nr)
                  P'*c2((nl+nq+nr)+(1:ns),:)];
            
            b = pr.b;
            P = blkdiag(speye(nf+nl+nq+nr),P);
            K = pr.K1;
            K.f = pr.numFree;
        end
    end
end
