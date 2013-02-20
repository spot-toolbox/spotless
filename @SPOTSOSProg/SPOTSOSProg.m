classdef SPOTSOSProg < SPOTSQLProg
    properties
        sosExpr = [];
    end
    methods (Access = private)
        function [flag,indet] = realPolyLinearInDec(pr,exp)
            [x,pow,Coeff] = decomp(exp);
            
            mtch = match(x,pr.variables);
            
            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)) | ...
                     any(imag(Coeff(:)) ~= 0));
            
            if ~flag, 
                indet = [];
            else
                mtch = match(pr.variables,x);
                indet = x(find(mtch==0));
                if any(istrig(indet))
                    flag = 0;
                    indet = [];
                end
            end
        end
        
        function [flag,tIn,pIn] = trigPolyLinearInDec(pr,expr)
            [x,pow,Coeff] = decomp(expr);
            
            mtch = match(x,pr.variables);
            
            % TODO: conj. symmetric check.
            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)));
            
            if ~flag, 
                indet = [];
            else
                mtch = match(pr.variables,x);
                indet = x(find(mtch==0));
                msk = istrig(indet);
                tIn = indet(find(msk));
                pIn = indet(find(~msk));
            end
        end
    end
    
    methods (Access = private, Static)
        function phi = buildGramBasis(expr,decvar)
            if ~spot_hasSize(expr,[1 1])
                error('buildGramBasis expects a scalar polynomial.');
            end
            
            [var,pow,M] = decomp(expr);
            mtch = match(var,decvar);
            b = 1:length(var);
            b(mtch(mtch ~= 0)) = [];
            indet = var(b);
            
            if length(indet) == 0
                phi = 1;
                return;
            end

            pow = pow(:,b);

            exponent_p_monoms = pow;
            csclasses={1:length(b)};
            exponent_m = monomialgeneration(exponent_p_monoms,csclasses);
    
            options = sdpsettings;
            options.verbose = 0;
            temp=sdpvar(1,1);
            tempops = options;
            tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
            tempops.verbose = 0;
            tempops.saveduals = 0;
            [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
            %disp('Reducing Monomials.');
            exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);
            exponent_m = exponent_m{1};
    
            phi = recomp(indet,exponent_m,eye(size(exponent_m,1)));
        end
    end
    
    methods (Access = protected)
        
        
        
        function [pr,Q,phi,y,basis] = buildSOSDecompPrimal(pr,expr)
            if ~spot_hasSize(expr,[1 1])
                error('buildSOSDecomp expects a scalar polynomial.');
            end

            decvar = pr.variables;

            phi = SPOTSOSProg.buildGramBasis(expr,decvar);

            [pr,Q] = pr.newPSD(length(phi));
    
            decvar = [decvar ; mss_s2v(Q)];
            sosCnst = expr-phi'*Q*phi;

            A = diff(sosCnst,decvar);
            b = subs(sosCnst,decvar,0*decvar);
            [var,pow,Coeff] = decomp([b A].');
            
            %[pr,y] = pr.withEqs(Coeff'*[1;decvar]);
            [pr] = pr.withEqs(Coeff'*[1;decvar]);
            y = [];
            basis = recomp(var,pow,eye(size(pow,1)));
        end
        
        function [pr,Q,phi,y,basis] = buildSOSDecompDual(pr,expr)
            if ~spot_hasSize(expr,[1 1])
                error('buildSOSDecomp expects a scalar polynomial.');
            end

            y = pr.variables;

            phi = SPOTSOSProg.buildGramBasis(expr,y);
            

            %  This next code requires that phi be the monomials
            %
            %  q in PSD,  C - A(y) in K,   D(y) + E(q) = f.
            %
            %  Note E has full row rank.  Pick q = q0 + G1.y + G2.z s.t.:
            %
            %  E(q0) = f,   D + E.G1 = 0,   E.G2 = 0.  
            %
            %  G2 w/ linearly indep. cols.
            %

            
            % Introuduce /dummy/ semidefinite variables.
            [~,Q] = pr.newPSD(length(phi));            
            q = mss_s2v(Q);
            
            decvar = [y ; q];
            sosCnst = phi'*Q*phi-expr;

            % A1*y + A2*q = b.
            A1 = diff(sosCnst,y);
            A2 = diff(sosCnst,q);
            b = -subs(sosCnst,[y;q],0*[y;q]);
            [var,pow,Coeff] = decomp([b A1 A2].');
            
            Coeff = Coeff.';
            
            f  = Coeff(:,1);
            D = Coeff(:,1+(1:length(y)));
            E = Coeff(:,1+length(y)+(1:length(q)));
        
            
            ny = length(y);
            m = length(f);
            nq = length(q);
            
            %  Now we'll use some of the structure of E.
            [row,col,s] = find(E);
            
            if ~all(col == (1:nq)') || length(unique(row)) ~= m
                error('Basis assumptions violated for SOS decomposition.');
            end
            
            % Scaling matrix.
            S = sparse(col,col,1./s,nq,nq);
            
            [row,I] = sort(row);
            S = S(:,I);

            
            % Find representatives for each column.
            [rr,ii,~] = unique(row);
            [~,I] = sort(rr);
            col_of_Id = col(ii(I));
            
            % E(i,col_of_ID(i)) = 1, so setting v(col_of_ID) = b, zeros o.w. makes
            % (Ev) = b.
            
            % Construct the particular solution.
            q0 = S*sparse(col_of_Id,ones(m,1),f,nq,1);
            
            % Construct G1 that is m-by-ny so that D = - E.G1
            G1 = S*sparse(repmat(col_of_Id,ny,1),...  
                        reshape(repmat(1:ny,m,1),[],1),...
                        -D(:),nq,ny);
            
            % Construct G2 with (nq - m) linearly independent rows.
            G2 = speye(nq);
            
            % Find when col(i) ~= col(i-1)
            paired = col_of_Id(row);
            %paired = [ col(2:end) ; col(end)];
            
            G2 = G2 - sparse(paired,col,ones(nq,1),nq,nq);
            G2(:,col_of_Id) = [];
            G2 = S*G2;

            if nq > m
                [pr,z] = pr.newFree(nq-m);
                pr = pr.withPSD(mss_v2s(q0 + G1*y + G2*z));
            else
                pr = pr.withPSD(mss_v2s(q0 + G1*y));
            end
        end
        
    end
    
    methods
        function pr = SPOTSOSProg(varargin)
            pr@SPOTSQLProg(varargin{:});
        end
        
        function [pr,tokens] = withSOS(pr,expr)
            if ~pr.realPolyLinearInDec(expr)
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                      'be linear in decision variables.']);
            end
            tokens = length(pr.sosExpr) + (1:prod(size(expr)));
            pr.sosExpr = [ pr.sosExpr ; expr(:)];
        end
        
        function [pr,y,basis] = withPolyEqs(pr,expr)
            if ~pr.realPolyLinearInDec(expr)
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                       'be linear in decision variables.']);
            end

            expr = expr(:);
            decvar = pr.variables;
            
            [indet,pow,M] = decomp(expr,decvar);
            
            monom = recomp(indet,pow,eye(size(pow,1)));
            
            [I,J,S] = find(M);
            
            [pr,y] = pr.withEqs(S);
            
            basis = monom(J);
        end
        
        function [pr,tokens] = withSOSMatrix(pr,expr)
            [lindec,indet] = pr.realPolyLinearInDec(expr);
            if ~lindec
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                      'be linear in decision variables.']);
            end
            
            if size(expr,1) ~= size(expr,2) || size(expr,1) == 0
                error('Expression must be a square non-empty matrix.');
            end
            
            x = msspoly('x',size(expr,1));
            expr = anonymize(expr,'y',indet);
            [pr,tokens] = withSOS(pr,x'*expr*x);
        end
        
        function [pr,tokens] = withUniTrigSOSMatrix(pr,expr)
            if size(expr,1) ~= size(expr,2) || size(expr,1) == 0
                error('Expression must be a square non-empty matrix.');
            end
            
            [lindec,trigIndet,polyIndet] = pr.trigPolyLinearInDec(expr);
            if ~lindec
                error(['Expression ' ...
                       'must be linear in decision variables.']);
            elseif ~isempty(polyIndet) || (length(trigIndet) ~= 1)
                error('Only a single, trigonometric indeterminate allowed.');
            end
            
            
            % How to handle these?  Special case? Maybe...

        end
        
        function n = numSOS(pr)
            n = length(pr.sosExpr);
        end
        
        function [pr,poly,coeff] = newFreePoly(pr,basis,n)
            if nargin < 3, n = 1; end
            [pr,coeff] = pr.newFree(length(basis)*n);
            poly = reshape(coeff,n,length(basis))*basis;
        end
        
        function sol = minimizeDualForm(pr,varargin)
            if nargin < 2, objective = msspoly(0); end

            Q = cell(pr.numSOS,1);
            phi = cell(pr.numSOS,1);
            y   = cell(pr.numSOS,1);
            basis   = cell(pr.numSOS,1);
            for i = 1:pr.numSOS
                pr = pr.buildSOSDecompDual(pr.sosExpr(i));
            end
            
            sol = minimizeDualForm@SPOTSQLProg(pr,varargin{:});
        end
        
        function sol = minimizePrimalForm(pr,varargin)
            if nargin < 2, objective = msspoly(0); end

            Q = cell(pr.numSOS,1);
            phi = cell(pr.numSOS,1);
            y   = cell(pr.numSOS,1);
            basis   = cell(pr.numSOS,1);
            for i = 1:pr.numSOS
                [pr,Q{i},phi{i},y{i},basis{i}] = pr.buildSOSDecompPrimal(pr.sosExpr(i));
            end

            sqlsol = minimizePrimalForm@SPOTSQLProg(pr,varargin{:});
            
            sol = SPOTSOSSoln(sqlsol,Q,phi,y,basis);
        end
    end
end