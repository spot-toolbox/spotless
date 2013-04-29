classdef spotsosprog < spotprog
    properties
        sosExpr = [];
        indeterminates = msspoly([]);
    end
    methods (Access = private)
        function flag = isPolyInIndet(pr,exp)
            mtch = match(pr.indeterminates,decomp(exp));
            flag =  all(mtch ~= 0);
        end
        
        function [flag,indet] = isRealPolyLinearInDec(pr,exp)
            [x,pow,Coeff] = decomp(exp);
            
            mtch = match(x,pr.variables);
            
            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)) | ...
                     any(imag(Coeff(:)) ~= 0));
            
            if ~flag,
                indet = [];
            else
                mtch = match(pr.variables,x);
                indet = x(find(mtch==0));
                if any(istrig(indet)) || ~pr.isPolyInIndet(indet)
                    flag = 0;
                    indet = [];
                end
            end
        end
        
        % function [flag,tIn,pIn] = trigPolyLinearInDec(pr,expr)
        %     [x,pow,Coeff] = decomp(expr);
            
        %     mtch = match(x,pr.variables);
            
        %     % TODO: conj. symmetric check.
        %     flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)));
            
        %     if ~flag, 
        %         indet = [];
        %     else
        %         mtch = match(pr.variables,x);
        %         indet = x(find(mtch==0));
        %         msk = istrig(indet);
        %         tIn = indet(find(msk));
        %         pIn = indet(find(~msk));
        %     end
        % end
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

            exponent_m = spot_build_gram_basis(pow);
    
            phi = recomp(indet,exponent_m,eye(size(exponent_m,1)));
        end
    end
    
    methods (Access = protected)
        function [pr,Q,phi,y,basis] = buildSOSDecompPrimal(pr,expr)
            if ~spot_hasSize(expr,[1 1])
                error('buildSOSDecomp expects a scalar polynomial.');
            end

            decvar = pr.variables;

            phi = spotsosprog.buildGramBasis(expr,decvar);

            [pr,Q] = pr.newPSD(length(phi));
    
            decvar = [decvar ; mss_s2v(Q)];
            sosCnst = expr-phi'*Q*phi;

            A = diff(sosCnst,decvar);
            b = subs(sosCnst,decvar,0*decvar);
            [var,pow,Coeff] = decomp([b A].');
            
            [pr,y] = pr.withEqs(Coeff'*[1;decvar]);
            basis = recomp(var,pow,eye(size(pow,1)));
        end
    end
    
    methods
        function pr = spotsosprog(varargin)
            pr@spotprog(varargin{:});
        end
        function [pr,x] = newIndeterminate(pr,name,number)
            x = msspoly(name,number);
            pr = pr.withIndeterminate(x);
        end
        
        function pr = withIndeterminate(pr,var)
            [vid,f] = isfree(var);
            if ~f,
                error('Indeterminates must be free msspolys.');
            end
            skip=match(pr.indeterminates,var);
            pr.indeterminates = [ pr.indeterminates ; var(skip == 0)];
        end
        
        function [pr,poly,coeff] = newSOSPoly(pr,basis,n)
            if nargin < 3, n = 1; end
            [pr,poly,coeff] = newFreePoly(pr,basis,n);
            pr = pr.withSOS(poly);
        end
        
        function [pr] = withSOS(pr,expr)
            if ~pr.isRealPolyLinearInDec(expr)
                error(['Coefficients must be real, indeterminates ' ...
                       'non-trigonometric, and expression must ' ...
                      'be linear in decision variables.']);
            end
            tokens = length(pr.sosExpr) + (1:prod(size(expr)));
            pr.sosExpr = [ pr.sosExpr ; expr(:)];
        end
        
        function [pr] = withPolyEqs(pr,expr)
            if ~pr.isRealPolyLinearInDec(expr)
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
        
        function [pr] = withSOSMatrix(pr,expr)
            [lindec,indet] = pr.isRealPolyLinearInDec(expr);
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
        
        function n = numSOS(pr)
            n = length(pr.sosExpr);
        end
        
        function [pr,poly,coeff] = newFreePoly(pr,basis,n)
            if nargin < 3, n = 1; end
            if ~isempty(basis) && n > 0
                if ~pr.isPolyInIndet(basis)
                    error('Basis must be polynomial in the indeterminates.');
                end
                [pr,coeff] = pr.newFree(length(basis)*n);
                poly = reshape(coeff,n,length(basis))*basis;
            else
                poly = [];
            end
        end
        
        
        function [sol,log] = minimize(pr,varargin)
            if nargin < 2, objective = msspoly(0); end

            Q = cell(pr.numSOS,1);
            phi = cell(pr.numSOS,1);
            y   = cell(pr.numSOS,1);
            basis   = cell(pr.numSOS,1);
            for i = 1:pr.numSOS
                [pr,Q{i},phi{i},y{i},basis{i}] = pr.buildSOSDecompPrimal(pr.sosExpr(i));
            end

            sol = minimize@spotprog(pr,varargin{:});
        end
    end
end