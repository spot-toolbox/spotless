classdef spotsosprg < spotsqlprg
    properties
       sosExpr = []; 
    end
    methods (Access = private)
        function flag = realPolyLinearInDec(pr,exp)
            [x,pow,Coeff] = decomp(exp);
            [~,xid] = isfree(x);
            [~,vid] = isfree(pr.variables);

            mtch = mss_match(xid,vid);

            flag = ~(any(any(pow(:,mtch(mtch~=0)) > 1)) | ...
                     any(imag(Coeff(:)) ~= 0));
        end
    end
    
    methods (Access = protected)
        function [pr,Q,phi] = buildSOSDecomp(pr,expr)
        % Here be dragons.
        %
            if ~spot_hasSize(expr,[1 1])
                error('buildSOSDecomp expects a scalar polynomial.');
            end

            decvar = pr.variables;

            [var,pow,M] = decomp(expr);
            [~,decvarid] = isfree(decvar);    
            [~,varid] = isfree(var);
            mtch = mss_match(varid,decvarid);
            b = 1:length(varid);
            b(mtch(mtch ~= 0)) = [];
            indet = var(b);
    
            pow = pow(:,b);

            exponent_p_monoms = pow;
            csclasses={1:length(b)};
            exponent_m = monomialgeneration(exponent_p_monoms,csclasses);
    
            options = sdpsettings;
            temp=sdpvar(1,1);
            tempops = options;
            tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
            tempops.verbose = 0;
            tempops.saveduals = 0;
            [aux1,aux2,aux3,LPmodel] = export(set(temp>0),temp,tempops);  
            disp('Reducing Monomials.');
            exponent_m = monomialreduction(exponent_m,exponent_p_monoms,options,csclasses,LPmodel);
            exponent_m = exponent_m{1};
    
            phi = recomp(indet,exponent_m,eye(size(exponent_m,1)));
            [pr,Q] = pr.newPSD(length(phi));
    
            decvar = [decvar ; mss_s2v(Q)];
            sosCnst = expr-phi'*Q*phi;

            A = diff(sosCnst,decvar);
            b = subs(sosCnst,decvar,0*decvar);
            [var,pow,Coeff] = decomp([b A].');
    
            pr = pr.withEqs(Coeff'*[1;decvar]);
        end
    end
    
    methods
        function pr = spotsosprg(varargin)
            pr@spotsqlprg(varargin{:});
        end
        
        function [pr,tokens] = withSOS(pr,expr)
            if ~pr.realPolyLinearInDec(expr)
                error(['Coefficients must be real, and expression ' ...
                       'must be linear in decision variables.']);
            end
            tokens = length(pr.sosExpr) + 1:prod(size(expr));
            pr.sosExpr = [ pr.sosExpr ; expr(:)];
        end
        
        function n = numSOS(pr)
            n = length(pr.sosExpr);
        end
        
        function [pr,poly,coeff] = newFreePoly(pr,basis,number)
            if nargin < 3, n = 1; end
            [pr,coeff] = pr.newFree(length(basis)*n);
            poly = reshape(coeff,n,length(basis))*basis;
        end
        
        function sol = optimize(pr,objective)
            Q = cell(pr.numSOS,1);
            phi = cell(pr.numSOS,1);
            
            for i = 1:pr.numSOS
                [pr,Q{i},phi{i}] = pr.buildSOSDecomp(pr.sosExpr(i));
            end
            
            sqlsol = optimize@spotsqlprg(pr,objective);
            
            sol = spotsossol(sqlsol,Q,phi);
        end
    end
end