classdef (InferiorClasses = {?double}) msspoly
    properties 
        % Developer's Notes about Internal Representation:
        % 
        % Currently a sparse matrix representation is used.
        %
        % A large table is stored with unique rows for the distinct
        % monomials appearing in any matrix term.
        %
        %  e.g. for [ 2*x*y^3  0  x^7]
        %
        %  dim = [ 1 3 ]
        %  sub = [ 1 1 ; 1 3 ]
        %  vars = [ id(x) id(y) ; id(x) emptyid ]
        %  pow = [ 1 3 ; 7 0]
        %  coeff = [ 2 ; 1]
        %
        %  Though there is a canonical sorting for these entries not
        %  representend in this example.
        
        
        % Dimensions are stored here: dim == [numrow numcol]
        dim = [ 0 0 ];
        
        % E denotes the number of non-zero matrix entries.
        % v maximum number of variables appearing in any term.
        
        % Subscripts are E-by-2, sub(i,:) = [ rowno colno ] 
        % rows should be unique and 1 <= rowno <= numrow
        sub = zeros(0,2); 
        
        % Variables
        % E-by-v  where v is the maximal number of variables in any
        % term. The elements indicate the "id number" for variables
        % with a set aside id number for "no variable".
        var = zeros(0,1);
        
        % Powers
        % E-by-v list of non-negative integer powers for each variable.
        pow = zeros(0,1);
        
        % Conjugation
        % E-by-v list of variables which have been conjugated.
        cnj = zeros(0,1);
        
        % Coeff
        % E-by-1 list of complex coefficients appearing in any variable.
        coeff = zeros(0,1);
    end
    
    
    
    
    methods
        function p=msspoly(x,y,z,a,b)
            switch nargin,
              case 0,
              case 1,
                switch class(x),
                  case 'msspoly',
                    p = x;
                  case 'double',
                    if any(isnan(x(:)) | isinf(x(:)))
                        error('infinite coefficients not permitted');
                    end
                    p.dim = size(x);
                    [ii,jj,cc]=find(x);
                    p.sub = [ii(:) jj(:)];
                    p.var = zeros(size(p.sub,1),0);
                    p.pow = zeros(size(p.sub,1),0);
                    p.coeff = cc(:); 
                    p = p.make_canonical();
                  case 'char',
                    p = msspoly(x,1);
                  otherwise
                    error(['conversion of ' class(x) ' to msspoly not supported'])
                end
              case 2,
                if ~msspoly.isName(x)
                    error( ['Variable names must match n or Tn where n ' ...
                            'is a string, length(n) in 1:' ...
                            num2str(msspoly.nameLength) ...
                            ' and n(i) in ' msspoly.nameChars]);
                end
                if ~isa(y,'double') || ...
                        length(y) < 1 || ...
                        length(y) > 2
                    error(['2nd argument must be a double of length ' ...
                           'one or two.']); 
                end
                y=max(0,round(y(:)'));  % ERROR SUPRESSED: fix me
                m = y(1);
                if length(y) == 1, y(2) = 0; end
                p.dim = [ m 1 ];
                p.sub = [ (1:m)' ones(m,1) ];
                p.pow = ones(m,1);
                p.coeff = ones(m,1);                
                p.var = msspoly.name_to_id(x,y(2)+(1:y(1))'-1);
              case 5,
                p.dim = full(x);
                p.sub = full(y);
                p.var = full(z);
                p.pow = full(a);
                p.coeff = full(b);
                p = p.make_canonical();
              otherwise,
                error('Unsupported arguments.');
            end
            
            [flg,emsg] = check_canonical(p);
            if flg, error(emsg); end
        end
        
        function y = det(x)
          if ~isequal(size(size(x)), [1 2]) || size(x,1) ~= size(x,2)
            error('Must be a square matrix to compute determinant')
          end
          n = size(x,1);
          sigma = perms(1:n);
          I = eye(n);
          y = det(I(:,sigma(1,:)))*prod(diag(x.indexinto(1:n,sigma(1,:))));
          for i=2:size(sigma,1),
            y = y + det(I(:,sigma(i,:)))*prod(diag(x.indexinto(1:n,sigma(i,:))));
          end
        end
        
        function y = prod(x)
          if min(size(x)) == 1,
            y = x.indexinto(1);
            for i=2:length(x),
              y = y * x.indexinto(i);
            end
          else
            error('Prod not implemented for matrices yet');
          end
        end
        
        function y = adjugate(x)
          if ~isequal(size(size(x)), [1 2]) || size(x,1) ~= size(x,2)
            error('Must be a square matrix to compute adjugate')
          end
          n = size(x,1);
          y = zeros(n,n)*x.indexinto(1,1);
          s.type = '()';
          for i=1:n,
            for j=1:n,
              s.subs = {j, i};
              y = subsasgn(y,s,(-1)^(i+j)*det(x.indexinto([1:i-1 i+1:n],[1:j-1 j+1:n])));
            end
          end
        end
    end
    
    methods (Static)
        
        function p = zeros(sz,m)
            if nargin == 2
                sz = [sz m];
            end
            if ~spot_hasSize(sz,[1 2])
                error('Argument must be 1-by-2');
            elseif ~spot_isIntGE(sz,0)
                error('Argument must be non-negative integers');
            end
            
            p = msspoly(sz,zeros(0,2),zeros(0,1),zeros(0,1),zeros(0,1));
        end
        
    end
    
    methods (Static)
        function [b,xn] = isfreemsspoly(y)
            if ~isa(y,'msspoly'), 
                b = 0; 
                xn = []; 
            else
                [b,xn] = isfree(y);
            end
                
        end
        
        
        function [s1,s2] = padZeros(s1,s2)
            m1 = size(s1,2);
            m2 = size(s2,2);
            
            if m1 > m2
                s2 = [s2 zeros(size(s2,1),m1-m2)];
            elseif m2 > m1
                s1 = [s1 zeros(size(s1,1),m2-m1)];
            end
        end
        
        function ik = relate(x,y)
        %  ik = relate(x,y)
        %
        %  ik is N-by-2 where each the rows ik(j,:) = [ i k]  correspond
        %               to the pairs of indices with x(i) == y(k).
        %
        %  Throws an error if x or y cannot be cast to a double.
            x=full(double(x(:)));
            y=full(double(y(:)));

            mx=length(x);
            my=length(y);
            if mx*my==0, 
                ik=[]; 
            else
                m=mx+my;
                g=sortrows([[x  -(1:mx)'];[y (1:my)']]);
                g=[[g(1:m-1,1)==g(2:m,1);0] g(:,2)];
                ik=spot_gset(g);
            end
        end
        
% User facing name check.

        function ch = nameChars(idx)
            ch = ['@#_.' 'a'+(0:25)];
            
            if nargin > 0, ch = ch(idx); end
        end
        function l = nameLength()
            l = 4;
        end
        
        function mnp = maxNamePart()
            mnp = (length(msspoly.nameChars)+1).^(msspoly.nameLength);
        end

        function f = isName(ch)
            if isempty(ch) || ~ischar(ch)
                f = 0;
            elseif ch(1) == 'T'
                f = msspoly.isName(ch(2:end));
            else
                f = length(ch) <= msspoly.nameLength && ...
                    all(mss_match(msspoly.nameChars,ch)>0);
            end
            
        end
        
        function msk = isTrigId(vs)
            msk = logical(mod(vs,2));% == 1;
        end
        
        function n=name_to_id(ch,m)
            if (length(ch) > 1 & ch(1) == 'T')
                trigFlag = 1;
                ch = ch(2:end);
            else
                trigFlag = 0;
            end
            % This line copy pasted in id_to_name :\
            exponents = (length(msspoly.nameChars)+1).^(length(ch)-1:-1:0);
            namePart = sum(exponents.*(mss_match(msspoly.nameChars,ch)));
            
            maxId = floor((2^50)/msspoly.maxNamePart);
            
            if nargin == 1
                m = 1;
            elseif any(m > maxId)
                error(['Variable IDs greater than ' num2str(maxId) ...
                      ' not supported']);
            end

            n = trigFlag + 2*( namePart + msspoly.maxNamePart*m);
        end
        
        function name=id_to_name(id)
            trigFlag = mod(id,2);
            namePart = mod(id/2,msspoly.maxNamePart);
            m = round(id/2/msspoly.maxNamePart);

            % This line copy pasted in name_to_id :\
            exponents = (length(msspoly.nameChars)+1).^(msspoly.nameLength-1:-1:0);
            
            nameNo = mod(floor((namePart)./exponents),(length(msspoly.nameChars)+1));
            
            nameNo = nameNo(find(nameNo ~= 0,1,'first'):end);
            
            if isempty(nameNo), nameNo = 0; end

            name = [msspoly.nameChars(nameNo) num2str(m+1)];
            
            if trigFlag, name = [ 'T' name ]; end
        end
        
        function s=degree_to_string(d)
            if d~=1, s=['^' num2str(d)];
            else, s=''; end
        end
        
        function s1=term_to_string(s,t)
            
            m=length(t);
            k=round((m-1)/2);
            a=length(find(t(1:k)));      % index of first non-zero element in t(1:k)
            if a>1,                      % there are at least two terms
                ss=[msspoly.id_to_name(t(1)) msspoly.degree_to_string(t(1+k))];
                for i=2:a,
                    ss=[ss '*' msspoly.id_to_name(t(i)) msspoly.degree_to_string(t(i+k))];
                end
            elseif a==1,
                ss=[msspoly.id_to_name(t(1)) msspoly.degree_to_string(t(1+k))];
            else
                ss='';
            end
            if isempty(s), format = '%.5g';
            else, format = '%.5g'; end
            
            if isreal(t(m))
                coeff_str = sprintf(['(' format ')'],t(m));
            elseif real(t(m)) == 0
                coeff_str = sprintf(['(' format 'i)'],imag(t(m)));
            else
                coeff_str = sprintf(['(' format '+' format 'i)'],...
                                    real(t(m)), ...
                                    imag(t(m)));
            end
            
            if isempty(ss),
                if isempty(s),
                    tm=t(m);
                    s1=coeff_str;
                else
                    s1=[s '+' coeff_str];
                end
            else
                if t(m)==1,
                    if isempty(s),
                        s1=ss;
                    else
                        s1=[s '+' ss];
                    end
                elseif t(m)==-1,
                    if isempty(s),
                        s1=['-' ss];
                    else
                        s1=[s '-' ss];
                    end
                else
                    if isempty(s),
                        s1=[coeff_str '*' ss];
                    else
                        s1=[s '+' coeff_str '*' ss];
                    end
                end
            end
        end
        
        function z=match_list(x,y)
        % z = match(x)
        %
        % z is 1 if all elements of x match, 0 otherwise.
        %
        % z = match(x,y)
        %
        % x -- k-by-p matrix of unique elements.
        % y -- n-by-m matrix.
        %
        % Throws an error if x or y is not castable to double.
        %
        % Returns z -- n-by-m matrix with z(i,j) = l when
        %              y(i,j)=x(l) and zero otherwise.
        %
            xu = unique(x(:));            
            if nargin == 1
                if length(xu) == 1
                    z=1; 
                else
                    z= 0;
                end
            else
                if length(xu) ~= length(x)
                    error(['Elements of first argument must be ' ...
                           'unque']);
                end
                
                x = double(x);
                y = double(y);
                
                ik = msspoly.relate(y(:),x);
                if isempty(ik), z = zeros(size(y));
                else
                    [i,j] = ind2sub(size(y),ik(:,1));
                    z = sparse(i,j,ik(:,2),size(y,1),size(y,2));
                end
            end
        end
        
        
    end
    
    
    methods (Access = private)
        
        function [flg,emsg] = check_dimensions(p)
        % Tests that the internal representation has the correct
        % dimensions and types.
            emsg = [];
            if ~isa(p.dim,'double') || ~spot_hasSize(p.dim,[1 2])
                emsg = 'Dimensions must be a 1-by-2 double';
            elseif ~spot_isIntGE(p.dim,0)
                emsg = 'Dimensions must be non-neg. integers';
            elseif ~all(size(p.sub,1) == ...
                        [size(p.pow,1)  size(p.var,1)  size(p.coeff,1)])
                    emsg = 'Row count of sub/pow/var/coeff do not match.';
            elseif size(p.sub,2) ~= 2
                emsg = 'sub must be E-by-2';
            elseif size(p.coeff,2) ~= 1
                emsg = 'coeff must be E-by-1';
            elseif size(p.pow,2) ~= size(p.var,2)
                emsg = ['pow and var must have the same number of ' ...
                        'columns.'];
            end
            
            flg = ~isempty(emsg);
        end
        
        function [flag,emsg] = check_canonical(p)
        %
        %  Checks if a variety of invariants hold for the
        %  polynomial representation.
        %
        %  flag == 0   the values are in "canonical" format.
        %              
        %  flag == 1   the values can be made to be
        %              canonical (e.g. by sorting some rows).
        %
        %  flag == 2,  the values are illegal in a way that cannot
        %              be repaird.

            [errno,check_flag] = spot_mex_msspoly_check_canonical(p.dim,p.sub,...
                                                      p.var,p.pow,...
                                                      p.coeff);
            
            if errno ~= 0
                if check_flag
                    flag = 2;
                    emsg = ['msspoly internal check failed: errno ' ...
                            num2str(errno) ', see spot_mex_msspoly_check_canonical.'];
                else
                    flag = 1;
                    emsg = ['msspoly not canonical: errno is ' ...
                            num2str(errno) '. See spot_mex_msspoly_check_canonical.'];
                end
            else
                flag = 0;
                emsg = [];
            end
            
        end
        
        
        
        function q = removeEntries(p,I)
            q = p;
            q.sub(I,:) = [];
            q.var(I,:) = [];
            q.pow(I,:) = [];
            q.coeff(I,:) = [];
            q = q.make_canonical();
        end
            
        
        function p = make_canonical(p)
        % Takes as input q for which check_dimensions and
        % check_values return 0.
        %
        % Say ci are coeff, vi variables and pi powers generically.
        %
        % Simplifies ci*...*vi^0*... by setting vi = 0.
        %
        % Sorts entries of var in each row.
        %
        % Simplifies c1*...v1^p1*v1^p2... terms to v1^(p1+p2)
        %
        % Combines coefficients with identical var/pow
        % combinations.
        %
        % Removes:  0*v1^p1... terms.
        %
        % Sorts rows for increasing sub.

            function q = keepEntries(p,I)
                q = p;
                q.sub = p.sub(I,:);
                q.var = p.var(I,:);
                q.pow = p.pow(I,:);
                q.coeff = p.coeff(I,:);
            end

            [flg,emsg] = check_canonical(p);
            if flg == 2
                error(emsg);
            elseif flg == 0
                return;
            end
            
            p = keepEntries(p,p.coeff ~= 0);
            
            
            % Only the zero polynomial is left.
            if size(p.coeff,1) == 0
                return;
            end
            p.var = p.var.*(p.pow ~= 0); % Remove zero power
                                        % variables, vi^0
            [p.var,I] = sort(p.var,2,'descend'); % Sort for dec. id
            ii = repmat((1:size(p.pow,1))',1,size(p.pow,2));
            % TODO clean up this line below (applies sort to pow)
            p.pow = reshape(p.pow(sub2ind(size(p.pow),ii(:),I(:))),size(p.pow));
            I = [];
            ii = [];


            if size(p.var,2) ~= 0
                % Combine powers for variables repeated in a row.
                [p.pow,p.var,ktest] = ...
                    spot_mex_msspoly_make_canonical_combine_powers(p.pow,p.var);
                p.pow = p.pow(:,1:ktest);
                p.var = p.var(:,1:ktest);
            end
            v = size(p.pow,2);
            
            % Compute terms with identical (i,j) 
            [s,I] = sortrows([p.sub p.var p.pow]);

            % As a quick note, the operations in this MEX file seem
            % inefficient as they act row-wise instead of
            % column-wise (memory locality).
            %
            % A simple experiment showed a 10% speedup by switching
            % this order of operations, but transposing s is far
            % more costly.
            %
            %

            [k,s,p.coeff] = ...
                spot_mex_msspoly_make_canonical_combine_coeffs(s,p.coeff(I));

            p.coeff = p.coeff(1:k,:);
            p.sub = s(1:k,1:2);
            p.var = s(1:k,2+(1:v));
            p.pow = s(1:k,2+v+(1:v));
            s = [];
            
            % cmb = [1==0;~all(s(1:end-1,:) - s(2:end,:) == 0,2)] ;
            % acc = cumsum(cmb);
            % rtn = [1;1+find(acc(2:end) > acc(1:end-1))];
            % p.coeff = accumarray([1+acc ones(size(acc))],p.coeff(I));
            % p.pow = p.pow(I(rtn),:);
            % p.sub = p.sub(I(rtn),:);
            % p.var = p.var(I(rtn),:);

            p = keepEntries(p,p.coeff ~= 0);

            % Remove all zero cols of pow
            msk = ~all(p.pow == 0,1);
            p.var = p.var(:,msk);
            p.pow = p.pow(:,msk);

            [flg,emsg] = p.check_canonical();
            if flg
                error(emsg);
            end
        end
        
        function q = iter_binary(p,n,f)
            if ~spot_isIntGE(n,0) || ~spot_hasSize(n,[1 1])
                error('Second argument must be a non-negative integer.');
            end
            
            switch n
              case 0,
                q = msspoly(ones(size(p)));
              case 1,
                q = p;
              case 2,
                q = f(p,p);
              case 3,
                q = f(p,f(p,p));
              otherwise
                q = p;
                for i = 1:n-1
                    q = f(q,p);
                end
            end
        end
        
        
        function x = cat(d,varargin)
            x = msspoly();
            for k = 1:size(varargin,2)
                if ~isempty(varargin{k})
                    if isempty(x), x = msspoly(varargin{k}); 
                    else
                        x = cat_helper(x,varargin{k},d);
                    end
                end
            end
        end
        
        function q = cat_helper(p1,p2,d)
            p1 = msspoly(p1);
            p2 = msspoly(p2);
            
            od = 3-d;
            
            if size(p1,od) ~= size(p2,od), error('Incompatible dimensions.'); end
            
            [var1,var2] = msspoly.padZeros(p1.var,p2.var);
            [pow1,pow2] = msspoly.padZeros(p1.pow,p2.pow);
            
            dim = p1.dim;
            dim(d) = p1.dim(d) + p2.dim(d);
            off = [ 0 0 ];
            off(d) = p1.dim(d);
            
            q = msspoly(dim,...
                       [ p1.sub ; p2.sub+repmat(off,size(p2.sub,1),1) ],...
                       [ var1 ; var2 ],...
                       [ pow1 ; pow2 ],...
                       [ p1.coeff ; p2.coeff]);
        end
    end
    
    
    
    methods (Static)
        function [] = runAllTests()
        % Tests fall into two broad categories.
        % First is testing dimensionality issues, i.e.
        %  the correct sizes of the results.  Handling
        %  empty matrices including 0-by-m and n-by-0 matrices.
        %
        % Second is numerical tests, i.e. that the correct
        %   polynomial is represented by the result of arithmetic.
        %
        % The later tests generally make use of the functions:
        %    decomp,recomp,dmsubs
        % Whose behavior can be tested explicitly by the defining formulas.
            msspoly.runInternalTests();
            
            % Other tests rely on these functions working well,
            % please leave them to be run right after the internal tests.
            %msspoly.testIsfree();
            msspoly.testDecomp();
            msspoly.testDmsubs(); % Relies on decomp, run decomp
                                 % test first.
            msspoly.testRecomp(); % Tests rely on dmsubs, run dmsubs
                                 % test first.
            
            
            % Other functions.

            msspoly.testMTimes();
            
            
            
        end
        
        function [] = runInternalTests()
            msspoly.testMatch();
        end
    end
    
    methods (Static, Access = private)
        function [] = testDecomp()
        % Test decomp first for dimension issues then numerical
        % correctness.
            x = msspoly('x',4);
            
            e = 'Test Failed (decomp): dimension test.';
            [xx,pp,MM]=decomp(x);
            if ~spot_hasSize(pp,[4 4]), error(e); end
            if ~spot_hasSize(MM,[4 4]), error(e); end
            if ~spot_hasSize(xx,[4 1]), error(e); end
            
            % Assumes indexinto is working.
            x13 = indexinto(x,1:3); x4 = indexinto(x,4);            
            % This assumes plus is working.
            [xx,pp,MM]=decomp(x13+x4);
            if ~spot_hasSize(pp,[4 4]), error(e); end
            if ~spot_hasSize(MM,[3 4]), error(e); end
            if ~spot_hasSize(xx,[4 1]), error(e); end

            % This assumes mtimes / diag are working.
            [xx,pp,MM]=decomp(x13*x4);
            if ~spot_hasSize(pp,[3 4]), error(e); end
            if ~spot_hasSize(MM,[3 3]), error(e); end
            if ~spot_hasSize(xx,[4 1]), error(e); end
            
            
            % TODO: test arithmetical correctness (difficult)
        end
        
        function [] = testDmsubs()
        % Toy example:
            n = 10;
            K = 100;
            x = msspoly('x',n);
            x13 = indexinto(x,1:3); x4 = indexinto(x,4);            

            X1 = poissrnd(1,n,K);
            V1 = dmsubs(x,x,X1);
            V2 = dmsubs(x13+x4,x,X1);
            perm = randperm(n); xperm = indexinto(x,perm);
            [~,iperm] = sort(perm);
            V3 = dmsubs(xperm,x,X1);
            V4 = dmsubs(x,xperm,X1);
            
            edim = 'Test failed (dmsubs): dimension test failed.';
            eart = 'Test failed (dmsubs): arithmetic test failed.';
            if ~spot_hasSize(V1,size(X1)), error(edim);  end
            if ~all(all(X1 == V1)), error(eart); end

            if ~spot_hasSize(V2,[3 K]), error(edim); end
            if ~all(all(X1(1:3,:)+repmat(X1(4,:),3,1) == V2)), error(eart); end
            

            if ~spot_hasSize(V3,size(X1)), error(edim); end
            if ~spot_hasSize(V4,size(X1)), error(edim); end
            if ~all(all(X1(perm,:) == V3)), error(eart); end
            if ~all(all(X1(iperm,:) == V4)), error(eart); end
        end
        
        function [] = testRecomp()
        % Test recomp first for dimension issues then numerical
        % correctness.
            
            
        % Numerical correctness checks.
            e = 'Test Failed (decomp): arithmetical correctness.';
            
            % Dense version
            N = 1; m = 5; n = 5;
            M = poissrnd(1,N,m);
            p = poissrnd(1,m,n);
            x = msspoly('x',n);

            q = recomp(x,p,M);

            msspoly.compareAtPoints(x,q,@(x) M*prod((repmat(x',m,1).^p),2),100,e);
            
        end
        
        function [] = testMatch()
            x = ones(5,1);
            
            e = 'Test Failed (match): Single argument test failed.';
            if ~msspoly.match_list(x), error(e); end
            if ~msspoly.match_list(x'), error(e); end
            if ~msspoly.match_list(x*x'), error(e); end
            if ~msspoly.match_list(rand(1)*x'), error(e); end
            if ~msspoly.match_list(rand(1)*x'), error(e); end
            if msspoly.match_list(rand(5,1)*x'), error(e); end
            if msspoly.match_list(x*rand(1,5)), error(e); end
            
            e = 'Test Failed (match): Two argument test failed.';
            y = eye(10);
            x = (-5:5)';
            
            j = find(x == 1);
            
            z = msspoly.match_list(x,y);
            if any(diag(z) ~= j), error(e); end
            if diag(diag(z)) ~= diag(diag(j)), error(e); end
            
            y = repmat(10+(5:-1:1),10,1);
            x = 10+(1:10);
            z = msspoly.match_list(x,y);

            if ~all(all(repmat(z(1,:),10,1) == z)), error(e); end
            if ~all(z(1,:) == 5:-1:1), error(e); end
        end
        function [] = testMTimes()
            n=4;
            x = msspoly('x',n);
            x1=indexinto(x,1);
            x2=indexinto(x,2);

            eart = 'Matrix vector arithmetic test failed.';
            edim = 'Matrix vector dimension test failed.';
            op = x*x';
            ip = x'*x;
            if ~spot_hasSize(op,[4 4]), error(edim); end
            if ~spot_hasSize(ip,[1 1]), error(edim); end
            if ~isequal(x'*op*x,ip*ip), error(eart); end
            
            e = 'Scalar multiplication arithmetic test failed.';
            % Test scalar left multiplication
            if ~isequal(x1*x2,x2*x1), error(e); end
            if  ~isequal(4*x1,x1*4), error(e); end
            if ~isequal(4*x,x*4), error(e); end
            
            e = 'Scalar multiplication dimension test failed.';
            if ~spot_hasSize(4*x,size(x)), error(e); end
            
        end
        
        function [] = compareAtPoints(x,p,f,N,err)
        % Relies on dmsubs working appropriately.
        %
        %  x -- n-by-1 free msspoly including all variables p depends on
        %  p -- k-by-1 msspoly
        %  f -- @(x) (n-by-1) double --> (k-by-1) double
        %  N -- number of points to compare
        %  err -- error message if test fails.
        %
        %  Tests if polynomial and f evaluate to be identical on
        %  N integer valued random points.  
        %
        %  NOTE: Exact matching only makes sense if p has integer coefficients.
        %

            X = poissrnd(1,length(x),N);
            P = dmsubs(p,x,X);

            for i = 1:N
                if ~all(f(X(:,i)) == P(:,i)),
                    f(X(:,i))
                    P(:,i)
                    error(err); 
                end
                
            end
        end
        
        function p = randPoly(n,dim,dens,int)
            if nargin < 3
                dens = 4;
            end
            if nargin < 4
                int = 1;
            end
            x = msspoly('x',n);
            [~,xn] = isfree(x);
        % Pick random points in the matrix.
            K = floor(prod(dim)*dens);
            sub = ceil(repmat(dim,K,1).*rand(K,2));
            
            var = repmat(xn',K,1);
            pow = poissrnd(1,K,n);
            if int
                coeff = poissrnd(2,K,1);
            else
                coeff = randn(K,1);
            end
                
            p = msspoly(dim,sub,var,pow,coeff);
            
        end
    end
end
