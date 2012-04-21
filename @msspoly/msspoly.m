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
        
        % Coeff
        % E-by-1 list of complex coefficients appearing in any variable.
        coeff = zeros(0,1);
    end
    
    
    
    
    methods
        function p=msspoly(x,y,z,a,b)
        % function p=msspoly(x,y,z)
        %
        % constructor for msspoly class object p 
        % An msspoly object has 3 fields:  n,m,s, and represents a polynomial
        % p.m-by-p.n matrix P. Each row s(i,:)=[i,j,k1,...,km,d1,...,dm,c] of s
        % corresponds to a single term c*(v1^d1)*(v2*d2)*...*(vm^dm) in the
        % (i,j) entry of P, where vi is the variable with number ki
        % 
        % with 0 arguments, 
        %    p is the empty msspoly object: p.n=p.m=0, p.s=[]
        % with 1 argument 
        %    (x 'double', 'msspoly', or 'char')
        %    p is the msspoly conversion of x
        % with 2 arguments 
        %    (x single character, y=[], y=a, or y=[a,b], a,b positive integers)
        %    x is a column vector of different independent variables;
        %    y=[]: p=x; y=a: p=[x0;...;x{a-1}]; y=[a,b]: p=[x{b};...;x{b+a-1)]
        % with 3 arguments (x,y positive integers, z an ms-by-ns,
        %    p.m=x,  p.n=y,  p.s=z
        
        % AM 09.01.09
            switch nargin,
              case 0,
              case 1,
                switch class(x),
                  case 'msspoly',
                    p = x;
                  case 'double',
                    if ~all(isfinite(x(:))),
                        error('infinite coefficients not permitted');
                    end
                    p.dim = size(x);
                    [ii,jj,cc]=find(x);
                    p.sub = [ii(:) jj(:)];
                    p.var = zeros(size(p.sub,1),0);
                    p.pow = zeros(size(p.sub,1),0);
                    p.coeff = cc(:); 
                  case 'char',
                    if ~msspoly.hasSize(x,[1 1])
                        error(['Variable names must be a single ' ...
                               'character.']);
                    end
                    p.dim = [ 1 1 ];
                    p.sub = [ 1 1 ];
                    p.var = msspoly.name_to_id(x);
                    p.pow = 1;
                    p.coeff = 1;
                  otherwise
                    error(['conversion of ' class(x) ' to msspoly not supported'])
                end
              case 2,
                if ~ischar(x) || ~msspoly.hasSize(x,[1 1])
                    error(['Variable names must be a single ' ...
                           'character.']);
                end
                if ~isa(y,'double') || ...
                        length(y) < 1 || ...
                        length(y) > 2
                    error(['2nd argument must be a double of length ' ...
                           'one or two.']); 
                end
                y=max(1,round(y(:)'));  % ERROR SUPRESSED: fix me
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
    end

    
    methods (Static)
        
        function p = zeros(sz,m)
            if nargin == 2
                sz = [sz m];
            end
            if ~msspoly.hasSize(sz,[1 2])
                error('Argument must be 1-by-2');
            elseif ~msspoly.isIntGE(sz,0)
                error('Argument must be non-negative integers');
            end
            
            p = msspoly(sz,zeros(0,2),zeros(0,1),zeros(0,1),zeros(0,1));
        end
        
    end
    
    methods (Static, Access = private)
        function [b,xn] = isfreemsspoly(y)
            if ~isa(y,'msspoly'), 
                b = 0; 
                xn = []; 
            else
                [b,xn] = isfree(y);
            end
                
        end
        
        function flg = isIntGE(var,bnd)
            if nargin < 2, bnd = -Inf; end

            if ~isa(var,'double')
                flg = 0;
            elseif ~all(round(var) == var) | any(var < bnd)
                flg = 0;
            else
                flg = 1;
            end
        end
        function flg = hasSize(v,sz)
            if length(size(v)) ~= length(sz)
                flg = 0;
            else
                flg = all(size(v) == sz);
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
                ik=mss_gset(g);
            end
        end
        
        function n=name_to_id(ch,m)
            base=1000000;
            if nargin<2,
                n=base*double(ch);
            elseif m>=base, 
                error('variable indexes so large are not supported'); 
            else
                n=(base*double(ch)+1)+m;
            end
        end
        
        function name=id_to_name(id)
            base=1000000;
            m=floor(id/base);
            k=id-m*base;
            name=[char(m) repmat(num2str(k-1),1,(k>0))];
        end
        
        function s=degree_to_string(d)
            if d>1, s=['^' num2str(d)];
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
            if isempty(ss),
                if isempty(s),
                    tm=t(m);
                    s1=sprintf('%.5g',t(m));
                else
                    s1=[s sprintf('%+.5g',t(m))];
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
                        s1=[sprintf('%.5g',t(m)) '*' ss];
                    else
                        s1=[s sprintf('%+.5g',t(m)) '*' ss];
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
            if ~isa(p.dim,'double') || ~msspoly.hasSize(p.dim,[1 2])
                emsg = 'Dimensions must be a 1-by-2 double';
            elseif ~msspoly.isIntGE(p.dim,0)
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
        
        function [flg,emsg] = check_values(p)
        % Assumes that check_dimensions returns 0.
        % Tests the contents of the internal representation
        % for type consistency.
            [flg,emsg] = check_dimensions(p);
            if flg 
                return;
            end
            emsg = [];
            if ~msspoly.isIntGE(p.pow,0)
                emsg = 'pow must be non-negative integers.';
            elseif ~msspoly.isIntGE(p.var,0) % Strengthen to test valid var. names.
                emsg = 'var must be non-negative integers.';
            elseif ~msspoly.isIntGE(p.sub,1)
                emsg = 'sub must be positive integers';
            elseif any(p.sub(:,1) > p.dim(1)) ||...
                    any(p.sub(:,2) > p.dim(2))
                emsg = 'sub do not lie in legal range for dimensions.';
            elseif ~all(isfinite(p.coeff) | ~isreal(p.coeff))
                emsg = ['coeff must not be NaN/Inf/Complex ' ...
                        'numbers.'];
            end
            
            flg = ~isempty(emsg);
        end

        
        function [flg,emsg] = check_canonical(p)
        %  Checks that the representation has the following
        %  canonical form.
        %
        %  - There are no zero entries in coeff.
        %
        %  - All rows of var are in order of decreasing id.
        %
        %  - The only repeated elements of var are 0s.
        %
        %  - If var has a zero then pow has a corresponding zero.
        %
        %  - pow has no all-zero columns
        %
        %  - The rows are sorted by increasing values for the subscript.
        %
        %  - The rows of [ p.sub p.var p.pow ] are unique.

            [flg,emsg] = check_values(p);
            if flg
                return;
            end
            
            E = size(p.sub,1);
            v = size(p.pow,2);
            emsg = [];
            
            if any(p.coeff == 0)
                emsg = 'Zero coefficients found.';
            elseif v > 0
                if any(any(p.var(:,1:v-1) - p.var(:,2:v) < 0))
                    emsg = 'var is not sorted row-wise.';
                elseif any((p.var(:,1:v-1) ~= 0 ) & ...
                           (p.var(:,1:v-1) == p.var(:,2:v)))
                    emsg = 'var has repeated entries';
                elseif any(p.var(:) == 0 & p.pow(:) ~= 0)
                    emsg = 'pow non-zero for some zero in var';
                elseif size(p.pow,1) ~= 0 && any(all(p.pow == 0))
                    emsg = 'pow has a zero column';
                elseif E > 0 % Question to answer, how should this
                             % be sorted?
                    if any(p.sub(1:end-1,1) > p.sub(2:end,1)) |...
                        any((p.sub(1:end-1,1) == p.sub(2:end,1)) &...
                            (p.sub(1:end-1,2) > p.sub(2:end,2)))
                        emsg = ['sub not sorted in increasing ' ...
                                'order.'];
                    end
                end
            end
            flg = ~isempty(emsg);
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
            
            [flg,emsg] = check_values(p);
            if flg
                error(emsg);
            end

            cnz = p.coeff ~= 0;
            p = keepEntries(p,cnz);
            
            % Only the zero polynomial is left.
            if size(p.coeff,1) == 0
                return;
            end
            
            p.var = p.var.*(p.pow > 0); % Remove zero power
                                        % variables, vi^0
            
            [p.var,I] = sort(p.var,2,'descend'); % Sort for dec. id
            ii = repmat((1:size(p.pow,1))',1,size(p.pow,2));
            % TODO clean up this line below (applies sort to pow)
            p.pow = reshape(p.pow(sub2ind(size(p.pow),ii(:),I(:))),size(p.pow));
            
            v = size(p.pow,2);
            
            if size(p.var,2) ~= 0
                % Combine powers for variables repeated in a row.
                ROW = repmat((1:size(p.var,1))',1,size(p.var,2));
                v0  = [zeros(size(p.var,1),1) p.var];
                COL = 1+cumsum(diff(v0,[],2) < 0,2);

                p.pow = accumarray([ ROW(:) COL(:) ], p.pow(:));
                p.var = accumarray([ ROW(:) COL(:) ], p.var(:), size(p.pow), @max);
            end

            % Compute terms with identical (i,j) 
            [s,I] = sortrows([p.sub p.var p.pow]);
            cmb = [1==0;~all(s(1:end-1,:) - s(2:end,:) == 0,2)] ;
            acc = cumsum(cmb);
            rtn = [1;1+find(acc(2:end) > acc(1:end-1))];
            p.coeff = accumarray([1+acc ones(size(acc))],p.coeff(I));
            p.pow = p.pow(I(rtn),:);
            p.sub = p.sub(I(rtn),:);
            p.var = p.var(I(rtn),:);
            
            cnz = p.coeff ~= 0;
            p = keepEntries(p,cnz);

            % Remove all zero rows of pow
            msk = ~all(p.pow == 0,1);
            p.var = p.var(:,msk);
            p.pow = p.pow(:,msk);

            [flg,emsg] = p.check_canonical();
            if flg
                error(emsg);
            end
        end
        
        function q = iter_binary(p,n,f)
            if ~msspoly.isIntGE(n,0) || ~msspoly.hasSize(n,[1 1])
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
            if ~msspoly.hasSize(pp,[4 4]), error(e); end
            if ~msspoly.hasSize(MM,[4 4]), error(e); end
            if ~msspoly.hasSize(xx,[4 1]), error(e); end
            
            % Assumes indexinto is working.
            x13 = indexinto(x,1:3); x4 = indexinto(x,4);            
            % This assumes plus is working.
            [xx,pp,MM]=decomp(x13+x4);
            if ~msspoly.hasSize(pp,[4 4]), error(e); end
            if ~msspoly.hasSize(MM,[3 4]), error(e); end
            if ~msspoly.hasSize(xx,[4 1]), error(e); end

            % This assumes mtimes / diag are working.
            [xx,pp,MM]=decomp(x13*x4);
            if ~msspoly.hasSize(pp,[3 4]), error(e); end
            if ~msspoly.hasSize(MM,[3 3]), error(e); end
            if ~msspoly.hasSize(xx,[4 1]), error(e); end
            
            
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
            if ~msspoly.hasSize(V1,size(X1)), error(edim);  end
            if ~all(all(X1 == V1)), error(eart); end

            if ~msspoly.hasSize(V2,[3 K]), error(edim); end
            if ~all(all(X1(1:3,:)+repmat(X1(4,:),3,1) == V2)), error(eart); end
            

            if ~msspoly.hasSize(V3,size(X1)), error(edim); end
            if ~msspoly.hasSize(V4,size(X1)), error(edim); end
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
            if ~msspoly.hasSize(op,[4 4]), error(edim); end
            if ~msspoly.hasSize(ip,[1 1]), error(edim); end
            if ~isequal(x'*op*x,ip*ip), error(eart); end
            
            e = 'Scalar multiplication arithmetic test failed.';
            % Test scalar left multiplication
            if ~isequal(x1*x2,x2*x1), error(e); end
            if  ~isequal(4*x1,x1*4), error(e); end
            if ~isequal(4*x,x*4), error(e); end
            
            e = 'Scalar multiplication dimension test failed.';
            if ~msspoly.hasSize(4*x,size(x)), error(e); end
            
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
