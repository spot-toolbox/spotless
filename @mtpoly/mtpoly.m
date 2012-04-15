classdef (InferiorClasses = {?double}) mtpoly
    properties
        m=0;
        n=0;
        s=[];
    end
    methods
        function p=mtpoly(x,y,z)
        % function p=msspoly(x,y,z)
        %
        % constructor for msspoly class object p 
        % An msspoly object has 3 fields:  n,m,s, and represents a polynomial
        % p.m-by-p.n matrix P. Each row s(i,:)=[i,j,k1,...,km,d1,...,dm,c] of s
        % corresponds to a single term c*(v1^d1)*(v2*d2)*...*(vm^dm) in the
        % (i,j) entry of P, where vi is the variable with number ki
        % 
        % with 0 arguments, 
        %    p is the empty mtpoly object: p.n=p.m=0, p.s=[]
        % with 1 argument 
        %    (x 'double', 'mtpoly', or 'char')
        %    p is the mtpoly conversion of x
        % with 2 arguments 
        %    (x single character, y=[], y=a, or y=[a,b], a,b positive integers)
        %    x is a column vector of different independent variables;
        %    y=[]: p=x; y=a: p=[x0;...;x{a-1}]; y=[a,b]: p=[x{b};...;x{b+a-1)]
        % with 3 arguments (x,y positive integers, z an ms-by-ns,
        %    p.m=x,  p.n=y,  p.s=z
        
        % AM 09.01.09
            switch nargin,
              case 0,
                m=0; 
                n=0; 
                s=zeros(0,3); 
              case 1,
                switch class(x),
                  case 'mtpoly',
                    m=x.m;
                    n=x.n;
                    s=x.s;
                  case 'double',
                    if ~all(isfinite(x(:))),
                        error('infinite coefficients not permitted');
                    end
                    m=size(x,1);
                    n=size(x,2);
                    [ii,jj,cc]=find(x);
                    s=[ii(:) jj(:) cc(:)];
                  case 'char',
                    m=1;
                    n=1;
                    s=[1 1 mssp_id2n(x) 1 1];
                  otherwise
                    error(['conversion of ' class(x) ' to mtpoly not supported'])
                end
              case 2,
                if ~ischar(x), error('1st argument must be a string'); end
                x=x(:);
                x=x(1);
                if ~isa(y,'double'), 
                    error('2nd argument must be a double'); 
                end
                y=max(1,round(y(:)'));
                n=1;
                switch length(y),
                  case 0,
                    m=1;
                    s=[1 1 mssp_id2n(x) 1 1];
                  case 1,
                    m=y;
                    s=[(1:m)' ones(m,1) mssp_id2n(x,(0:m-1)') ones(m,2)];
                  case 2,
                    m=y(1);
                    s=[(1:m)' ones(m,1) mssp_id2n(x,(y(2):y(2)+m-1)') ones(m,2)];
                  otherwise
                    error('2nd argument must have less than 3 elements');
                end
              case 3,
                m=x;
                n=y;
                s=z;
            end
            p.m = m; 
            p.n = n; 
            p.s = s;
            [p,emsg]=mt_chk(p);
            if ~isempty(emsg),
                error(emsg);
            end
        end
    end
    methods (Access = private)
        function [p,emsg]=mt_chk(p)
        % function [p,emsg]=mssp_chk(m,n,s)
        %
        % checks that the triplet (m,n,s) defines a mssp polynomial
            m = p.m;
            n = p.n;
            s = p.s;

            emsg=[];
            if ~isa(m,'double'), error('number of rows must be a double'); end
            if ~isa(n,'double'), error('number of columns must be a double'); end
            if ~isa(s,'double'), error('structure matrix must be a double'); end
            if m~=max(0,round(m)), error('number of rows must be integer>=0'); end
            if n~=max(0,round(n)), error('number of columns must be integer>=0'); end
            s=full(s);
            [ms,ns]=size(s);
            if isempty(s), p.m=m;p.n=n;p.s=zeros(0,3); return; end
            k=max(0,round((ns-3)/2));   % maximal number of variables per term
            if ns~=2*k+3, 
                emsg='structure matrix must have odd>1 number of columns'; return;
            end
            if m<max(s(:,1)), 
                emsg='designated number of rows is too small'; return
            end
            if n<max(s(:,2)), 
                emsg='designated number of columns is too small'; return;
            end
            if ns<3, 
                emsg='structure matrix must have at least 3 columns'; return;
            end
            
            
            s=s(s(:,ns)~=0,:);                  % remove terms with zero coefficients
            s(:,3:2+k)=s(:,3:2+k).*(s(:,3+k:2+2*k)>0);  % match vi/ki zeros   [!]
            s(:,3+k:2+2*k)=s(:,3+k:2+2*k).*(s(:,3:2+k)>0);
            s=mssp_srows(s);                    % order rows variable-wise
            for i=1:k-1,                        % merge repeated variables term-wise
                b=(s(:,2+i)==s(:,3+i))&(s(:,2+i)>0);
                s(b,3+i+k)=s(b,3+i+k)+s(b,2+i+k);
                s(b,2+i)=0;
                s(b,2+i+k)=0;
            end
            % [!] re-run matching vi/ki zeros
            s=mssp_srows(s);                    % order rows variable-wise
            s=s(:,max(abs(s),[],1)>0);          % remove all-zero columns from s
            [ms,ns]=size(s);                    % re-calculate the size of s
            if ns<3,
                p.m=m; p.n=n; p.s=zeros(0,3); return;
            elseif ns==3,                           % constant polynomial
                [ii,jj,ss]=find(sparse(s(:,1),s(:,2),s(:,3)));
                s=[ii(:) jj(:) ss(:)];
            else                                % combine terms
                s=sortrows(s);
                b=[all(s(1:ms-1,1:ns-1)==s(2:ms,1:ns-1),2);1==0];
                y=mss_gsum([b s(:,ns)]);        % [!] re-write for complex
                                                %y(:,2)
                                                %s(b==0,1:ns-1)
                s=[s(b==0,1:ns-1) y(:,2)];
                %for i=1:ms-1,
                %    if isequal(s(i,1:ns-1),s(i+1,1:ns-1)),
                %        s(i+1,ns)=s(i+1,ns)+s(i,ns);
                %        s(i,ns)=0;
                %    end
                %end
            end
            s=s(s(:,ns)~=0,:);                  % remove terms with zero coefficients
            s=s(:,max(abs(s),[],1)>0);          % remove all-zero columns from s
            if isempty(s), s=zeros(0,3); end
            
            p.m=m;
            p.n=n;
            p.s=s;
            
        end
    end
end