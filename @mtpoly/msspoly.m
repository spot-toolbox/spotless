function p=msspoly(x,y,z)
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

superiorto('double')
if (exist('TrigPoly'))% note: this is gross, but the only work-around short of converting this to the new class style (>= matlab v7.6)
  inferiorto('TrigPoly')  
end

switch nargin,
    case 0,
        m=0; 
        n=0; 
        s=zeros(0,3); 
    case 1,
        switch class(x),
            case 'msspoly',
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
                error(['conversion of ' class(x) ' to msspoly not supported'])
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
[p,emsg]=mssp_chk(m,n,s);
if isempty(emsg),
    p=class(p,'msspoly');
else
    error(emsg);
end
