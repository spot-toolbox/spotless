function q=complex(p)
% function q=complex(p)
%
% simplification support for complex polynomials:
% utilizes identities $k*&k=1 and !^2=-1

[ms,ns]=size(p.s);
k=round((ns-3)/2);
vs=p.s(:,3:k+2);        % id array
ds=p.s(:,k+3:2*k+2);    % degree array
cs=p.s(:,ns);           % coefficients

% dealing with "!"=sqrt(-1)
Nj=mssp_id2n('!');      % sqrt(-1) symbol id number
b=(vs==Nj);  % presence of !^d
ds(b)=mod(ds(b),4);
e=any(b&(ds>1),2);
cs(e)=-cs(e);
ds(b)=mod(ds(b),2);

% dealing with "$" and "&"
N0=1000000;
N1=mssp_id2n('$');
N2=mssp_id2n('&');
vs1=vs-N1;
b1=(vs1>-1)&(vs1<N0);
vs2=vs-N2;
b2=(vs2>-1)&(vs2<N0);
f1=mss_unique(vs1(b1));
f2=mss_unique(vs2(b2));
ik=mss_relate(f1,f2);
ff=f1(ik(:,1));


q=msspoly(p.m,p.n,[p.s(:,1:2) vs ds cs]);