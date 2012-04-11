function [R,dd]=dsgn_iq2p(a,b,d,n)
% function [R,dd]=dsgn_iq2p(a,b,d,n)
%
% given a in (0,1) and b in [0,1), find a polynomial 
%      r(x,y)=R*prod(repmat([x y],m,1).^dd,2) of degree d=[dx,dy],
% of degree d=[dx,dy] which approximates n^2  "uniform" samples of
% f(x,y)=angle(x+j*y)/(2*pi) on the set {(x,y): a<x<1, b*x<y<x}
% when n>0, minimizes L2 norm
% when n<0, minimizes Linf norm

if nargin<1, a=0.75; end
if nargin<2, b=0.5; end
if nargin<3, d=[3 3]; end
if nargin<4, n=-50; end

if (a<=0)||(a>=1), error('input 1 not in (0,1)'); end
if (b<0)||(b>=1), error('input 2 not in (0,1)'); end

N=n^2;
dd=mint_ch(mint_down(d));
m=size(dd,1);
x=repmat(linspace(a,1,abs(n))',1,abs(n));     
y=x.*repmat(linspace(b,1,abs(n)),abs(n),1);
x=x(:);                              % x samples
y=y(:);                              % y samples
f=(1/(2*pi))*angle(x+1i*y);            % f samples
M=repmat(x,1,m).^repmat(dd(:,1)',N,1);
M=M.*repmat(y,1,m).^repmat(dd(:,2)',N,1);
if n>0,
    R=(M\f)';
else
    R=ltid_Linf(M,f)';
end

r=M*R';
e=r-f;
fprintf(' (I,Q)->ph: max error=%e,  msq error=%e\n',max(abs(e)),std(e))
if nargout == 0,
    close(gcf);hist(e,abs(n))
end

