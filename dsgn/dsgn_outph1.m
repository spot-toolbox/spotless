function p=dsgn_outph1(a,b,m,n)
% function p=dsgn_outph(a,b,m,n)
%
% given 0<b<a, find polynomial p of degree less than abs(m) which gives 
% a good approximation to n samples (uniformly spaced on [0,a+b])  
% of the function  f(t)=acos(t/(2*a))/(2*pi)
%
% when m>0, minimizes L2 norm, 
% when m<0 minimizes Linf norm

if nargin<1, a=1; end
if nargin<2, b=0.9*a; end
if nargin<3, m=-4; end
if nargin<4, n=1000; end

if (a<=0)||(a<=b), error('incompatible first two inputs'); end

t=linspace(0,a+b,n)';                % argument samples
ft=(1/(2*pi))*acos(t/(2*a));   % samples of f
M=repmat(t,1,abs(m)).^repmat((abs(m)-1:-1:0),n,1);
if m>0,
    p=(M\ft)'; 
else
    p=ltid_Linf(M,ft)';
end
pt=polyval(p,t);                        % samples of p
er=max(abs(ft-pt));
fprintf('order=%d:  error = %e\n',abs(m)-1,er)
if nargout==0,
    close(gcf);
    subplot(2,1,1);plot(t,ft,t,pt);grid
    subplot(2,1,2);plot(t,ft-pt);grid
end

