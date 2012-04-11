function p=dsgn_sqrt(a,m,n)
% function p=dsgn_sqrt(a,m,n)
%
% given a in (0,1), find polynomial p of degree less than abs(m) which
% gives a good approximation to n samples (uniformly spaced on [a,1])  
% of the function  f(t)=sqrt(t)
%
% when m>0, minimizes L2 norm, 
% when m<0 minimizes Linf norm

if nargin<1, a=0.5; end
if nargin<2, m=-4; end
if nargin<3, n=1000; end

if (a<=0)||(a>=1), error('input 1 not in the (0,1) range'); end

t=linspace(a,1,n)';                % argument samples
ft=sqrt(t);                        % samples of f
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

