function sim_ppdiff_test(m,d,n)
% function sim_ppdiff_test(m,d,n)
% 
% test sim_ppdiff.m on a pp vector of dimension m, 
% with degree d and n+1 breakpoints

if nargin<3, n=10; end
if nargin<2, d=4; end
if nargin<1, m=3; end

p.form='pp';
p.breaks=cumsum(rand(1,n+1));

p.coefs=randn(m*n,d+1);
p.pieces=n;
p.order=d+1;
p.dim=m;

q=sim_ppdiff(sim_ppint(p));

N=200; 
t0=p.breaks(1);
t1=p.breaks(n+1);
t=linspace(t0,t1,N+1);
pt=ppval(p,t);
qt=ppval(q,t);

for i=1:m,
    close(gcf);
    subplot(2,1,1);plot(t,qt(i,:),t,pt(i,:),'.');grid;
    subplot(2,1,2);plot(t,qt(i,:)-pt(i,:));grid
    pause
end
