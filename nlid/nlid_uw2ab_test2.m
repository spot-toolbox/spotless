function nlid_uw2ab_test2(T,n,d)
% function nlid_uw2ab_test2(T,n,d)
%
% test nlid_uw2ab.m by approximating the relation y=max{0,1-T|u|}
% (default T=3) for n samples (default n=300) on the interval {u}=[-1,1]  
% by trigonometric rational function of order d (default 4)

if nargin<1, T=3; end
if nargin<2, n=300; end
if nargin<3, d=4; end

uu=linspace(-1,1,n)';
vv=max(0,1-T*abs(uu));

tt=acos(uu);
z=msspoly('z',[1 1]);

[aa,bb,dd,L]=nlid_uw2ab(vv,tt,(1+z)^d,0);

vh=nlid_abdx2u(aa,bb,dd,uu);
vt=nlid_tpls(tt,vv,2*d,1);
fprintf('\n nlid_uw2ab error: bound=%f, true=%f  (%d terms) vs. %f\n', ...
        L,sqrt(mean((vh-vv).^2)),2*size(dd,2)-1,sqrt(mean((vt-vv).^2)))
close(gcf);
subplot(2,1,1);plot(uu,vv,uu,vh);grid
subplot(2,1,2);plot(uu,vv,uu,vt);grid
