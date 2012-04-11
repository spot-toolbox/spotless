function nlid_yw2ab_test2(T,n,d)
% function nlid_yw2ab_test2(T,n,d)
%
% test nlid_yw2ab.m by approximating the relation y=max{0,1-T|u|}
% (default T=3) for n samples (default n=300) on the interval {u}=[-1,1]  
% by trigonometric rational function of order d (default 4)

if nargin<1, T=3; end
if nargin<2, n=300; end
if nargin<3, d=4; end

uu=linspace(-1,1,n)';
vv=max(0,1-T*abs(uu));
tt=acos(uu);
z=msspoly('z',[1 1]);

%[aa,bb,dd,L]=nlid_yw2ab(vv,tt,zeros(0,1),zeros(0,1),(1+z)^d,0);
[aa,bb,dd,L]=nlid_yw2ab(vv,tt,1,acos(0),(1+z)^d,0);

vh=nlid_abdx2u(aa,bb,dd,uu);
fprintf('\n nlid_yw2ab error: bound=%f, true=%f  (%d terms)\n', ...
        L,max(abs(vh-vv))^2,2*size(dd,2)-1)
close(gcf);
subplot(2,1,1);plot(uu,vv,uu,vh);grid
subplot(2,1,2);plot(vv,vh);grid
