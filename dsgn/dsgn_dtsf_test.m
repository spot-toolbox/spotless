function dsgn_dtsf_test(u0,w0,m,s,n)

if nargin<1, u0=pi/2; end
if nargin<2, w0=3*pi/4; end
if nargin<3, m=4; end
if nargin<4, s=2; end
if nargin<5, n=100; end

u=linspace(0,u0,n)';
w=linspace(w0,pi,n)';
d=1:(s*m);
d=[0 d(d~=s*round(d/s))];

[h,r]=dsgn_dtsf(d,u,w);
t=linspace(0,pi,10*n)';
ht=(cos(t*d)*h).^2;
p0=-10*log10(sum(ht(t>w0))/sum(ht));
p1=-10*log10(max(ht(t>w0))/min(ht(t<u0)));
p2=-20*log10(r);
fprintf('\n %2.1fdB (%2.1fdB~%2.1fdB) for %d taps\n',p0,p1,p2,2*m)
h=h/(2*h(1));
k=length(h);
h=[h(k:-1:2);2*h(1);h(2:k)];
d=[-d(k:-1:2) 0 d(2:k)];
%ii=(1:(m*s));g=(pi/s)*ii;g=sin(g)./g;g=[g(length(g):-1:1) 1 g];
%ii=[-ii(length(ii):-1:1) 0 ii];

close(gcf)
subplot(2,1,1); bar(d,h); grid
subplot(2,1,2); semilogy(t/pi,ht,u/pi,ones(size(u)), ...
    w/pi,repmat(r^2,n,1)); grid