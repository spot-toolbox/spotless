function ltid_vxyu2Gdt_test(n,k,m,N)

if nargin<1, n=5; end
if nargin<2, k=3; end
if nargin<3, m=2; end
if nargin<4, N=100*n; end
A=zeros(2*n);
for i=1:n,
    a=0.1+rand;
    b=10*rand;
    A(2*i-1:2*i,2*i-1:2*i)=[-a b; -b -a];
end
B=randn(2*n,m);
C=randn(k,2*n);
D=randn(k,m);
%G0=c2d(ss(A,B,C,D),1,'tustin');
G0=ltid_rand(n,k,m);
[A,B,C,D]=ssdata(G0);
U=randn(m,N);
X=randn(2*n,N);
V=A*X+B*U+0.01*randn(2*n,N);
Y=C*X+D*U+0.01*randn(k,N);
G=ltid_vxyu2Gdt(V,X,Y,U,1e-5,1);
fprintf(' error: %f out of %f\n',norm(G-G0,Inf),norm(G0,Inf))
w=linspace(0,pi,200);
for i=1:k,
    for j=1:m,
        g0=squeeze(freqresp(G0(i,j),w));
        g=squeeze(freqresp(G(i,j),w));
        close(gcf)
        subplot(2,1,1);plot(w,real(g0),'.',w,real(g)); grid
        subplot(2,1,2);plot(w,imag(g0),'.',w,imag(g)); grid
        pause
    end
end
