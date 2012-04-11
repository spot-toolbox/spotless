function ltid_vwat2Gh_test(n,k,N)
% function aa=ltid_vwat2Gh_test(n,k,N)
%
% testing ltid_vwat2GH.m on a randomly generated k-by-k
% passive symmetric transfer matrix G of order 2*n, with N samples

if nargin<1, n=3; end
if nargin<2, k=2; end
if nargin<3, N=100; end

STREAM = RandStream.getDefaultStream;
reset(STREAM);

z=tf('z');
L=ltid_rand_pas(n,k);
L=L+L.';
S1=0.1*randn(k+1,k);
S2=0.1*randn(k+1,k);
L=L+(S1'*S1)*ss((z^2-1)/(z^2-2*cos(0.1)*z+1));
L=L+(S2'*S2)*ss((z^2-1)/(1-2*cos(2.8)*z+z^2));
w=linspace(0.5,2.5,N)';
v=reshape(permute(freqresp(L,w),[3 1 2]),N,k^2);
v=v+0.0001*randn(size(v));

aa=zeros(2*n+1,k);
for i=1:k,
    u=real(v(:,i*k-k+1:i*k));
    [aaa,bbb]=ltid_vw2ab(u,w,2*n,-1);
    aa(:,i)=aaa;
    e=u-ltid_wabc2v(w,aaa,bbb);
    fprintf(' column %d: error = %f\n',i,sqrt(max(sum(e.^2,2))))
end

%roots(ltid_a2q(aa))

%[G,H]=ltid_vwat2GH(v,w,aa);
t=[0;0.2;2.7;2.9];
[G,h]=ltid_vwat2Gh(v,w,aa,t,1);
H=ltid_th2H(t,h);
%a=ssdata(G);eig(a)
v1=reshape(permute(freqresp(G+H,w),[3 1 2]),N,k^2);
fprintf(' ltid_vwat2Gh max error: %f\n', ...
    max(abs(v1(:)-v(:))))

for i=1:k,
    for j=i:k,
        ij=i+k*(j-1);
        close(gcf)
        subplot(2,1,1);plot(w,real(v(:,ij)),'.',w,real(v1(:,ij)));grid
        subplot(2,1,2);plot(w,imag(v(:,ij)),'.',w,imag(v1(:,ij)));grid
        pause
    end
end