function ltid_vwqt2ph_test(n,k,N)
% function aa=ltid_vwqt2ph_test(n,k,N)
%
% testing ltid_vwqt2GH.m on a randomly generated k-by-k
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

[aa,bb]=ltid_vw2ab_psd(v,w,2*n,1);
v1=ltid_wabc2v(w,aa,bb);
fprintf(' ltid_vw2ab_psd max error: %f\n', ...
    max(abs(real(v1(:)-v(:)))))

q=ltid_a2q(aa);
t=[0;0.2;2.7;2.9];
[p,h]=ltid_vwqt2ph(v,w,q,t,1);
G=ltid_reshape(ltid_pq2G(p',q),k,k);
H=ltid_th2H(t,h);

v1=reshape(permute(freqresp(G+H,w),[3 1 2]),N,k^2);
fprintf(' ltid_vwat2GH max error: %f\n', ...
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