function ltid_imfit_test(m,n,N,f)
% function ltid_imfit_test(m,n,N,f)
%
% test ltid_imfit.m on the frequency response of
%    H(s)=H1*(1/(s+j/5)+1/(s-j/5))+H2*(1/(s+j9/5)+1/(s-j9/5))
% where Hi are m-by-m real symmetric, using n samples from [1/2,3/2],
% and N poles from [0,1/4] and [7/4,2]

if nargin<1, m=2; end
if nargin<2, n=200; end
if nargin<3, N=4; end
if nargin<4, f=1; end

w=linspace(1/2,3/2,n)';
s=1i*w;
t=[linspace(0,1/4,N) linspace(7/4,2,N)]';
H1=randn(m); H1=reshape(H1'*H1,1,m^2);
H2=randn(m); H2=reshape(H2'*H2,1,m^2);
g=(1./(s+1i/5)+1./(s-1i/5))*H1+(1./(s+1i*9/5)+1./(s-1i*9/5))*H2;
g=g+(1i*0.2)*randn(size(g));

[c,v]=ltid_imfit(imag(g),w,t,f);

d=1;
for i=1:(2*N), d=min(d,min(eig(reshape(c(i,:),m,m)))); end
fprintf(' positivity check: %f>0\n',d)

ss=repmat(s,1,2*N);
pp=repmat(1i*t',n,1);
gh=0.5*((1./(ss-pp)+1./(ss+pp))*c);
fprintf(' zero real part check: %f=0\n',max(abs(real(gh(:)))))
fprintf(' max pointwise error: %f\n',max(abs(gh(:)-g(:))))
for i=1:m,
    for j=i:m,
        ij=i+m*(j-1);
        close(gcf);plot(w,imag(g(:,ij)),'.', ...
            w,imag(gh(:,ij)),w,v(:,ij));grid;pause
    end
end


