function M = dsgn_lut1( f,m )
% function M = dsgn_lut1( f,m )
% 
% given a function f:[0,1)->R (applicable point-wise to data columns)
% and integer m>1, produce 3-by-N integer matrix M (where N=2^m)
% such that f(t)~M(1,i)/N+M(2,i)*(k-M(3,i))/N^3 
%       for i=1+fix(N*t),  k=1+fix(N*(N*t-i+1))

if nargin<1, f=@(x)(x./(1+x)); end % default f(t)=t/(1+t)
if nargin<2, m=8; end              % default m=8
N=2^m;                            
NN=N^2;                            % NN=N^2
M=zeros(3,N);
t=(0:(NN-1))'/NN;                  % input samples
y=f(t);                            % output samples
yy=reshape(y,N,N);                 % N-length segments of y, column-wise
tt=[repmat(1/N,N,1) t(1:N)];     
cc=tt\yy;                          % least squares fit
ee=yy-t(1:N)*cc(2,:);              % update the fit
cc(1,:)=(0.5*N)*(max(ee,[],1)+min(ee,[],1));
M(1,:)=round(cc(1,:));
M(3,:)=1+round(N*(M(1,:)-cc(1,:))./cc(2,:));
M(2,:)=round(N*cc(2,:));

i=1+fix(N*t);  
k=1+fix(N*(N*t-i+1));
yh=M(1,i)'/N+M(2,i)'.*(k-M(3,i)')/N^3;
h=[min(M,[],2) max(M,[],2)];
fprintf('LUT1: %d bits (0..%d range) max.error %e\n',m,N-1,max(abs(y-yh)))
fprintf('      principal value range:  %d..%d\n',h(1,:))
fprintf('      coefficient range:      %d..%d\n',h(2,:))
fprintf('      offset range:           %d..%d\n',h(3,:))

if nargout==0,
    close(gcf);plot(t,y-yh);grid
end

