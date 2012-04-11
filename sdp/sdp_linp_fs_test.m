function sdp_linp_fs_test

n=100;
N=500;
a=[0.5+rand(n,1) randn(n,1) randn(n,1)];
[t,s]=sdp_linp_fs(a);

ss=linspace(0,1,N);
y=min(repmat(a(:,1)-a(:,2)*t,1,N)-a(:,3)*ss,[],1);
close(gcf);plot(ss,y,[s s],[min(y) max(y)]);grid



