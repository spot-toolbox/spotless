function ltid_reshape_test

N=1+round(5*rand);
m=1+round(3*rand);
n=1+round(4*rand);
G0=ltid_rand(N,m,n);
G=ltid_reshape(G0,m*n,1);
G=ltid_reshape(G,m,n);
fprintf(' reshape: m=%d, n=%d, error = %f\n',m,n,norm(G-G0))