function [A,B,C]=ltid_g2abc(g,m,o)
% function [A,B,C]=ltid_g2abc(g,m,o)
%
% INPUTS:
%   g  -  n-by-1 real
%   m  -  positive integer
%   o  -  option: "fire"(default), 'fir', 'fft'
%
% OUTPUTS:
%   A  -  m-by-m real Schur matrix
%   B  -  m-by-1 real
%   C  -  1-by-m real
%
% C(A^{t-1})B approximates g(t) 

if nargin<2, error('2 inputs required'); end
g=g(:);
n=length(g);
if nargin<3, o='fire'; end
n0=floor(n/2);
w=(0:n0)'*(2*pi/n);
z=exp(1i*w);

switch o,
    case 'fft',
        v=fft(g);
        v=v(1:n0+1)./z;
        G=ltid_vw2G(v,w,m,-2);
        [A,B,C]=ssdata(G);
        B=(eye(size(A,1))-A^n)\B;
        va=squeeze(freqresp(ss(A,B,C,0,-1),w));
        close(gcf);
        subplot(2,1,1);plot(w,real(v),'.',w,real(va));grid
        subplot(2,1,2);plot(w,imag(v),'.',w,imag(va));grid
        pause
    case 'fir',
        A=[zeros(1,n); eye(n-1) zeros(n-1,1)];
        B=[1;zeros(n-1,1)];
        C=g';
        Gr=reduce(ss(A,B,C,0,-1),m);
    case 'fire',
        q.type='()';
        q.subs={n:2*n-1};
        op.disp=0;
        op.issym=1;
        f=@(x)subsref(conv(g,x(n:-1:1)),q);
        %[V,D]=eigs(f,n,m,'LM',op);
        [V,D]=eigs(@(x)f(f(x)),n,m,'LA',op);
        A=V'*[zeros(1,m);V(1:n-1,:)];
        B=V(1,1:m)';
        C=g'*V;
    otherwise
        error(['method ' o ' not supported'])
end