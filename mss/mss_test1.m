function h=mss_test1(c)
% function h=mss_test1(c)
%
% find the minimum of polynomial p(t)=c(1)*t+...+c(n)*t^n on [-0.5,0.5]

if nargin<1, c=[0.1 1 0 -3 0 3]; end;        % default input
c=c(:); n=length(c);                  % problem dimension
x=msspoly('x');                       % independent variables
h=msspoly('h');
q=0; for i=n:-1:1, q=q+c(i)*((1+x^2)^(n-i))*(x^i); end
pr=mssprog;                           % initialize mss program
pr.free=h;                            % register h as a free decvar
pr.sos=q+h*(1+x^2)^n;
pr.sedumi=h;
h=-pr({h});

if nargout<1,
    tt=linspace(-0.5,0.5,200);            % verify the result graphically
    cc=polyval([c(n:-1:1)' 0],tt);
    close(gcf);plot(tt,cc,[-0.5 0.5],[h,h]);grid
end


