function h=mss_test2
% function mss_test2
%
% find the minimum of polynomial q(x0,x1)=x0^4+x1^4-4*x0*x1

x=msspoly('x',2) ;
h=msspoly('h');
q=x(1)^4+x(2)^4-4*x(1)*x(2);

pr=mssprog;
pr.free=h;
pr.sos=q+h;
pr.sedumi=h;
h=pr({h});
if nargout<1, fprintf('\n Check: 2=%f\n',h); end