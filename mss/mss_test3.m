function h=mss_test3
% function mss_test3
%
% tests @mssprog/@msspoly environments by finding the minimal r
% such that p(x,y)=4*x^4*y^6+r*x^2-x*y^2+y^2 is a sum of squares

x=msspoly('x');                             % define the variables
y=msspoly('y');
r=msspoly('r');
pr=mssprog;                                 % initialize an MSS program
pr.free=r;                                  % register r as free decvar
pr.sos=4*(x^4)*(y^6)+r*(x^2)-x*(y^2)+y^2;   % register the sos constraint
pr.sedumi=r;                                % minimize r
h=pr({r});                                  % get the optimal r
