function [p1,p2]=dsgn_outph2(x,m,n)
% function [p1,p2]=dsgn_outph2(x,m,n)
%
% given x such that 0 <= x(2)-x(1) <= x(3) < x(4) <= x(2)+x(1),
% find polynomials p1,p2 of order <|m| to approximate functions 
%       f1,f2: [(x(2)-x(1))^2,(x(2)+x(1))^2]->[0,pi/2]      
% defined by  
%       sqrt(t)=x(1)*exp(2*pi*j*f1(t))+x(2)*exp(-2*pi*j*f2(t))
% at n samples of t uniformly spaced on [(x(3)+1e-6)^2,x(4)^2]
%
% when m>0, minimizes L2 norm, 
% when m<0 minimizes Linf norm

if nargin<1, x=[1,2,0.9*(1+1),0.9*(1+2)]; end
if nargin<2, m=-5; end
if nargin<3, n=1000; end

if (x(1)>x(2))||(x(3)<x(2)-x(1))||(x(4)<=x(3))||(x(4)>x(2)+x(1)),
    error('inadmissible input 1')
end

t=linspace((x(3)+1e-6)^2,x(4)^2,n)';           % argument samples
h=(x(1)^2-x(2)^2+t)./(2*sqrt(t));
f1t=acos(h/x(1))/(2*pi);  % samples of f
f2t=acos((sqrt(t)-h)/x(2))/(2*pi);
%t=linspace(x(3)+1e-6,x(4),n)';           % argument samples
%h=(x(1)^2-x(2)^2+t.^2)./(2*t);
%f1t=acos(h/x(1))/(2*pi);  % samples of f
%f2t=acos((t-h)/x(2))/(2*pi);
M=repmat(t,1,abs(m)).^repmat((abs(m)-1:-1:0),n,1);
if m>0,
    p1=(M\f1t)'; 
    p2=(M\f2t)'; 
else
    p1=ltid_Linf(M,f1t)';
    if x(1)~=x(2),
        p2=ltid_Linf(M,f2t)';
    else
        p2=p1;
    end
end
p1t=polyval(p1,t);                        % samples of p
p2t=polyval(p2,t);
z=exp(1i*2*pi*p1t)*x(1)+exp(-1i*2*pi*p2t)*x(2);
er=max(abs(z-sqrt(t)));
%er=max(abs(z-t));
fprintf('\n outph2:    order=%d:  error = %e\n',abs(m)-1,er)
if nargout==0,
    close(gcf);
    subplot(2,1,1);plot(t,f1t,t,p1t);grid
    subplot(2,1,2);plot(t,f2t,t,p2t);grid
end



