function [bb,L]=ltid_vwa2b(v,w,aa,o)
% function [bb,L]=ltid_vwa2b(v,w,aa,o)
%
% fitting real part of v to trigonometric ratio b/a 
%
% INPUTS:
%   v   - mv-by-nv complex 
%   w   - mv-by-1 from [0,pi]
%   aa  - (m+1)-by-1 real
%   o   - display parameter
%
% OUTPUTS:
%   bb  - (m+1)-by-1 real
%   L   - real>=0: lower bound on point-wise approximation quality |e|
%   
% for a(t)=cos(t*(0:m))*aa, b(t)=cos(t*(0:m))*bb, minimizes 
%    L=max_k|v(k,:)-b(w(k))/a(w(k))|


if nargin<3, error('3 inputs required'); end
if nargin<4, o=0; end
if ~isa(v,'double'), error('input 1 not a double'); end
[mv,nv]=size(v);
if ~isa(w,'double'), error('input 2 not a double'); end
if ~isreal(w), error('input 2 not real'); end
[mw,nw]=size(w);
if nw~=1, error('input 2 not a column'); end
if mw~=mv, error('inputs 1,2 have different number of rows'); end
if ~isa(aa,'double'), error('input 2 not a double'); end
if size(aa,2)~=1, error('input 3 not a column'); end
m=size(aa,1)-1;

vmx=max(abs(v(:)));
if vmx>0, v=v/vmx; end            % normalize
cs=cos(w*(0:m));                  % samples of trigonometric functions
aw=cs*aa;                         % aw(i)=a(w(i))
rva=real(v).*repmat(aw,1,nv);

pr=mssprog;
B=msspoly('B',nv*(m+1));          % coefficients of b
B=reshape(B,m+1,nv);
pr.free=B;                        % B is unconstrained
bw=cs*B;                          % bw(i,:)=b(w(i))
L=msspoly('L');
pr.free=L;
x=msspoly('x',mv*(nv+1));         % rotated Lorentz cone variables
x=reshape(x,nv+1,mv);             % columns are individual cones          
pr.lor=x;                         % x(1,i)>|x(2:nv+1,i)| 
pr.eq=L*aw-x(1,:)';               % x(1,i)=L*aw(i)
pr.eq=rva-bw-x(2:nv+1,:)';        % x(2:nv+1,i)'=v(i,:)*aw(i)-bw(i,:)
pr=sedumi(pr,L,o>0);              % optimize using SeDuMi
bb=vmx*pr({B});
L=vmx*pr({L});