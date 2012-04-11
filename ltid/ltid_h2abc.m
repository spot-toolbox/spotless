function [A,B,C]=ltid_h2abc(h,m,s,t)
% function [A,B,C]=ltid_h2abc(h,m,s,t)
%
% INPUTS:
%   h  -  n-by-k real
%   m  -  positive integer
%   s  -  1-by-d complex, with negative real part
%   t  -  n-by-1 positive increasing (default t=(0:n-1)')
%
% OUTPUTS:
%   A  -  m-by-m real Schur matrix
%   B  -  m-by-k real
%   C  -  1-by-m real
%
% C*exp(t(i)*A)*B approximates h(i) 

if nargin<2, error('2 inputs required'); end
[n,k]=length(h);
if nargin<, s=
if nargin<4, t=(0:n-1)'; end
