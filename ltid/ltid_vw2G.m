function [G,L]=ltid_vw2G(v,w,n,o,m)
% function [G,L]=ltid_vw2G(v,w,n,o,m)
%
% INPUTS:
%   v   - mv-by-nv complex 
%   w   - mv-by-1 from [0,pi]
%   n   - desired model order
%   m   - desired number of inputs (default m=1)
%   o   - algorithm/display flag (o<0 means silence, default o=1)
%
% OUTPUTS:
%   G  - stable DT system with m inputs and k=nv/m outputs
%   L  - lower bound for the error
%        
% |o|=1: no passivity, cheaper  - use vw2abc, a2q, vwq2G, reshape
% |o|=2: no passivity, advanced - use vw2ab, vwq2G, reshape
% |o|=3: symmetric passive      - use vw2ab_psd, ab2G, reshape

if nargin<3, error('3 inputs required'); end
if ~isa(v,'double'), error('input 1 not a double'); end
[mv,nv]=size(v);
if ~isa(w,'double'), error('input 2 not a double'); end
if ~isreal(w), error('input 2 not real'); end
[mw,nw]=size(w);
if nw~=1, error('input 2 not a column'); end
if mw~=mv, error('inputs 1,2 have different number of rows'); end
if ~isa(n,'double'), error('input 3 not a double'); end
n=max(1,round(real(n(1))));
if nargin<4, o=1; end
if nargin<5, 
    if abs(o)<3, 
        m=1; 
    else
        m=round(sqrt(nv));
    end
end
if ~isa(m,'double'), error('input 3 not a double'); end
m=max(1,round(real(m(1))));
k=round(nv/m);
if nv~=k*m, error('inputs 1,4 incompatible'); end
ix=1:nv;
switch abs(o),
    case 1,
        [aa,bb,cc,L]=ltid_vw2abc(v,w,n,o>0);
        q=ltid_a2q(aa);
        p=ltid_vwq2p(v,w,q);
    case 2,
        [aa,bb,L]=ltid_vw2ab(v,w,n,o>0);
        q=ltid_a2q(aa);
        p=ltid_vwq2p(v,w,q);
    case 3,
        if k~=m, error('not a square transfer matrix'); end
        [aa,bb,L]=ltid_vw2ab_psd(v,w,n,o>0);
        [p,q]=ltid_ab2pq(aa,bb);
        ix=mss_s2v(reshape(1:k^2,k,k));
end
G=ltid_pq2G(p,q);
g=squeeze(freqresp(G,w));
G=ltid_reshape(G,k,m);
if o>0, 
    L=sqrt(L);
    amin=ltid_tpmin(aa);
    er=sqrt(max(sum(abs(v(:,ix)-g(:,ix)).^2,2)));
    ere=sqrt(max(sum(real(v(:,ix)-g(:,ix)).^2,2)));
    eim=sqrt(max(sum(imag(v(:,ix)-g(:,ix)).^2,2)));
    fprintf(' error: %f>%f (re=%f, im=%f) check %f>0\n',er,L,ere,eim,amin)
    ltid_chk_vwg(v(:,ix),w,g(:,ix)); 
end