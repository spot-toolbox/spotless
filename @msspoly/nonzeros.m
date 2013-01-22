function S = nonzeros(p)
%  S = nonzeros(p)
%
%  p -- n-by-m msspoly
%  S -- k-by-1 msspoly
%
%  S consists of the non-zero entries of p.
    
    [~,~,S] = find(p);
end