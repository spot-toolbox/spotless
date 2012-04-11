function n=mssp_id2n(ch,m)
% function n=mssp_id2n(ch,m)
%
% generate mssp polynomial variable number from its letter and order number

base=1000000;
if nargin<2,
    n=base*double(ch);
elseif m>=base, 
    error('variable indexes so large are not supported'); 
else
    n=(base*double(ch)+1)+m;
end

% old encoding scheme
%if nargin==2,
%    n=1000*m+10*(double(ch)-50);
%else
%    n=double(ch)-50;
%end