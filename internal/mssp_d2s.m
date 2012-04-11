function s=mssp_d2s(d)
% function s=mssp_d2s(d)
%
% degree to string conversion for mssp display

if d>1,
    s=['^' num2str(d)];
else
    s='';
end