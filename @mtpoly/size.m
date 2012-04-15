function [x,y]=size(p,n)
% function [x,y]=size(p,n)

% AM 09.01.09

if nargin<2,
    if nargout==2,
        x=p.m;
        y=p.n;
    else
        x=[p.m p.n];
        y=[];
    end
else
    switch n,
        case 1,
            x=p.m;
        case 2,
            x=p.n;
        otherwise
            error('option not supported')
    end
end
