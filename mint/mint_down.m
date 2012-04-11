function y=mint_down(x)
% function y=mint_down(x)
%
% all rows of double integer x>=0, and their projections onto the
% coordinate planes

if ~mint_isint(x), error('input not a double integer'); end
if any(x(:)<0), error('input has negative elements'); end
[mx,nx]=size(x);
if mx*nx==0, error('input is empty'); end

my=sum(prod((x>0)+1,2));
y=zeros(my,nx);
jj=1:nx;
iy=1;                    % current row in y
ix=0;                    % current row in x
xx=zeros(1,nx);          % modified row
while 1,                 % go over the rows of x
    j=max((xx>0).*jj);   % maximal j such that xx(j)>0     
    if j==0,             % done with this xx
        ix=ix+1;
        if ix>mx,        % completely done
            break;
        else
            xx=x(ix,:);
        end
    else
        y(iy,:)=xx;
        iy=iy+1;
        xx(j)=0;
        xx(jj>j)=x(ix,jj>j);
    end
end
y=mss_unique(y);     % remove repeated rows