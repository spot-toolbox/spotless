function [p,emsg]=mssp_chk(m,n,s)
% function [p,emsg]=mssp_chk(m,n,s)
%
% checks that the triplet (m,n,s) defines a mssp polynomial

p=[];
emsg=[];
if ~isa(m,'double'), error('number of rows must be a double'); end
if ~isa(n,'double'), error('number of columns must be a double'); end
if ~isa(s,'double'), error('structure matrix must be a double'); end
if m~=max(0,round(m)), error('number of rows must be integer>=0'); end
if n~=max(0,round(n)), error('number of columns must be integer>=0'); end
s=full(s);
[ms,ns]=size(s);
if isempty(s), p.m=m;p.n=n;p.s=zeros(0,3); return; end
k=max(0,round((ns-3)/2));   % maximal number of variables per term
if ns~=2*k+3, 
    emsg='structure matrix must have odd>1 number of columns'; return;
end
if m<max(s(:,1)), 
    emsg='designated number of rows is too small'; return
end
if n<max(s(:,2)), 
    emsg='designated number of columns is too small'; return;
end
if ns<3, 
    emsg='structure matrix must have at least 3 columns'; return;
end


s=s(s(:,ns)~=0,:);                  % remove terms with zero coefficients
s(:,3:2+k)=s(:,3:2+k).*(s(:,3+k:2+2*k)>0);  % match vi/ki zeros   [!]
s(:,3+k:2+2*k)=s(:,3+k:2+2*k).*(s(:,3:2+k)>0);
s=mssp_srows(s);                    % order rows variable-wise
for i=1:k-1,                        % merge repeated variables term-wise
    b=(s(:,2+i)==s(:,3+i))&(s(:,2+i)>0);
    s(b,3+i+k)=s(b,3+i+k)+s(b,2+i+k);
    s(b,2+i)=0;
    s(b,2+i+k)=0;
end
                                    % [!] re-run matching vi/ki zeros
s=mssp_srows(s);                    % order rows variable-wise
s=s(:,max(abs(s),[],1)>0);          % remove all-zero columns from s
[ms,ns]=size(s);                    % re-calculate the size of s
if ns<3,
    p.m=m; p.n=n; p.s=zeros(0,3); return;
elseif ns==3,                           % constant polynomial
    [ii,jj,ss]=find(sparse(s(:,1),s(:,2),s(:,3)));
    s=[ii(:) jj(:) ss(:)];
else                                % combine terms
    s=sortrows(s);
    b=[all(s(1:ms-1,1:ns-1)==s(2:ms,1:ns-1),2);1==0];
    y=mss_gsum([b s(:,ns)]);        % [!] re-write for complex
    %y(:,2)
    %s(b==0,1:ns-1)
    s=[s(b==0,1:ns-1) y(:,2)];
    %for i=1:ms-1,
    %    if isequal(s(i,1:ns-1),s(i+1,1:ns-1)),
    %        s(i+1,ns)=s(i+1,ns)+s(i,ns);
    %        s(i,ns)=0;
    %    end
    %end
end
s=s(s(:,ns)~=0,:);                  % remove terms with zero coefficients
s=s(:,max(abs(s),[],1)>0);          % remove all-zero columns from s
if isempty(s), s=zeros(0,3); end

p.m=m;
p.n=n;
p.s=s;
