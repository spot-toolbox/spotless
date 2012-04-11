function display(p)
% function display(p)
%
% full matrix display for mss polynomials

% AM 09.01.09

if p.m*p.n==0,
    fprintf('   Empty msspoly matrix: %d-by-%d\n',p.m,p.n);
    return
end
S{p.m,p.n}=[];             % strings of entries
[ns,ms]=size(p.s);
for i=1:ns,
    ii=p.s(i,1);
    jj=p.s(i,2);
    S{ii,jj}=mssp_t2s(S{ii,jj},p.s(i,3:ms));
end
L=zeros(1,p.n);
for i=1:p.n,
    L(i)=3;
    for j=1:p.m,
        k=length(S{j,i});
        if k==0,
            S{j,i}=' 0';
        else
            L(i)=max(L(i),k);
        end
    end
end
fprintf('\n')
for i=1:p.m,
    fprintf('[  ');
    for j=1:p.n,
        fprintf('%s%s  ',repmat(' ',1,L(j)-length(S{i,j})),S{i,j});
    end
    fprintf(']\n');
end
fprintf('\n')
        