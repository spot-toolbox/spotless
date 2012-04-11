function s1=mssp_t2s(s,t)
% function s1=mssp_t2s(s,t)
%
% term to string conversion for mssp polynomial display

m=length(t);
k=round((m-1)/2);
a=length(find(t(1:k)));      % index of first non-zero element in t(1:k)
if a>1,                      % there are at least two terms
    ss=[mssp_n2id(t(1)) mssp_d2s(t(1+k))];
    for i=2:a,
        ss=[ss '*' mssp_n2id(t(i)) mssp_d2s(t(i+k))];
    end
elseif a==1,
    ss=[mssp_n2id(t(1)) mssp_d2s(t(1+k))];
else
    ss='';
end
if isempty(ss),
    if isempty(s),
        tm=t(m);
        s1=sprintf('%.5g',t(m));
    else
        s1=[s sprintf('%+.5g',t(m))];
    end
else
    if t(m)==1,
        if isempty(s),
            s1=ss;
        else
            s1=[s '+' ss];
        end
    elseif t(m)==-1,
        if isempty(s),
            s1=['-' ss];
        else
            s1=[s '-' ss];
        end
    else
        if isempty(s),
            s1=[sprintf('%.5g',t(m)) '*' ss];
        else
            s1=[s sprintf('%+.5g',t(m)) '*' ss];
        end
    end
end
        