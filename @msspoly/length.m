function y=length(p)
if min(size(p)) == 0,
    y = 0;
else
    y = max(size(p));
end
end
