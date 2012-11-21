function q=assign(p,value,index)
    s.type = '()';
    s.subs = {index};
    q = subsasgn(p,s,value);
end
