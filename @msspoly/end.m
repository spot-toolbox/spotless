function index = end(p,pos,n)
    switch n,
      case 1,
        index = prod(size(p));
      case 2,
        index = size(p,pos);
      otherwise,
        error('msspoly accepts only two indices.');
    end
end