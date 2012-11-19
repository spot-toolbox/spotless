function q=subsref(p,s)

switch char(s.type),
  case '()',
    q = indexinto(p,s.subs{:});
  case '.',
    switch s.subs,
      case 'm',
        q=p.dim(1);
      case 'n',
        q=p.dim(2);
      case 'dim',
        q=p.dim;
      case 'sub',
        q=p.sub;
      case 'var',
        q=p.var;
      case 'pow',
        q=p.pow;
      case 'coeff',
        q=p.coeff;
      otherwise
        error('option not supported')
    end
  otherwise
    error('option not supported')
end
