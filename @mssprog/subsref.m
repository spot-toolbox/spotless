function y=subsref(pr,s)

switch s.type,
    case '()',
        switch class(s.subs{1}),
            case 'msspoly',
                y=get(pr,s.subs{1});
            case 'cell',
                if isempty(pr.x),
                    y=[];
                else
                    y=double(pr,s.subs{1}{1});
                end
            otherwise
                error('only msspoly or cell arguments supported')
        end
    case '.',
        switch s.subs,
            case 'v',
                y=pr.v;
            case 't',
                y=pr.t;
            case 'o',
                y=pr.o;
            case 'n',
                y=pr.n;
            case 's',
                y=pr.s;
            case 'q',
                y=pr.q;
            case 'r',
                y=pr.r;
            case 'e',
                y=pr.e;
            case 'x',
                y=pr.x;
            case 'feasible',
                y=isempty(pr.x);
            otherwise
                error('unsupported option')
        end
    otherwise
        error('unsupported option')
end