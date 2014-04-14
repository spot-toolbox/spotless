function pr=subsasgn(pr0,s,w)

switch s.type,
    case '.',
        switch s.subs,
            case 'free',
                pr=new(pr0,w,'free');
            case 'pos',
                pr=new(pr0,w,'pos');
            case 'lor',
                pr=new(pr0,w,'lor');
            case 'rlor',
                pr=new(pr0,w,'rlor');
            case 'psd',
                pr=new(pr0,w,'psd');
            case 'eq',
                pr=eq(pr0,w);
            case 'eqs',
                if (size(w,1)>1)&&(size(w,1)==size(w,2)),
                    pr=eq(pr0,mss_smm(w));
                end
            case 'sos',
                pr=sos(pr0,w);
            case 'sss',
                pr=sss(pr0,w);
            case 'sedumi',
                pr=sedumi(pr0,w);
            case 'sdpt',
                pr=sdpt3(pr0,w);
            otherwise
                error('option not supported');
        end
    otherwise
        error('option not supported');
end
