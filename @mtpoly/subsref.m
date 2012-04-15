function q=subsref(p,s)

switch char(s.type),
  case '()',
    q = indexinto(p,s.subs{:});
%         sp=p.s;
%         mp=p.m;
%         np=p.n;
%         sp(:,1)=sp(:,1)+mp*(sp(:,2)-1);    % essentially, p(:)
%         Np=mp*np;
%         pp=1:Np;
%         P=reshape(pp,mp,np);
%         Q=subsref(P,s);
%         [mq,nq]=size(Q);
%         qq=Q(:);
%         iq=repmat((1:mq)',1,nq); iq=iq(:);
%         jq=repmat(1:nq,mq,1); jq=jq(:);
%         ik=mss_relate(qq,sp(:,1));
%         sq=sp(ik(:,2),:);
%         sq(:,1:2)=[iq(ik(:,1)) jq(ik(:,1))];
%         q=mtpoly(mq,nq,sq);
    case '.',
        switch s.subs,
            case 'm',
                q=p.m;
            case 'n',
                q=p.n;
            case 's',
                q=p.s;
            otherwise
                error('option not supported')
        end
    otherwise
        error('option not supported')
end
