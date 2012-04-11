function q=indexinto(p,varargin)
    sp=p.s;
    mp=p.m;
    np=p.n;
    sp(:,1)=sp(:,1)+mp*(sp(:,2)-1);    % essentially, p(:)
    Np=mp*np;
    pp=1:Np;
    P=reshape(pp,mp,np);
    Q=subsref(P,struct('type','()','subs',{varargin}));
    if isempty(Q),
        q = msspoly(0,0,[]); return;
    end    
    [mq,nq]=size(Q);
    qq=Q(:);
    iq=repmat((1:mq)',1,nq); iq=iq(:);
    jq=repmat(1:nq,mq,1); jq=jq(:);
    ik=mss_relate(qq,sp(:,1));
    sq=sp(ik(:,2),:);
    sq(:,1:2)=[iq(ik(:,1)) jq(ik(:,1))];
    q=msspoly(mq,nq,sq);
end