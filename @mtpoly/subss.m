% Slow Variable Substitution
% q = subss(p,x,r)
%
% p -- p(x,...) a mtpoly
% x -- n-by-1 free mtpoly
% r -- n-by-1 mtpoly
%
% Returns
% q -- p(r,...)
%
function q = subss(p,x,r)
    if size(x,1) ~= size(r,1), 
        error(['2nd and 3rd argument do not agree.']); 
    end
    r = mtpoly(r);
    anon = mtpoly('#',length(x));
    p = subs(p,x,anon);
    x = anon;
    if length(x) > 1
        q = p;
        for i = 1:length(x)
            q = subss(q,indexinto(x,i),indexinto(r,i));
        end
    else
        [R,pw] = pdecomp(p,x);
        if length(pw) < 1
            q = p;
        else
            % Substitution by Horner's Rule
            q = indexinto(R,':',size(R,2));
            
            for i = length(pw):-1:2
                q = q*r^(pw(i)-pw(i-1)); % We should memoize powers of r.
                q = q + indexinto(R,':',i-1);
                % q = q + R(:,i-1);
            end
            if double(pw(1)) ~= 0
                q = q*r^pw(1);
            end
            q = reshape(q,size(p,1),size(p,2));
        end
    end
end
