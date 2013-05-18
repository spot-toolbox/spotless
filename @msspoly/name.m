function [nm,err] = name(p)
    [f,id] = isfree(p);
    if f
        for i = 1:length(p)
            nm{i} = msspoly.id_to_name(id(i));
        end
    elseif nargin < 2
        error('msspoly is not free.');
    else
        err = 1; nm = [];
    end
    
end