function q=subsasgn(p1,s,p2)
%
% TODO: Support logical indexing.
%
%
errindextype = ['Subscript indices must either be real positive ' ...
                'integers (logicals not currently supported).'];

    function tst = iscolon(arg)
    tst = length(arg) == 1 && ...
          isa(arg,'char') && ...
          strmatch(arg,':');
    end


switch s.type
  case '()',
    p1 = msspoly(p1);
    p2 = msspoly(p2);
    
    % Case with ':'.
    if iscolon(s.subs{1}) & size(s.subs,2) == 1
        if spot_hasSize(p2,[1 1])
            p2 = repmat(p2,size(p1));
        elseif prod(size(p2)) ~= prod(size(p1))
            error(['In an assignment  A(:) = B, the number of elements ' ...
                      'in A and B must be the same.']);
        end
        q = reshape(p2,size(p1));
        return;
    else

        switch length(s.subs)
          case 1,
            p2 = indexinto(p2,':');
            ind = s.subs{1};
            ind = ind(:);
            
            if isa(ind,'logical'),
                ind = find(ind);
            end
    

            
            if spot_hasSize(p2,[1 1])
                p2 = repmat(p2,length(ind),1);
            elseif length(p2) ~= length(ind)
                error(['In an assignment  A(I) = B, the number of elements ' ...
                       'in B and  I must be the same.']);
            end            
                        
            if ~spot_isIntGE(ind,1), 
                error(errindextype);
            end
            
            if any(ind(:) > prod(size(p1)))
                if any(size(p1) == 1)
                    d = find(size(p1) ~= 1);
                    if isempty(d), d = 2; end
                    p1.dim(d) = max(ind(:));
                else
                    error(['In an assignment  A(I) = B, a matrix A cannot ' ...
                           'be resized.']);
                end
            end
            dim = p1.dim;
            
          case 2,
            % We are given two subscript arrays.
            % First lets normalize the subscript arrays.
            if iscolon(s.subs{1}), is = (1:p1.dim(1))';
            else, is = s.subs{1}(:); end
            
            if iscolon(s.subs{2}), js = (1:p1.dim(2))'; 
            else, js = s.subs{2}(:); end
            
            if spot_hasSize(p2,[1 1])
                p2 = repmat(p2,[length(is) length(js)]);
            elseif ~spot_hasSize(p2,[length(is) length(js)])
                error(['Subscripted assignment dimension ' ...
                       'mismatch.']);
            end
            
            if ~spot_isIntGE(is,1) || ~spot_isIntGE(js,1)
                error(errindextype);
            end
            
            p1.dim(1) = max(max(is),p1.dim(1));
            p1.dim(2) = max(max(js),p1.dim(2));
%             if any(is > p1.dim(1)) || any(js > p1.dim(2))
%                 error('Index exceeds matrix dimensions.');
%             end

            % Now we have a list of row/col subscripts.
            % We want all pairs in order, one column at a time.
            Js = repmat(js',length(is),1);  
            Is = repmat(is,1,length(js));   

            p2 = indexinto(p2,':');
            ind = sub2ind(p1.dim,Is(:),Js(:));
          otherwise,
            error('Only two dimensional indexing supported.');
        end
        
        % Last instance of an index to appear is the value to assign.
        [I,uI] = unique(ind);
        
        if length(I) ~= length(ind)
            p2    = indexinto(p2,uI);
            ind   = I;
        end
        % At this point: p2 is k-by-1 msspoly and ind is k-by-1 unique
        % indices.

        p1ind = sub2ind(size(p1),p1.sub(:,1),p1.sub(:,2));
        ik = msspoly.relate(p1ind,ind);

        if ~isempty(ik)
            p1 = p1.removeEntries(ik(:,1)); % Equivalent to zeroing
                                                    % those entries.
        end
        
        % Remove entries of ind which correspond to p2 == 0.
        [i,j] = ind2sub(size(p1),ind);
        ij = [ i j ];
        
        [var1,var2] = msspoly.padZeros(p1.var,p2.var);
        [pow1,pow2] = msspoly.padZeros(p1.pow,p2.pow);

        q = msspoly(p1.dim,...
                   [ p1.sub ; ij(p2.sub(:,1),:) ],...
                   [ var1 ; var2 ],...
                   [ pow1 ; pow2 ],...
                   [ p1.coeff ; p2.coeff ]);
    end
  otherwise
    error('Unsuppported assignment type.');
end
    
p = p1;
end
