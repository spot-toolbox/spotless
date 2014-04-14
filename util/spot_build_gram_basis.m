function [mpow,coeffIndxPos,coeffIndxZero] = spot_build_gram_basis(pow,mpow)
  
    pow = full(pow);
    
    if nargin < 2 
       mpow = spot_exponent_bound_polytope(pow);  
    else
       mpow = full(mpow); 
    end
       
    numHyperPlanes = 2000;
    [mpow] = RandomPrune(pow,mpow,numHyperPlanes);
    
    [mpow,coeffIndxPos,coeffIndxZero] = DiagConsistent(pow,mpow);

end
    
  
    
%Remove monomials m if 2*m is outside convhull(pow) using
%2*numHyperPlanes random seperating hyperplanes
function mpow = RandomPrune(pow,mpow,numHyperPlanes)

    M = size(pow,2);
    
    for i=1:numHyperPlanes
        
        w = randn(M,1);
        thres1 = max(pow*w);
        thres2 = min(pow*w);
 
        y = 2*mpow*w;
        
        %keep mpow if 2*mpow is in convhull(pow)
        mpow = mpow(y <= thres1 & y >= thres2,:);
               
        %Quit if we've pruned mpow enough.
        if (size(mpow,1) < 100)
            break;
        end
        
    end
        
end

%Get the "diagonally" consistent monomials. 2*mpow must be in pow
%or 2*mpow must be cancelled by a cross term
function [mpow,coeffIndxPos,coeffIndxZero] = DiagConsistent(pow,mpow)

    powUnSort = pow;

    %sorting and casting seems to speed up intersect()
    dataType = 'int16';
    pow = sortrows(cast(pow,dataType));
    mpow = sortrows(cast(mpow,dataType));
   
    %Pure square terms in f(x)
    [~,indxSqr] = intersect(2*mpow,pow,'rows');
    mpowSqr = mpow(indxSqr,:);

    N = size(mpow,1);
    while (1)

        %Compute cross terms 
        crossTerms = CrossTerms(mpow);

        %Find square terms that can be cancelled by cross terms
        if ~isempty(crossTerms)
            [~,indxCanc] = intersect(2*mpow,crossTerms,'rows');
        else
            indxCanc = [];
        end

        %Keep square terms if they are in pow or can be cancelled
        mpow = [mpowSqr;mpow(indxCanc,:)];
        mpow = unique(mpow,'rows');

        if (size(mpow,1) == N);
            break;
        else
            N = size(mpow,1);
        end

    end

    
    %These coeffs must be non-negative.
    if ~isempty(crossTerms)
        [~,coeffIndxPos] = setdiff(powUnSort,crossTerms,'rows');
    else
        coeffIndxPos = 1:size(pow,1);
    end

    %Coefficients that must vanish.
    if ~isempty(mpow)
        [~,coeffIndxZero] = setdiff(powUnSort,MinkSum(mpow),'rows');
    else
        coeffIndxZero = 1:size(pow,1);
    end

    mpow = double(mpow);
      
end

%Compute A+A, excluding term A_i + A_j if i = j
function [S] = CrossTerms(A)


    numTerm = size(A,1);
    numVar = size(A,2);
    
    %Get index position of upper triangular part, minus diagonal
    posKeep = find(triu(ones(numTerm,numTerm,'uint8'),1));

    if ~isempty(posKeep)
        %Compute A(i,:)+A(j,:) for all i \ne j
        S=zeros(length(posKeep),numVar,class(A));

        for k=1:numVar
           temp=bsxfun(@plus,A(:,k),A(:,k)');
           S(:,k) = temp(posKeep);
        end
        
    else
        S = [];
    end

    %Remove duplicates
    %[S]=unique(S,'rows');

end


%Computes A + A
function [S] = MinkSum(A)

    numTerm = size(A,1);
    numVar = size(A,2);
    
    %Get index position of upper triangular part
    posKeep = find(triu(ones(numTerm,numTerm,'uint8'),0));

    if ~isempty(posKeep)
        %Compute A(i,:)+A(j,:) for all i, j \in 1..numTerm
        S=zeros(length(posKeep),numVar,class(A));

        for k=1:numVar
           temp=bsxfun(@plus,A(:,k),A(:,k)');
           S(:,k) = temp(posKeep);
        end
        
    else
        S = [];
    end

    %Remove duplicates
    [S]=unique(S,'rows');

end







