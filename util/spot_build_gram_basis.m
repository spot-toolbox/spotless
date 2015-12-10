function [m_partition,coeffIndxPos,coeffIndxZero] = spot_build_gram_basis(pow,mpow)
  
    pow = full(pow);
    
    if nargin < 2 
       mpow = spot_exponent_bound_polytope(pow);  
    else
       mpow = full(mpow); 
    end
       
    numHyperPlanes = 2000;
    [mpow] = RandomPrune(pow,mpow,numHyperPlanes);  
    [mpow,coeffIndxPos,coeffIndxZero] = DiagConsistent(pow,mpow);
    [~,m_partition] = BlockDiagonalize(pow,mpow);
 
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
    
    %Compute A(i,:)+A(j,:) for all i, j \in 1..numTerm
    S = zeros(numTerm*numTerm,numVar,class(A)); %#ok<ZEROLIKE>

    for k=1:numVar
       temp=bsxfun(@plus,A(:,k),A(:,k)');
       S(:,k) = temp(:);
    end
        
end


function [pow,m_partition] = BlockDiagonalize(pow,mpow)

    n = size(mpow,1);
    [minkSum] = MinkSum(mpow);
  
    M = zeros(n);
    while (1)

        size_pow = size(pow,1);
        [indxKeep] = find(ismember(minkSum,pow,'rows'));
        M(indxKeep) = 1;
        [M,partition] = TransitiveClosure(M);
        pow = minkSum(M(:) == 1,:);
        if size(pow,1) == size_pow, break, end
        
    end
  
    for i = 1:length(partition)
        m_partition{i} = mpow( partition{i},:);
    end
  

end



function [M,cliques] = TransitiveClosure(M)

[r,~] = find(M);
if isempty(r)
    cliques = {}; return
end

r = unique(r);
[~,~,cliques] = conncomp(M(r,r));

for i=1:length(cliques)
    cliques{i} = r(cliques{i});
    M(cliques{i},cliques{i}) = 1;
end

end


function [nComponents,sizes,members] = conncomp(A)
% Number of nodes
N = size(A,1);
% Remove diagonals
A(1:N+1:end) = 0;
% make symmetric, just in case it isn't
A=A+A';
% Have we visited a particular node yet?
isDiscovered = zeros(N,1);
% Empty members cell
members = {};
% check every node
for n=1:N
    if ~isDiscovered(n)
        % started a new group so add it to members
        members{end+1} = n;
        % account for discovering n
        isDiscovered(n) = 1;
        % set the ptr to 1
        ptr = 1;
        while (ptr <= length(members{end}))
            % find neighbors
            nbrs = find(A(:,members{end}(ptr)));
            % here are the neighbors that are undiscovered
            newNbrs = nbrs(isDiscovered(nbrs)==0);
            % we can now mark them as discovered
            isDiscovered(newNbrs) = 1;
            % add them to member list
            members{end}(end+1:end+length(newNbrs)) = newNbrs;
            % increment ptr so we check the next member of this component
            ptr = ptr+1;
        end
    end
end
% number of components
nComponents = length(members);
for n=1:nComponents
    % compute sizes of components
    sizes(n) = length(members{n});
end

[sizes,idx] = sort(sizes,'ascend');
members = members(idx);



% Copyright (c) 2014, Daniel Larremore
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
end

