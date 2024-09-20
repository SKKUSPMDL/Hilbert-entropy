function value = permEn(signal, m, t, type, isNorm)
%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculation of Permutation Entropy
%   Version [24/01/30] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal      : signal [1 x N]
%   m           : length of sequences (m < N)
%   t           : delay
%   type        : "unique" or "tied"
%   isNorm      : (optional) normalization flag. Off: 0 / On: 1(default)
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   value       : Permutation Entropy
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] C.Bandt & B.Pompe, Physical Review Letter 88(17), 174102 (2002)
%   [2] V.A.Unakafova & K.Keller, Entropy 15(10), 4392-4415 (2013)
%   [3] M.Mor & A.S.Fraenkel, Discrete Mathematics, 48, 101-112 (1984)
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal = [85 80 89 85 80 89 85 80 89 85 80 89];
%   m = 3;
%   t = 1;
%   type = "tied";
%   value1 = permEn(signal, m, t, type, 0); % Permutation Entropy
%   value2 = permEn(signal, m, t, type); % normalized Permutation Entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    isNorm = 1;
end

signal = signal(:)';
N = length(signal);

sequenceNum = N-(1+(m-1)*t)+1;
tempMat = (1:sequenceNum)' + (0:t:(m-1)*t);
sequenceMat = signal(tempMat);
sequenceMat = reshape(sequenceMat, size(tempMat));

switch type
    case "unique"
        [~, sortIdx] = sort(sequenceMat,2);
        [~, permMat] = sort(sortIdx,2);
    case "tied"
        permMat = zeros(size(sequenceMat));
        for j = 1:size(sequenceMat,1)
            permMat(j,:) = tiedrank(sequenceMat(j,:));
        end
end
[~, ~, ic] = unique(permMat, 'rows');
counts = groupcounts(ic);
prob = counts/sequenceNum;
value = -sum(prob.*(log(prob)));

if isNorm == 1
    switch type
        case "unique"
            maxPermNum = factorial(m);
        case "tied"
            cPerm = zeros(1, m+1);
            cPerm(1) = 1; % f(0) = 1
            for k = 1:m
                temp = 0;
                for i = 1:k
                    temp = temp + nchoosek(k, i) * cPerm(k+1-i);
                end
                cPerm(k+1) = temp;
            end
            maxPermNum = cPerm(end);
    end
    maxProb = ones(maxPermNum,1)/maxPermNum;
    maxEntropy = -sum(maxProb.*(log(maxProb)));
    value = value/maxEntropy;
end

end