function value = LZnEn(signal, isNorm)
%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculation of Generalized n-order LZ Entropy
%   Version [24/01/30] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal      : signal [1 x N]
%   isNorm      : (optional) normalization flag. Off: 0 / On: 1(default)
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   value       : n-order LZ Entropy
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] A.Lempel & J.Ziv, IEEE Transactions on Information Theory 22(1), 75
%       -81 (1976)
%   [2] Q.Thai, [MATLAB] calc_lz_complexity (2019)
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal = [85 80 89 85 80 89 85 80 89 85 80 89];
%   value1 = LZEn(signal, 0); % LZn complexity
%   value2 = LZEn(signal); % LZn entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    isNorm = 1;
end

signal = signal(:)';
N = length(signal);

eigIdxMat = zeros(N, 2);

vocaIdx = 1;
eigIdxMat(vocaIdx, 1) = 1;

eigCnt = 0;
for i = 1:N
    S = signal(1:eigCnt);
    Q = signal(eigCnt+1:i);
    SQ = [S, Q];
    SQpi = SQ(1:end-1);

    tempMat = (1:length(SQpi)-length(Q)+1)' + (0:length(Q)-1);
    SQpiList = SQpi(tempMat);
    SQpiList = reshape(SQpiList, size(tempMat));
    matches = all(bsxfun(@eq, SQpiList, Q), 2);

    if i~=N && ~any(matches)
        eigIdxMat(vocaIdx, 2) = i;
        vocaIdx = vocaIdx + 1;
        eigIdxMat(vocaIdx, 1) = eigIdxMat(vocaIdx-1, 2) + 1;
        eigCnt = i;
    end

end
eigIdxMat(vocaIdx, 2) = i;

eigIdxMat = nonzeros(eigIdxMat);
eigIdxMat = reshape(eigIdxMat, [], 2);

voca = arrayfun(@(x) signal(eigIdxMat(x, 1):eigIdxMat(x, 2)), ...
    1:size(eigIdxMat, 1), 'UniformOutput', false)';
value = length(voca);

if isNorm == 1
    base = length(unique(signal));
    if base == 1
        value = 0;
    else
        maxComplexity = N / (log(N)/log(base));
        value = value/maxComplexity;
    end
end

end