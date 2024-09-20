function value = LZEn(signal, isNorm)
%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculation of LZ Entropy
%   Version [24/01/30] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal      : digital signal [1 x N]
%   isNorm      : (optional) normalization flag. Off: 0 / On: 1(default)
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   value       : LZ Entropy
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
%   signal = (signal > mean(signal));
%   value1 = LZEn(signal, 0); % LZ complexity
%   value2 = LZEn(signal); % LZ entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    isNorm = 1;
end

signal = signal(:)';
N = length(signal);
signal = num2str(signal, '%g');

eigIdxMat = zeros(N, 2);

eigIdxMat(1, :) = 1;
endIdx = 1;
count = 1;

for i = 2:N
    S = signal(1:endIdx);
    Q = signal(endIdx+1:i);
    SQ = [S, Q];
    SQpi = SQ(1:end-1);

    if ~contains(SQpi, Q)
        endIdx = i;
        startIdx = eigIdxMat(count, 2) + 1;
        count = count + 1;
        eigIdxMat(count, 2) = endIdx;
        eigIdxMat(count, 1) = startIdx;
    end
end
if endIdx < N
    eigIdxMat(count+1, 2) = N;
    eigIdxMat(count+1, 1) = endIdx+1;
end

eigIdxMat = nonzeros(eigIdxMat);
eigIdxMat = reshape(eigIdxMat, [], 2);

voca = arrayfun(@(x) signal(eigIdxMat(x, 1):eigIdxMat(x, 2)), ...
    1:size(eigIdxMat, 1), 'UniformOutput', false)';
value = length(voca);

if isNorm == 1
    numAlphabet = unique(signal);
    if length(numAlphabet) == 1
        value = 0;
    else
        maxComplexity = N/(log2(N));
        value = value/maxComplexity;
    end
end

end