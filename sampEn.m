function value = sampEn(signal, m, r, type, isNorm)
%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculation of Sample Entropy
%   Version [24/01/30] SPMDL
%           [24/03/10] Optimized for Low-Resource Environments
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal      : signal [1 x N]
%   m           : length of sequences (m < N)
%   r           : tolerance
%   type        : "light"(L.res, H.eff) / "heavy"(H.res, L.eff; default)
%   isNorm      : (optional) normalization flag. Off: 0 / On: 1(default)
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   value       : Sample Entropy
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] J.S.Richman & J.R.Moorman, American Journal of Physiology - Heart a
%       nd Circulatory Physiology, 278(6), H2039-H2049 (2000)
%   [2] V.Martínez-Cagigal, [MATLAB] sampen (2018)
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal = [85 80 89 85 80 89 85 80 89 85 80 89];
%   m = 2;
%   r = 0.2;
%   type = "heavy";
%   value1 = sampEn(signal, m, r, type, 0); % Sample Entropy
%   value2 = sampEn(signal, m, r, type); % normalized Sample Entropy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5
    isNorm = 1;
end

if nargin < 4
    type = "heavy";
end

signal = signal(:)';
N = length(signal);
signal = double(signal);
sigma = std(signal);

maxEntropy = -log(2/((N-m-1)*(N-m)));

tempMat = (1:N-m)' + (0:m);
sequenceMat = signal(tempMat);
sequenceMat1 = reshape(sequenceMat, size(tempMat));
sequenceMat0 = sequenceMat1(:,1:end-1);

switch type
    case "heavy"
        dm0 = pdist(sequenceMat0, 'chebychev');
        dm1 = pdist(sequenceMat1, 'chebychev');
        numSim0 = sum(dm0 <= r*sigma); % B
        numSim1 = sum(dm1 <= r*sigma); % A
    case "light"
        dm0 = zeros(N-m, 1);
        dm1 = zeros(N-m, 1);
        for i = 1:N-m
            dm0(i) = sum(max(abs(sequenceMat0(i,:)-sequenceMat0),[],2) <= r*sigma);
            dm1(i) = sum(max(abs(sequenceMat1(i,:)-sequenceMat1),[],2) <= r*sigma);
        end
        numSim0 = sum(dm0) - (N-m); % B
        numSim1 = sum(dm1) - (N-m); % A
end

value = -log(numSim1/numSim0);
if numSim0*numSim1 == 0
    value = maxEntropy;
end

if isNorm == 1
    value = value/maxEntropy;
end

end