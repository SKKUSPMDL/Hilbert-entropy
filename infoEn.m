function value = infoEn(signal, isNorm)
%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculation of Shannon's Information Entropy
%   Version [24/02/14] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal      : digital signal [1 x N]
%   isNorm      : (optional) normalization flag. Off: 0 / On: 1(default)
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   value       : Information Entropy
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] C.E.Shannon, The Bell System Technical Journal 27(3), 379-423 (1948
%       )
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   signal = [0 1 0 0 1 0 1 1 1 1 0];
%   value1 = infoEn(signal, 0); % infromation entropy with base-e log
%   value2 = infoEn(signal); % infromation entropy with base-n log
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    isNorm = 1;
end

signal = signal(:)';
N = length(signal);
[~, ~, ic] = unique(signal);
cnt = groupcounts(ic);
prob = cnt/N;
value = -sum(prob.*(log(prob)));

if isNorm == 1 && value ~= 0
    maxProb = ones(size(prob))/length(cnt);
    maxEntropy = -sum(maxProb.*(log(maxProb)));
    value = value/maxEntropy;
end

end