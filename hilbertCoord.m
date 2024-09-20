function hCoord = hilbertCoord(order, dim, orgOrder)

%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generalized Hilbert curve for 2D and 3D.
%   Version [24/01/20] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   order       : Order of the Hilbert curve
%   dim         : Dimension of the space (2 or 3)
%   (orgOrder)  : Internal use only for recursive calls
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   hCoord      : A matrix of coordinates representing the Hilbert curve
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] F.Forte, [MATLAB] hilbert (2000)
%   [2] I.Martynov, [MATLAB] hilbert3 (2009)
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1]
%   order = 8;
%   dim = 2;
%   hCoord = hilbertCoord(order,dim);
%   Hx = hCoord(:,1);
%   Hy = hCoord(:,2);
%   Hz = zeros(2^(dim*order),1);
%   % Hz = hCoord(:,3); 
%   figure();
%   hold on;
%   lineColor = 1:2^(dim*order);
%   surf([Hx Hx], [Hy Hy], [Hz Hz], [lineColor(:) lineColor(:)], ...
%       'FaceColor', 'none', ...
%       'EdgeColor', 'interp', ...
%       'LineWidth', 1);
%
%   [2]
%   order = 2;
%   dim = 2;
%   hCoord = hilbertCoord(order,dim);
%   linIdx = sub2ind([2^order, 2^order], 2^order+1-hCoord(:,2), hCoord(:,1));
%   hIdx = zeros(2^order, 2^order);
%   hIdx(linIdx) = 1:size(hCoord, 1);
%
%   [3]
%   imgMat = randi(100,[8,8]);      % [=] 2^order x 2^order
%   order = log2(size(imgMat,1));
%   dim = ndims(imgMat);
%   hCoord = hilbertCoord(order, dim);
%   linIdx = sub2ind([2^order, 2^order], 2^order+1-hCoord(:,2), hCoord(:,1));
%   % linIdx = sub2ind([2^order, 2^order, 2^order], hCoord(:,1), hCoord(:,2), hCoord(:,3));
%   imgVec = imgMat(linIdx);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    orgOrder = order;
end

if order <= 0
    hCoord = zeros(1, dim);
    return;
end

prevHCoord = hilbertCoord(order-1, dim, orgOrder);

switch dim
    case 2
        x0 = prevHCoord(:,1);
        y0 = prevHCoord(:,2);
        x = .5*[-.5+y0; -.5+x0; +.5+x0; +.5-y0];
        y = .5*[-.5+x0; +.5+y0; +.5+y0; -.5-x0];
        hCoord = [x, y];
    case 3
        x0 = prevHCoord(:,1);
        y0 = prevHCoord(:,2);
        z0 = prevHCoord(:,3);
        x = .5*[-.5+z0; -.5+y0; +.5+y0; +.5-x0; +.5-x0; +.5-y0; -.5-y0; -.5+z0];
        y = .5*[-.5+x0; -.5+z0; -.5+z0; -.5+y0; +.5+y0; +.5-z0; +.5-z0; +.5-x0];
        z = .5*[-.5+y0; +.5+x0; +.5+x0; -.5-z0; -.5-z0; +.5+x0; +.5+x0; -.5-y0];
        hCoord = [x, y, z];
end

if order == orgOrder
    hCoord = (hCoord + 0.5) * 2^orgOrder + 0.5;
end

end
