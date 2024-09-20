function hCoord = gHilbertCoord2D(width, height, x, y, ax, ay, bx, by)

%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generalized Hilbert curve for 2D.
%   Version [24/02/10] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   width       : Dimension.1 length of the Hilbert curve
%   height      : Dimension.2 length of the Hilbert curve
%   (x,y)       : Internal use only for recursive calls
%   (ax,ay)     : Internal use only for recursive calls
%   (bx,by)     : Internal use only for recursive calls
%
%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   hCoord      : A matrix of coordinates representing the Hilbert curve
%
%%% Reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1] J. Červený, [GitHUB] gilbert (2018)
%
%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   [1]
%   width = 8;
%   height = 12;
%   hCoord = gHilbertCoord2D(width,height);
%   Hx = hCoord(:,1);
%   Hy = hCoord(:,2);
%   Hz = zeros(width*height,1);
%   figure();
%   hold on;
%   lineColor = 1:width*height;
%   surf([Hx Hx], [Hy Hy], [Hz Hz], [lineColor(:) lineColor(:)], ...
%       'FaceColor', 'none', ...
%       'EdgeColor', 'interp', ...
%       'LineWidth', 1);
%
%   [2]
%   width = 8;
%   height = 12;
%   hCoord = gHilbertCoord2D(width,height);
%   linIdx = sub2ind([height,width], height+1-hCoord(:,2), hCoord(:,1));
%   hIdx = zeros(height,width);
%   hIdx(linIdx) = 1:size(hCoord, 1);
%
%   [3]
%   imgMat = randi(100,[8,12]);
%   width = size(imgMat,1);
%   height = size(imgMat,2);
%   hCoord = gHilbertCoord2D(width,height);
%   linIdx = sub2ind([height,width], height+1-hCoord(:,2), hCoord(:,1));
%   imgVec = imgMat(linIdx);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    if width >= height
        [x, y, ax, ay, bx, by] = deal(0, 0, width, 0, 0, height);
    else
        [x, y, ax, ay, bx, by] = deal(0, 0, 0, height, width, 0);
    end
end

w = abs(ax + ay);
h = abs(bx + by);

if h == 1
    hCoord = zeros(w, 2);
    for i = 0:w-1
        hCoord(i + 1, :) = [x + i * sign(ax), y + i * sign(ay)];
    end
elseif w == 1
    hCoord = zeros(h, 2);
    for i = 0:h-1
        hCoord(i + 1, :) = [x + i * sign(bx), y + i * sign(by)];
    end
else

    ax2 = floor(ax / 2); ay2 = floor(ay / 2);
    bx2 = floor(bx / 2); by2 = floor(by / 2);

    w2 = abs(ax2 + ay2);
    h2 = abs(bx2 + by2);

    if mod(w2, 2) == 1 && w > 2
        ax2 = ax2 + sign(ax);
        ay2 = ay2 + sign(ay);
    end
    if mod(h2, 2) == 1 && h > 2
        bx2 = bx2 + sign(bx);
        by2 = by2 + sign(by);
    end

    if 2 * w > 3 * h
        points1 = gHilbertCoord2D(width, height, x, y, ax2, ay2, bx, by);
        points2 = gHilbertCoord2D(width, height, x + ax2, y + ay2, ax - ax2, ay - ay2, bx, by);
        hCoord = [points1; points2];
    else
        points1 = gHilbertCoord2D(width, height, x, y, bx2, by2, ax2, ay2);
        points2 = gHilbertCoord2D(width, height, x + bx2, y + by2, ax, ay, bx - bx2, by - by2);
        points3 = gHilbertCoord2D(width, height, x + ax - sign(ax) + bx2 - sign(bx2), y + ay - sign(ay) + by2 - sign(by2), -bx2, -by2, -(ax - ax2), -(ay - ay2));
        hCoord = [points1; points2; points3];
    end
end

if size(hCoord,1) == width*height
    hCoord = hCoord + 1;
end
end
