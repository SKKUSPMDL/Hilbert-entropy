function hCoord = gHilbertCoord3D(width, height, depth, x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz)

%%% Description %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Generalized Hilbert curve for 3D.
%   Version [24/02/10] SPMDL
%
%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   width       : Dimension.1 length of the Hilbert curve
%   height      : Dimension.2 length of the Hilbert curve
%   depth       : Dimension.3 length of the Hilbert curve
%   (x,y,z)     : Internal use only for recursive calls
%   (ax,ay,az)  : Internal use only for recursive calls
%   (bx,by,bz)  : Internal use only for recursive calls
%   (cx,cy,cz)  : Internal use only for recursive calls
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
%   depth = 6;
%   hCoord = gHilbertCoord3D(width,height,depth);
%   Hx = hCoord(:,1);
%   Hy = hCoord(:,2);
%   Hz = hCoord(:,3);
%   figure();
%   hold on;
%   lineColor = 1:width*height*depth;
%   surf([Hx Hx], [Hy Hy], [Hz Hz], [lineColor(:) lineColor(:)], ...
%       'FaceColor', 'none', ...
%       'EdgeColor', 'interp', ...
%       'LineWidth', 1);
%
%   [2]
%   width = 8;
%   height = 12;
%   depth = 6;
%   hCoord = gHilbertCoord3D(width,height,depth);
%   linIdx = sub2ind([width,height,depth], hCoord(:,1), hCoord(:,2), hCoord(:,3));
%   hIdx = zeros(width,height,depth);
%   hIdx(linIdx) = 1:size(hCoord, 1);
%
%   [3]
%   imgMat = randi(100,[8,12,6]);
%   width = size(imgMat,1);
%   height = size(imgMat,2);
%   depth = size(imgMat,3);
%   hCoord = gHilbertCoord3D(width,height,depth);
%   linIdx = sub2ind([width,height,depth], hCoord(:,1), hCoord(:,2), hCoord(:,3));
%   imgVec = imgMat(linIdx);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    if width >= height && width >= depth
        [x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz] = deal(0, 0, 0, width, 0, 0, 0, height, 0, 0, 0, depth);
    elseif height >= width && height >= depth
        [x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz] = deal(0, 0, 0, 0, height, 0, width, 0, 0, 0, 0, depth);
    else % depth >= width && depth >= height
        [x, y, z, ax, ay, az, bx, by, bz, cx, cy, cz] = deal(0, 0, 0, 0, 0, depth, width, 0, 0, 0, height, 0);
    end
end

w = abs(ax + ay + az);
h = abs(bx + by + bz);
d = abs(cx + cy + cz);

if h == 1 && d == 1
    hCoord = zeros(w, 3);
    for i = 0:w-1
        hCoord(i + 1, :) = [x + i * sign(ax), y + i * sign(ay), z + i * sign(az)];
    end
elseif w == 1 && d == 1
    hCoord = zeros(h, 3);
    for i = 0:h-1
        hCoord(i + 1, :) = [x + i * sign(bx), y + i * sign(by), z + i * sign(bz)];
    end
elseif w == 1 && h == 1
    hCoord = zeros(d, 3);
    for i = 0:d-1
        hCoord(i + 1, :) = [x + i * sign(cx), y + i * sign(cy), z + i * sign(cz)];
    end
else

    ax2 = floor(ax / 2); ay2 = floor(ay / 2); az2 = floor(az / 2);
    bx2 = floor(bx / 2); by2 = floor(by / 2); bz2 = floor(bz / 2);
    cx2 = floor(cx / 2); cy2 = floor(cy / 2); cz2 = floor(cz / 2);

    w2 = abs(ax2 + ay2 + az2);
    h2 = abs(bx2 + by2 + bz2);
    d2 = abs(cx2 + cy2 + cz2);

    if mod(w2, 2) && w > 2
        ax2 = ax2 + sign(ax);
        ay2 = ay2 + sign(ay);
        az2 = az2 + sign(az);
    end
    if mod(h2, 2) && h > 2
        bx2 = bx2 + sign(bx);
        by2 = by2 + sign(by);
        bz2 = bz2 + sign(bz);
    end
    if mod(d2, 2) && d > 2
        cx2 = cx2 + sign(cx);
        cy2 = cy2 + sign(cy);
        cz2 = cz2 + sign(cz);
    end

    if 2 * w > 3 * h && 2 * w > 3 * d
        points1 = gHilbertCoord3D(width, height, depth, x, y, z, ax2, ay2, az2, bx, by, bz, cx, cy, cz);
        points2 = gHilbertCoord3D(width, height, depth, x + ax2, y + ay2, z + az2, ax - ax2, ay - ay2, az - az2, bx, by, bz, cx, cy, cz);
        hCoord = [points1; points2];
    elseif 3 * h > 4 * d
        points1 = gHilbertCoord3D(width, height, depth, x, y, z, bx2, by2, bz2, cx, cy, cz, ax2, ay2, az2);
        points2 = gHilbertCoord3D(width, height, depth, x + bx2, y + by2, z + bz2, ax, ay, az, bx - bx2, by - by2, bz - bz2, cx, cy, cz);
        points3 = gHilbertCoord3D(width, height, depth, x + ax - sign(ax) + bx2 - sign(bx), y + ay - sign(ay) + by2 - sign(by), z + az - sign(az) + bz2 - sign(bz), -bx2, -by2, -bz2, cx, cy, cz, -(ax - ax2), -(ay - ay2), -(az - az2));
        hCoord = [points1; points2; points3];
    elseif 3 * d > 4 * h
        points1 = gHilbertCoord3D(width, height, depth, x, y, z, cx2, cy2, cz2, ax2, ay2, az2, bx, by, bz);
        points2 = gHilbertCoord3D(width, height, depth, x + cx2, y + cy2, z + cz2, ax, ay, az, bx, by, bz, cx - cx2, cy - cy2, cz - cz2);
        points3 = gHilbertCoord3D(width, height, depth, x + ax - sign(ax) + cx2 - sign(cx), y + ay - sign(ay) + cy2 - sign(cy), z + az - sign(az) + cz2 - sign(cz), -cx2, -cy2, -cz2, -(ax - ax2), -(ay - ay2), -(az - az2), bx, by, bz);
        hCoord = [points1; points2; points3];
    else
        points1 = gHilbertCoord3D(width, height, depth, x, y, z, bx2, by2, bz2, cx2, cy2, cz2, ax2, ay2, az2);
        points2 = gHilbertCoord3D(width, height, depth, x + bx2, y + by2, z + bz2, cx, cy, cz, ax2, ay2, az2, bx - bx2, by - by2, bz - bz2);
        points3 = gHilbertCoord3D(width, height, depth, x + bx2 - sign(bx) + cx - sign(cx), y + by2 - sign(by) + cy - sign(cy), z + bz2 - sign(bz) + cz - sign(cz), ax, ay, az, -bx2, -by2, -bz2, -(cx - cx2), -(cy - cy2), -(cz - cz2));
        points4 = gHilbertCoord3D(width, height, depth, x + ax - sign(ax) + bx2 + cx - sign(cx), y + ay - sign(ay) + by2 + cy - sign(cy), z + az - sign(az) + bz2 + cz - sign(cz), -cx, -cy, -cz, -(ax - ax2), -(ay - ay2), -(az - az2), bx - bx2, by - by2, bz - bz2);
        points5 = gHilbertCoord3D(width, height, depth, x + ax - sign(ax) + bx2 - sign(bx), y + ay - sign(ay) + by2 - sign(by), z + az - sign(az) + bz2 - sign(bz), -bx2, -by2, -bz2, cx2, cy2, cz2, -(ax - ax2), -(ay - ay2), -(az - az2));
        hCoord = [points1; points2; points3; points4; points5];
    end
end

if size(hCoord,1) == width*height*depth
    hCoord = hCoord + 1;
end
end