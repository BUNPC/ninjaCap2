function [idx, v, dist] = polygonNearestBoundary(poly, pt)
%POLYGONNEARESTBOUNDARY Summary of this function goes here
%   Detailed explanation goes here

% get points and rotated points
pts = poly.Vertices;
pts2 = [pts(2:end, :); pts(1, :)];
dif = pts2 - pts;

% figure out l2's
l2 = sum(dif .^ 2, 2);

% make t
t = (pt(1) - pts(:, 1)) .* (pts2(:, 1) - pts(:, 1)) + (pt(2) - pts(:, 2)) .* (pts2(:, 2) - pts(:, 2));
t(l2 == 0) = 0;
t(l2 ~= 0) = t ./ l2(l2 ~= 0);

% constrain
t = max(0, min(t, 1));

% calculate distance
compare = pts + t .* dif;

% get closest vertex
[v, i, dist] = closestVertex(compare, pt);

% return boundary indices
if i == size(pts, 1)
    idx = [i 1];
else
    idx = [i i + 1];
end

end

