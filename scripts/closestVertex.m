function [v, i, dist] = closestVertex(vertices, pt)
%CLOSESTVERTEX Summary of this function goes here
%   Detailed explanation goes here

[dist, i] = min(sum(bsxfun(@minus, pt, vertices) .^ 2, 2));
dist = sqrt(dist);
v = vertices(i, :);

end

