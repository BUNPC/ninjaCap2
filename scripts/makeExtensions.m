function extensions = makeExtensions(poly, pts, amount)
%MAKEEXTENSIONS Summary of this function goes here
%   Detailed explanation goes here

vertices = poly.Vertices;

extensions = nan(size(pts, 1) * 3, 2);
for i = 1:size(pts, 1)
    [boundary, v] = polygonNearestBoundary(poly, pts(i, :));
    
    % swap points if needed
    if norm(v - vertices(boundary(1), :), 2) < eps
        boundary = boundary([2 1]);
    end
    
    if poly.isinterior(lineMoveOrthogonal(v, vertices(boundary(1), :), eps))
        p = lineMoveOrthogonal(v, vertices(boundary(1), :), 0 - amount);
    else
        p = lineMoveOrthogonal(v, vertices(boundary(1), :), amount);
    end
    
    out_index = (i - 1) * 3;
    extensions(out_index + 1, :) = v;
    extensions(out_index + 2, :) = p;
end

end

