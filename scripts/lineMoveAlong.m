function p = lineMoveAlong(p0, p1, distance)
%LINEMOVEALONG Summary of this function goes here
%   Move along line p0 towards p1 by distance

% generate slope
slope = p1 - p0;
slope = slope ./ norm(slope);

p = p0 + slope .* distance;

end

