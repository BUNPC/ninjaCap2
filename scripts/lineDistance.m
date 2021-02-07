function dist = lineDistance(l, pt)
%LINEDISTANCE Summary of this function goes here
%   Detailed explanation goes here

% https://stackoverflow.com/a/1501725/383539

% get l2 norm
l2 = sum((l(1, :) - l(2, :)) .^ 2);

% no length? distance from either one
if l2 == 0
    dist = norm(pt - l(1, :), 2);
    return;
end

t = ((pt(1) - l(1, 1)) * (l(2, 1) - l(1, 1)) + (pt(2) - l(1, 2)) * (l(2, 2) - l(1, 2))) / l2;

% constrain to line segment (remove for distance from line)
t = max(0, min(t, 1));

dist = norm(pt - (l(1, :) + t * (l(2, :) - l(1, :))), 2);

end
