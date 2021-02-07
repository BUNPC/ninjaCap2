function p = lineIntersect(l1, l2)
%LINEINTERSECT Summary of this function goes here
%   l1 is a 2x2 matrix: [x1 y1; x2 y2]
%   l2 is a 2x2 matrix: [x3 y3; x4 y4]

% REFERENCE: http://mathworld.wolfram.com/Line-LineIntersection.html

x = det([det(l1) l1(1, 1) - l1(2, 1); ...
    det(l2) l2(1, 1) - l2(2, 1)]) / ...
    det([l1(1, :) - l1(2, :); l2(1, :) - l2(2, :)]);

y = det([det(l1) l1(1, 2) - l1(2, 2); ...
    det(l2) l2(1, 2) - l2(2, 2)]) / ...
    det([l1(1, :) - l1(2, :); l2(1, :) - l2(2, :)]);

p = [x y];

end

