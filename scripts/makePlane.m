function plane = makePlane(pts)
%MAKEPLANE Summary of this function goes here
%   Detailed explanation goes here

% make deltas
d1 = pts(2, :) - pts(1, :);
d2 = pts(3, :) - pts(1, :);

% get normal
normal = cross(d1, d2);

% make basis
b1 = cross(normal, d2);
b1 = b1 / norm(b1);

b2 = cross(normal, b1);
b2 = b2 / norm(b2);

plane = [pts(1, :) b1 b2];

end

