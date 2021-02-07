function distance = sphereDistance(vertices, p1, p2)
%SPHEREDISTANCE Summary of this function goes here
%   Detailed explanation goes here

% get mid point
p_mid = closestVertex(vertices, (p1 + p2) ./ 2);
%p1 = closestVertex(vertices, p1);
%p2 = closestVertex(vertices, p2);

% length and height of arc
l = sqrt(sum((p1 - p2) .^ 2));

% https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
if length(p1) == 3
    % project to 2d
    plane = makePlane([p1; p2; p_mid]);
    [~, pts] = projectPlane(plane, [p1; p2; p_mid]);
else
    pts = [p1; p2; p_mid];
end
h = abs((pts(2, 2) - pts(1, 2)) * pts(3, 1) - (pts(2, 1) - pts(1, 1)) * pts(3, 2) + pts(2,1) * pts(1, 2) - pts(2, 2) * pts(1, 1)) / sqrt(sum((pts(2, :) - pts(1, :)) .^ 2));

% from: https://en.wikipedia.org/wiki/Arc_(geometry)
% and: https://en.wikipedia.org/wiki/Circular_segment

% circle radius
r = l ^ 2 / (8 * h) + h / 2;

% arc angle
ang = 2 * asin(l / (2 * r));

% calculate distance
distance = ang * r;

% debug, print both
%disp([sqrt(sum((p1 - p2) .^ 2)) distance]);

end

