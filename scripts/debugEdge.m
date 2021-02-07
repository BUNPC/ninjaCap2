function debugEdge(poly, vertices, dist)
%DEBUGEDGE Summary of this function goes here
%   Detailed explanation goes here

hold on;
plot(poly);
clr = mat2gray(dist);
scatter(vertices(:, 1), vertices(:, 2), 36, clr, '*');
colormap jet;
hold off;

end

