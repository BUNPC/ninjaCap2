function alignPanelView(panel, im, transform, points, labels)
%ALIGNPANELVIEW Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    panel_name = panel;

    % load image
    im = imread(sprintf('template/%s.jpg', panel_name));
    
    % load template
    load(sprintf('template/transform-%s.mat', panel_name), 'transform');
    
    % load mapping
    load(sprintf('template/map-%s.mat', panel_name), 'panel', 'points', 'labels');
end

%% show transform
h = figure(1);
ax = axes;
hold on;
plot(ax, panel); axis equal;
debugDrawPoints(points, labels);
hold off;
drawnow;

%% show image
h = figure(2); set(0, 'CurrentFigure', h);
ax = axes; set(h, 'CurrentAxes', ax);
imshow(im, 'Parent', ax); axis equal;
drawnow;

hold on;
v = union(panel(:));
v = v.Vertices;
v = transformPointsInverse(transform, v);
plot(ax, v(:, 1), v(:, 2), 'LineWidth', 2);
debugDrawPoints(transformPointsInverse(transform, points), labels);
hold off;

end
