function alignPanel(panel_name, expected_labels)
%ALIGNPANEL Summary of this function goes here
%   Detailed explanation goes here

% default
if ~exist('expected_labels', 'var')
    expected_labels = {};
end

state = load('template/state.mat');

% get panel image
v = [panel_name 'Panel'];
panel = state.(v);

% get outline
v = [panel_name 'Outline'];
outline = state.(v);

% load image
im = imread(sprintf('template/%s.jpg', panel_name));

%% get transform
fl = sprintf('template/transform-%s.mat', panel_name);
if exist(fl, 'file')
    load(fl, 'transform');
else
    transform = alignImage(panel, im);
    save(fl, '-v7.3', 'transform');
end

%% get approximate points
fl = sprintf('template/map-%s.mat', panel_name);

if exist(fl, 'file')
    mapping = load(fl);
    approx_panel = mapping.panel;
    approx_labels = mapping.labels;
    approx_points = mapping.points;
else
    old_mapping = load('mapGrommets.mat');
    if strcmp(panel_name, 'top')
        v = 'topOutline';
    else
        v = 'sideOutline';
    end
    approx_panel = old_mapping.(v);
    v = [panel_name 'GrommetPositions'];
    approx_points = old_mapping.(v);
    v = [panel_name 'GrommetLabels'];
    approx_labels = old_mapping.(v);
end

% add any missing labels
for i = 1:length(expected_labels)
    if ~ismember(expected_labels{i}, approx_labels)
        approx_labels{end + 1} = expected_labels{i};
        approx_points = [approx_points; nan nan];
    end
end

[points, labels] = alignRefPts(panel, im, transform, approx_panel, approx_points, approx_labels);
save(fl, '-v7.3', 'panel', 'outline', 'points', 'labels');

%% show result
alignPanelView(panel, im, transform, points, labels);

end
