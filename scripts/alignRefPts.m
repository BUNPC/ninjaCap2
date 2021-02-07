function [points, labels] = alignRefPts(panel, im, transform, approx_panel, approx_points, approx_labels)
%ALIGNREFPTS Summary of this function goes here
%   Detailed explanation goes here

% provide instructions
fprintf('For each reference point, click the location in the image.\n');
fprintf('To use the current point (shown in red), press ENTER.\n');
fprintf('If the point is not on the panel, press ESC.\n');

% transform panel to lines that can be super imposed
p = union(panel(:));
p = transformPointsInverse(transform, p.Vertices);

% close all figures
close all;

% draw figures
h1 = figure(1);
h1.Position = h1.Position .* [1 1 2 2];
ax1 = axes;
imshow(im, 'Parent', ax1); axis equal; 

% draw the panel
hold(ax1, 'on');
plot(p(:, 1), p(:, 2));

% restore axis functionality
set(ax1, 'HitTest', true, 'PickableParts', 'all');
set(allchild(ax1), 'HitTest', false);

h2 = figure(2);
ax2 = axes;

% draw the approximate points
hold(ax2, 'on');
if isnumeric(approx_panel)
    plot(approx_panel(:, 1), approx_panel(:, 2));
else
    plot(approx_panel);
end
axis equal;
debugDrawPoints(approx_points, approx_labels);

points = [];
labels = {};

% for each reference point
for i = 1:length(approx_labels)
    % print label
    fprintf('Locate %s...\n', approx_labels{i});
    title(ax1, sprintf('Click %s', approx_labels{i}));
    
    % show current point
    cur = transformPointsInverse(transform, approx_points(i, :));
    h_cur = plot(ax1, cur(1, 1), cur(1, 2), 'ro', 'MarkerSize', 24);
    
    % draw active point
    h_active = plot(ax2, approx_points(i, 1), approx_points(i, 2), 'r.', 'MarkerSize', 36);
    drawnow;
    
    while true
        % set active figure
        set(0, 'CurrentFigure', h1);
        set(h1, 'CurrentAxes', ax1);
        keydown = waitforbuttonpress;
        
        % figure closed
        if ~ishghandle(h1)
            error('Figure closed before completion.');
        end

        % check figure
        figchildren = allchild(0);
        if ~isempty(figchildren)
            ptr_fig = figchildren(1);
            if ptr_fig ~= h1
                continue;
            end
        else
            continue;
        end

        % keypress
        if keydown
            char = get(h1, 'CurrentCharacter');
            if char == 13 % enter
                % use existing point
                pt = cur;
            elseif abs(char) == 27 %escape
                break;
            else
                continue;
            end
        else
            button = get(h1, 'SelectionType');
            if strcmp(button,'open')
                button = 1;
            elseif strcmp(button,'normal')
                button = 1;
            elseif strcmp(button,'extend')
                button = 2;
            elseif strcmp(button,'alt')
                button = 3;
            else
                error('Unable to interpret input.');
            end
            
            % unexpected button
            if button ~= 1
                continue;
            end

            % get point
            pt = get(ax1, 'CurrentPoint');  
            pt = pt(1, [1 2]);
        end
        
        % add point and label
        plot(ax1, pt(1), pt(2), '.b', 'MarkerSize', 14);
        text(pt(1), pt(2), approx_labels{i}, 'VerticalAlignment', 'bottom');
        
        points = [points; pt]; %#ok<AGROW>
        labels{end + 1} = approx_labels{i}; %#ok<AGROW>
        break;
    end
    
    % remove active point
    delete(h_cur);
    delete(h_active);
end

% transform points
points = transformPointsForward(transform, points);

% close everything
close(h1);
close(h2);

end

