function [pt_panel, pt_image] = alignImagePoint(panel, im)
%ALIGNIMAGEPOINT Summary of this function goes here
%   Detailed explanation goes here

% placeholders
pt_panel = [];
pt_image = [];

% open figure
h = figure(1);
h.Position = h.Position .* [1 1 2 1];

ax1 = subplot(1, 2, 1);
plot(panel); axis equal; hold on;
set(allchild(ax1), 'HitTest', false);
plt_panel = [];

ax2 = subplot(1, 2, 2);
imshow(im, 'Parent', ax2); axis equal; hold on;
set(ax2, 'HitTest', true, 'PickableParts', 'all'); % restore axis functionality
set(allchild(ax2), 'HitTest', false);
plt_image = [];

% listen for clicks
set(ax1, 'ButtonDownFcn', {@cb_clickPanel});
set(ax2, 'ButtonDownFcn', {@cb_clickImage});

% wait for window to close
waitfor(h);

    % handle clicks
    function cb_clickPanel(~, event)
        ax = event.Source;
        
        % get position
        pos = get(ax, 'CurrentPoint');
        pos = pos(1, 1:2); % convert from 3d to 2d
        
        % store it
        pt_panel = pos;
        
        % draw it
        if ~isempty(plt_panel)
            delete(plt_panel);
        end
        plt_panel = scatter(ax, pos(1), pos(2), 36, 'filled');
    end

    
    function cb_clickImage(~, event)
        ax = event.Source;
        
        % get position
        pos = get(ax, 'CurrentPoint');
        pos = pos(1, 1:2); % convert from 3d to 2d
        
        % store it
        pt_image = pos;
        
        % draw it
        if ~isempty(plt_image)
            delete(plt_image);
        end
        plt_image = scatter(ax, pos(1), pos(2), 36, 'filled');
    end
end

