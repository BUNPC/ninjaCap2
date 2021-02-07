function transform = alignImage(panel, im)
%ALIGNIMAGE Summary of this function goes here
%   Detailed explanation goes here

transform_type = 'projective'; % 'similarity' (3 pts) or 'projective' (4 pts)

% list of points
points_im = [];
points_panel = [];

while true
    % request points
    [cur_point_panel, cur_point_im] = alignImagePoint(panel, im);
    
    % either empty? stop
    if isempty(cur_point_im) || isempty(cur_point_panel)
        break;
    end
    
    % append to lists
    points_im = [points_im; cur_point_im]; %#ok<AGROW>
    points_panel = [points_panel; cur_point_panel]; %#ok<AGROW>
end

% enough points?
if strcmp(transform_type, 'projective')
    if size(points_im, 1) < 4
        error('Need at least 4 points for pojective transform.');
    end
elseif strcmp(transform_type, 'similarity')
    if size(points_im, 1) < 3
        error('Need at least 3 points for similarity transform.');
    end
else
    error('Unknown trasnsform type.');
end

% do transform
transform = fitgeotrans(points_im, points_panel, transform_type);

end
