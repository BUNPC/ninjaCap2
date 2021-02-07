function [idx, dist] = getVertNearLine(vertices, line, threshold)
%GETVERTNEARLINE Summary of this function goes here
%   Detailed explanation goes here

% default threshold
if ~exist('threshold', 'var')
    threshold = 0.1;
end

% line poly
poly = polybuffer(line, 'lines', threshold); % ,'JointType','miter','MiterLimit',lim);

% get vertices in line
idx = isinterior(poly, vertices);

% optionally: calculates the distance along the line (assume line(1, :) is
% the origin)
if nargout < 2
    dist = [];
else
    % nan for non intersecting points
    dist = nan(size(idx));
    
    % create poly shape for each line to figure out where it intersects
    lines_num = size(line, 1) - 1;
    lines_poly = polyshape();
    lines_dist = zeros(lines_num, 1);
    for i = 1:lines_num
        lines_poly(i) = polybuffer(line(i:(i + 1), :), 'lines', threshold);
        lines_dist(i) = sqrt(sum((line(i + 1, :) - line(i, :)) .^ 2));
    end
    
    % iterate over vertices
    for i = find(idx)'
        d = 0;
        for j = 1:lines_num
            if isinterior(lines_poly(j), vertices(i, :))
                d = d + sqrt(sum((vertices(i, :) - line(j, :)) .^ 2));
                break;
            end
            
            d = d + lines_dist(j);
        end
        dist(i) = d;
    end
end

end

