function pts = sliceGetEdge(slice, p0, p1, nearest)
%SLICEGETEDGE Summary of this function goes here
%   Detailed explanation goes here

% default nearest
if ~exist('nearest', 'var')
    nearest = [];
end

% convert to connectivity
[vertices, ~, connectivity] = edgesToConnectivity(slice);

% find closest point
[~, i0] = closestVertex(vertices, p0);
[v1, i1] = closestVertex(vertices, p1);

% visited
visited = false(1, size(vertices, 1));

% return
pts = expandVertex(i0, visited);
if isempty(pts)
    error('No path found.');
end


    function path = expandVertex(i, visited)
        % mark as visited
        visited(i) = true;
        
        % is at final point?
        if i == i1
            path = v1;
            return;
        end
        
        % find connections (that have not been visited
        path = [];
        score = nan;
        other_is = find(connectivity(i, :) & ~visited);
        for other_i = other_is
            potential_path = expandVertex(other_i, visited);
            if isempty(potential_path)
                continue;
            end
            
            if isempty(nearest) || length(other_is) == 1
                path = [vertices(i, :); potential_path];
                break;
            else
                cur_score = min(sum(bsxfun(@minus, potential_path, nearest) .^ 2, 2));
                if isnan(score) || cur_score < score
                    score = cur_score;
                    path = [vertices(i, :); potential_path];
                end
            end
        end
    end

end

