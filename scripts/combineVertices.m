function [vertices, dist] = combineVertices(vertices, dist, threshold)
%COMBINEVERTICES Summary of this function goes here
%   Detailed explanation goes here

% make distance matrix
d = abs(dist' - dist);

% clear lower section
d(logical(tril(ones(size(d))))) = nan;

% clear super threshold
d(d > threshold) = nan;

% merge vertices
[merge_dist, merge_into] = min(d, [], 2);

% only merge once
moved = false(size(vertices, 1), 1);

% do merging
for i = find(~isnan(merge_dist))'
    % merge into
    j = merge_into(i);
    
    % only move each vertex once
    if any(moved([i j]))
        continue;
    end
    moved([i j]) = true;
    
    % get new vertex
    v = mean(vertices([i, j], :), 1);
    
    % assign new vertex
    vertices([i j], :) = [v; v];
    
    % get new distance
    d = mean(dist([i j]));
    
    % merge distances
    dist([i j]) = d;
end

end
