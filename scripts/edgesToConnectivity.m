function [vertices, connections, connection_matrix] = edgesToConnectivity(edges, tolerance)
%EDGESTOCONNECTIVITY Summary of this function goes here
%   Detailed explanation goes here

% tolerance
if ~exist('tolerance', 'var') || isempty(tolerance)
    tolerance = 1e-5;
end

% convert slice to tree
vertices = zeros(0, 2);
connections = zeros(0, 2);
for i = 1:3:size(edges, 1)
    [d1, v1] = min(sqrt(sum(bsxfun(@minus, edges(i, :), vertices) .^ 2, 2)));
    [d2, v2] = min(sqrt(sum(bsxfun(@minus, edges(i + 1, :), vertices) .^ 2, 2)));
    
    % did not find vertices? add them
    if isempty(d1) || d1 > tolerance
        vertices = [vertices; edges(i, :)];
        v1 = size(vertices, 1);
    end
    if isempty(d2) || d2 > tolerance
        vertices = [vertices; edges(i + 1, :)];
        v2 = size(vertices, 1);
    end
    
    connections = [connections; v1 v2];
end

% convert to connectivity matrix
if nargout >= 3
    connection_matrix = false(size(vertices, 1), size(vertices, 1));
    connection_matrix(sub2ind(size(connection_matrix), connections(:, 1), connections(:, 2))) = true;
    connection_matrix(sub2ind(size(connection_matrix), connections(:, 2), connections(:, 1))) = true; % symmetric
end

end

