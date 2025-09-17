function eOut_oriented = orientBoundaryEdges(vOut, eOut, polygon_coords)
%orientBoundaryEdges Reorders the columns of a boundary edge list.
%
%   This function ensures that for each edge in eOut, the vertex that is
%   INSIDE the specified polygon is listed in the first column, and the
%   vertex that is OUTSIDE is listed in the second column.
%
%   INPUTS:
%   vOut           - Nx2 matrix of [x, y] boundary vertex coordinates.
%   eOut           - Mx2 matrix of boundary edges, with indices corresponding to vOut.
%   polygon_coords - Px2 matrix of [x, y] coordinates defining the polygon.
%
%   OUTPUTS:
%   eOut_oriented  - The Mx2 reordered ("oriented") edge list.

% If there are no edges, there's nothing to do.
if isempty(eOut)
    eOut_oriented = [];
    return;
end

%% Step 1: Determine which vertices in vOut are inside the polygon
% 'is_inside' will be a logical vector (true for inside, false for outside).
is_inside = inpolygon(vOut(:,1), vOut(:,2), polygon_coords(:,1), polygon_coords(:,2));

%% Step 2: Loop through each edge and reorder if necessary
eOut_oriented = eOut; % Initialize with the original to handle ambiguous cases

for i = 1:size(eOut, 1)
    % Get the indices of the two vertices for the current edge
    v1_idx = eOut(i, 1);
    v2_idx = eOut(i, 2);
    
    % Check if the first vertex is outside and the second is inside
    if ~is_inside(v1_idx) && is_inside(v2_idx)
        % The edge is ordered [outside, inside], so we must swap it.
        eOut_oriented(i, :) = [v2_idx, v1_idx];
    end
    % Note: If the edge is already [inside, outside], we do nothing.
    % If both are inside or both are outside (an ambiguous case for a
    % boundary edge), we also do nothing and keep the original order.
end

end