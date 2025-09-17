function [vOut_clean, eOut_clean] = cleanBoundaryData(vOut, eOut)
%cleanBoundaryData Removes orphan vertices from boundary data.
%   This function takes a set of boundary vertices (vOut) and edges (eOut)
%   and removes any vertex from vOut that is not referenced in eOut. It then
%   re-indexes the edge list to match the new, smaller vertex list.
%
%   INPUTS:
%   vOut - Nx2 matrix of boundary vertex coordinates.
%   eOut - Mx2 matrix of boundary edges, with indices corresponding to vOut.
%
%   OUTPUTS:
%   vOut_clean - Cleaned list of boundary vertex coordinates.
%   eOut_clean - Cleaned and re-indexed list of boundary edges.

% If there are no edges, there are no valid vertices to keep.
if isempty(eOut)
    vOut_clean = [];
    eOut_clean = [];
    return;
end

% Find all unique vertex indices that are actually used in the edge list
used_indices = unique(eOut(:));

% Create the final, clean vertex list by keeping only the used vertices
vOut_clean = vOut(used_indices, :);

% To re-index the edges, create a map from the old vOut indices to the new ones
% The map will have a size equal to the original number of vertices
map_old_to_new = zeros(size(vOut, 1), 1);

% Populate the map: the old index (e.g., used_indices(i)) now maps to the
% new index (e.g., i)
map_old_to_new(used_indices) = 1:length(used_indices);

% Apply the map to the original edge list to get the new, clean edge list
eOut_clean = map_old_to_new(eOut);

end