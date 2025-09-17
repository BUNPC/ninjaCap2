function [V, E] = connectionsToVE(connections)
    valid_coords = connections(~any(isnan(connections), 2), :);
    if isempty(valid_coords)
        V = [];
        E = [];
        return;
    end
    [V, ~, ic] = unique(valid_coords, 'rows', 'stable');
    E = reshape(ic, 2, [])';
end