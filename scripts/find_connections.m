function connections = find_connections(centers, positions, k, max_distance)
    connections = []; 
    if isempty(centers) || isempty(positions), return; end
    for i = 1:size(centers, 1)
        center = centers(i, :);
        distances = pdist2(center, positions);
        [sorted_distances, sorted_indices] = sort(distances);
        closest_count = 0;
        for j = 1:length(sorted_distances)
            if sorted_distances(j) <= max_distance
                if closest_count < k
                    position_index = sorted_indices(j);
                    connections = [connections ; center; positions(position_index, :); [NaN NaN]];
                    closest_count = closest_count + 1;
                else
                    break;
                end
            else
                break;
            end
        end
    end
end