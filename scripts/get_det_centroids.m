function centroids = get_det_centroids(detector_pos, ideal_edge_length, search_radius )

    num_detectors = size(detector_pos, 1);
    used_detectors = false(num_detectors, 1); % Tracker for used detectors
    triplets = []; % To store the results
    centroids = []; % To store the centroids of the found triplets
    triplet_count = 0;
%     search_radius = 15;
%     ideal_edge_length = 12;
    
    %%
    for i = 1:num_detectors
        % If this detector is already in a triplet, skip it
        if used_detectors(i)
            continue;
        end

        % current detector
        p1 = detector_pos(i, :);

        % available detectors
        available_indices = find(~used_detectors);
        available_indices = setdiff(available_indices, i); 
        available_pos = detector_pos(available_indices, :);

        % Calculate distances from p1 to all other available detectors
        distances = vecnorm(available_pos - p1, 2, 2);

        % Get the indices of neighbors within the radius
        neighbor_mask = distances <= search_radius;
        neighbor_indices = available_indices(neighbor_mask);

        % We need at least 2 neighbors to form a triangle
        if length(neighbor_indices) < 2
            continue;
        end

        % --- Search for the best pair among neighbors to form a triangle ---
        best_pair_indices = [];
        min_score = inf; % Lower score is better

        % Iterate through all combinations of two detectors from the neighbors
        combos = nchoosek(neighbor_indices, 2);

        for k = 1:size(combos, 1)
            j_idx = combos(k, 1);
            k_idx = combos(k, 2);

            % The other two vertices of the candidate triangle
            p2 = detector_pos(j_idx, :);
            p3 = detector_pos(k_idx, :);

            % --- Calculate the score based on ideal edge length ---
            % The score is the sum of squared differences from the ideal edge length.
            % This penalizes triangles that are not equilateral and not of the desired size.
            % A perfect score is 0.
            side_lengths = [norm(p1-p2), norm(p2-p3), norm(p3-p1)];

            % If any of the three side lengths is greater than the search radius,
            % this triangle is invalid, so we skip to the next combination.
            if any(side_lengths > search_radius)
                continue;
            end

            score = sum((side_lengths - ideal_edge_length).^2);

            % If this triangle is better than the best one found so far, update
            if score < min_score
                min_score = score;
                best_pair_indices = [j_idx, k_idx];
            end
        end

         % --- Finalize the Triplet ---
        % If a valid pair was found for detector 'i'
        if ~isempty(best_pair_indices)
            triplet_count = triplet_count + 1;

            final_triplet_indices = [i, best_pair_indices(1), best_pair_indices(2)];

            % Store the results
            triplets(triplet_count).indices = final_triplet_indices;
            triplets(triplet_count).score = min_score;

            % Calculate and store the centroid
            triplet_points = detector_pos(final_triplet_indices, :);
            centroids(triplet_count, :) = mean(triplet_points, 1);

            % Mark these three detectors as used
            used_detectors(final_triplet_indices) = true;
        end
    end