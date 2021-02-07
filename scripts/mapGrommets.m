function [grommets,  Equidistnace_Outline_2D, Equidistnace_Outline_3D, refGrommetPositions, refGrommetLabels] = mapGrommets(grommets, vHead, fHead, refpts, panelName, outline, mapPoints, outline3DVertices, debug)
%MAPGROMMETS Summary of this function goes here
%   Detailed explanation goes here

if ~exist('debug', 'var')
    debug = false;
end

% buffer (if a grommet is within this distance, place on this panel
buffer = 10;%10%5; % mm

% scaling factor
scaling = 1.09;

% load mapping
mapping = load(sprintf('template/map-%s-10-5.mat', panelName));
refGrommetLabels = mapping.labels;
refGrommetPositions = mapping.points;
refOutline = mapping.outline;

% transform outline to match
if mapPoints > 0
    tform_cpp = fitgeotrans(refOutline(1:mapPoints, :), outline(1:mapPoints, :), 'nonreflectivesimilarity');
elseif size(refOutline, 1) == size(outline, 1)
    tform_cpp = fitgeotrans(refOutline, outline, 'similarity');
else
    error('Outlines do not have matching size, so they can not be mapped.');
end
refOutline = transformPointsForward(tform_cpp, refOutline);
refGrommetPositions = transformPointsForward(tform_cpp, refGrommetPositions);

% only use reference points that appear in both
if size(refpts.pos, 1) < 50
    warning('Grommet placement may be inaccurate. Recommend using 10-10 as reference points.');
end
idx = ismember(refGrommetLabels, refpts.labels);
refGrommetLabels = refGrommetLabels(idx);
refGrommetPositions = refGrommetPositions(idx, :);

% eeg labels and positions
idx = ismember(refpts.labels, refGrommetLabels);
eegLabels = refpts.labels(idx);
eegPositions = refpts.pos(idx, :);

% make panel
poly = polyshape(outline);
if buffer > 0
    poly = polybuffer(poly, buffer);
end

% map each grommet
for j = 1:length(grommets)
    % already placed
%     if ~isempty(grommets(j).panel)
%         continue;
%     end
    
    % place grommet
    [posPanel, mse, surroundingPts, Equidistnace_Outline_2D,  surroundingPts3D] = mapGrommet(grommets(j).posHead);
%     mse
%     isinterior(poly, posPanel)
%     if isinterior(poly, posPanel) && mse < 50
%         grommets(j).posPanel = posPanel;
%         grommets(j).panel = panelName;
%     end
%     if mse < 50
        grommets(j).posPanel = posPanel;
        grommets(j).panel = panelName;
        grommets(j).surroundingPts = surroundingPts;
        grommets(j).surroundingPts3D = surroundingPts3D;
%     end
end

% draw
if debug
    figure;
    plot(outline(:, 1), outline(:, 2), refOutline(:, 1), refOutline(:, 2));
    hold on;
    debugDrawPoints(refGrommetPositions, refGrommetLabels);
    idx_draw = arrayfun(@(x) ~isempty(x.panel) && strcmp(x.panel, panelName), grommets);
    if any(idx_draw)
        debugDrawPoints(cat(1, grommets(idx_draw).posPanel), cat(1, grommets(idx_draw).flags), true);
    end
    hold off;
    axis equal;
end

    function position = getRefPosition(label)
        for i = 1:length(refGrommetLabels)
            if strcmpi(refGrommetLabels{i}, label)
                position = refGrommetPositions(i, :);
                return;
            end
        end
        error('Unable to find position: %s.', label);
    end

    function positions = getRefPositions(labels)
        positions = zeros(length(labels), 2);
        for i = 1:length(labels)
            positions(i, :) = getRefPosition(labels{i});
        end
    end

    function [coord, mse, positions, Equidistnace_Outline_2D, positions_3D] = mapGrommet(grommet)
        
        % Use surrounding EIGHT 10-5 points to determine grommet position
        positions = [];
        distances = [];
        positions_3D = [];
        covered_idx = [];
        ind = find(isnan(outline3DVertices(:,1)) == 0);
        [Equidistnace_Outline_3D,~] = getEquiDistPtsAlongOutline(outline3DVertices(ind,:),1);
        if strcmp(panelName,'top')
            [Equidistnace_Outline_2D,~] = getEquiDistPtsAlongOutline(outline,1);
        else
            [Equidistnace_Outline_2D,~] = getEquiDistPtsAlongOutline(outline(ind,:),1);
        end
        Equidistnace_Outline_3D = Equidistnace_Outline_3D(1:20:end,:);
        Equidistnace_Outline_2D = Equidistnace_Outline_2D(1:20:end,:);
        %%
        eeg_and_outline_points = [eegPositions; Equidistnace_Outline_3D]; 
        % determine directions to use
        % [d1, d2] = find_directions(grommet, Equidistnace_Outline_3D, Equidistnace_Outline_3D);
        %%
        
        nV = size(vHead,1);
        h = 20;
        % position of the grommet
        p = grommet;

        % find nearest vertex on the head
        lst = find( abs(vHead(:,1)-p(1)*ones(nV,1))<h & abs(vHead(:,2)-p(2)*ones(nV,1))<h & abs(vHead(:,3)-p(3)*ones(nV,1))<h );
        rho = sum( (ones(length(lst),1)*p - vHead(lst,:)).^2, 2 ).^0.5;
        [foo,idx] = min(rho);

        iV = lst(idx);

        % get surface normal
        [ir,ic]=find(fHead==iV);

        iF = ir(1);
        if ic(1)==1
            p1 = vHead(fHead(iF,2),:);
            p2 = vHead(fHead(iF,3),:);
        elseif ic(1)==2
            p1 = vHead(fHead(iF,1),:);
            p2 = vHead(fHead(iF,3),:);
        elseif ic(1)==3
            p1 = vHead(fHead(iF,1),:);
            p2 = vHead(fHead(iF,2),:);
        end

        n = cross( p1-p, p2-p );
        n = n / norm(n);

        % get tangential axes
        a1 = p1-p; 
        a1 = a1 / norm(a1);
        a2 = cross( a1, n );
        a2 = a2 / norm(a2);
%         lstSub = find( abs(eeg_and_outline_points(:,1)-p(1)*ones(nV,1))<h & abs(eeg_and_outline_points(:,2)-p(2)*ones(nV,1))<h & abs(eeg_and_outline_points(:,3)-p(3)*ones(nV,1))<h );
        %%
        eeg_length = size(eegPositions,1);
        for uu = 1:4
            if uu == 1
                search_idx = find((eeg_and_outline_points-p)*a1' > 0 );
%                 search_idx = find(eeg_and_outline_points(:,2) <= grommet(:,2));
                search_points = eeg_and_outline_points(search_idx,:);
            end
            if uu == 2
                search_idx = find((eeg_and_outline_points-p)*a1' <= 0 );
%                 search_idx = find(eeg_and_outline_points(:,2) > grommet(:,2));
                search_points = eeg_and_outline_points(search_idx,:);
            end
            if uu == 3
                search_idx = find((eeg_and_outline_points-p)*a2' > 0 );
%                 search_idx = find(eeg_and_outline_points(:,3) <= grommet(:,3));
                search_points = eeg_and_outline_points(search_idx,:);
            end
            if uu == 4
                search_idx = find((eeg_and_outline_points-p)*a2' <= 0 );
%                 search_idx = find(eeg_and_outline_points(:,3) > grommet(:,3));
                search_points = eeg_and_outline_points(search_idx,:);
            end
            num = size(search_points, 1);
            temp_distances = zeros(num, 1);
            for i = 1:num
                temp_distances(i) = surfaceDist(search_points(i,:), grommet);
            end
            N = 2;
            [~, idx1] = sort(temp_distances);
%             idx = idx(1:8);
            idx2 = search_idx(idx1);
            idx3 = setdiff(idx2, covered_idx,'stable');
            if length(idx3) >= 2
                idx3 = idx3(1:N);
            end
            covered_idx = [covered_idx; idx3];
            idx4 = idx3(idx3 <= eeg_length);
            idx5 = idx3(idx3 > eeg_length);
            positions_3D = [positions_3D; eeg_and_outline_points(idx3,:)];
%             distances = [distances; temp_distances(idx)];
            if ~isempty(idx4)
                positions = [positions; getRefPositions(eegLabels(idx4))];
                idx6 = find(ismember(idx2,idx4)==1);
                distances = [distances; temp_distances(idx1(idx6))];
            end
            if ~isempty(idx5)
                 positions = [positions; Equidistnace_Outline_2D(idx5-eeg_length,:)];
                 idx6 = find(ismember(idx2,idx5)==1);
                 distances = [distances; temp_distances(idx1(idx6))];
            end
        end
        [distances, d_idx] = sort(distances);
        distances = distances/scaling;
%         distances = distances/1.19;
        positions = positions(d_idx,:);
%%
%         % Using nearest SIX 10-5 points to determine grommet position
%         num = size(eegPositions, 1);
%         distances = zeros(num, 1);
%         for i = 1:num
%             distances(i) = surfaceDist(eegPositions(i, :), grommet);
%         end
%         % use N closest points
%         N=6;
%         [~, idx] = sort(distances);
%         idx = idx(1:N);
%         labels = eegLabels(idx);
%         distances = distances(idx) ./ scaling;
%         
%         % get reference positions
%         positions = getRefPositions(labels);
%         
%         ind = find(isnan(outline3DVertices(:,1)) == 0);
%         [Equidistnace_Outline_3D,~] = getEquiDistPtsAlongOutline(outline3DVertices(ind,:),1);
%         if strcmp(panelName,'top')
%             [Equidistnace_Outline_2D,~] = getEquiDistPtsAlongOutline(outline,1);
%         else
%             [Equidistnace_Outline_2D,~] = getEquiDistPtsAlongOutline(outline(ind,:),1);
%         end
%         % IN ADDITION: use distance from closest outline point, if closer
%         num2 = size(Equidistnace_Outline_3D, 1);
%         distancesOutline = ones(num2, 1)*1e9;
%         for i = 1:15:num2
%             distancesOutline(i) = surfaceDist(Equidistnace_Outline_3D(i, :), grommet);
%         end
%         distancesOutline = distancesOutline ./ scaling;
% 
%         if ~isempty(distancesOutline)
% %             [closest, idx] = min(distancesOutline);
%             [~, idx] = sort(distancesOutline);
%             idx = idx(1:N);
%             distancesOutline = distancesOutline(idx);
%             outlinePositions = Equidistnace_Outline_2D(idx, :);
%             distances = [distances; distancesOutline];
%             positions = [positions; outlinePositions];
%             [~, idx] = sort(distances);
%             idx = idx(1:N);
%             positions = positions(idx,:);
%             distances = distances(idx,:);
% %                 if closest < distances(end)
% %                     positions = [positions; outline(idx, :)];
% %                     distances = [distances; closest];
% %                 end
%         end
%%
        % start at the closest point
        coord0 = positions(1, :);
        
        options = optimset('MaxFunEvals', 50000);
        % find best position
        
        [coord, mse] = fminsearch(@(x) outlineDistance(x, positions, distances), coord0, options );
        
        % debugging: print out fit information
        %disp(mse);
        %disp([distances sqrt(sum(bsxfun(@minus, coord, outline) .^ 2, 2))]);
        
        % debugging: plot
%         figure;
%         plot(outline(:, 1), outline(:, 2), 'k');
%         hold on;
%         scatter(positions(:, 1), positions(:, 2));
%         viscircles(positions, distances);
%         scatter(coord(1), coord(2), 'filled');
%         hold off;
        
        % is closest point more than 25, inflate the error
        if distances(1) > 50
            mse = inf;
        end
    end

    function err = outlineDistance(coord, positions, desired_distances)
        distances = sqrt(sum(bsxfun(@minus, coord, positions) .^ 2, 2));
        err = mean((distances - desired_distances) .^ 2);
    end

    function dist = surfaceDist(vertex, point)
        % OPTION 1: straight euclidian distance
        %dist = sqrt(sum((vertex - point) .^ 2));
        
        % OPTION 2: sphere distance (assume a sphere, and get the arc
        % length (slightly better than euclidean distance)
        dist = sphereDistance(vHead, vertex, point);
        
        % OPTION 3: surface distance
        % model does not work well for this (ears and other things can
        % create artifacts)
        %dist = surfaceDistance(vHead, fHead, vertex, point);
        %if isinf(dist)
        %    dist = nan;
        %end
    end

%     function [d1, d2] = find_directions(grommet, Equidistnace_Outline_3D, Equidistnace_Outline_3D)
%         positions = [];
%         distances = [];
%         covered_idx = [];
%         eeg_and_outline_points = [eegPositions; Equidistnace_Outline_3D]; 
%          for uu = 1:6
%             if uu == 1
%                 search_idx = find(eeg_and_outline_points(:,1) <= grommet(:,1));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%             if uu == 2
%                 search_idx = find(eeg_and_outline_points(:,1) > grommet(:,1));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%             if uu == 3
%                 search_idx = find(eeg_and_outline_points(:,2) <= grommet(:,2));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%             if uu == 4
%                 search_idx = find(eeg_and_outline_points(:,2) > grommet(:,2));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%              if uu == 5
%                 search_idx = find(eeg_and_outline_points(:,3) <= grommet(:,3));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%             if uu == 6
%                 search_idx = find(eeg_and_outline_points(:,3) > grommet(:,3));
%                 search_points = eeg_and_outline_points(search_idx,:);
%             end
%             num = size(search_points, 1);
%             temp_distances = zeros(num, 1);
%             for i = 1:num
%                 temp_distances(i) = surfaceDist(search_points(i,:), grommet);
%             end
%             N = 2;
%             [~, idx1] = sort(temp_distances);
% %             idx = idx(1:8);
%             idx2 = search_idx(idx1);
%             idx3 = setdiff(idx2, covered_idx,'stable');
%             if length(idx3) >= 2
%                 idx3 = idx3(1:N);
%             end
%             covered_idx = [covered_idx; idx3]; 
%             positions = [positions; eeg_and_outline_points(idx3,:)];
%          end
%          diff_vectors = positions - grommet;
% %          unit_vectors = diff_vectors./(sum(diff_vectors.^2,2))
%     end

end
