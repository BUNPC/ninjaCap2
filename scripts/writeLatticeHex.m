function vertices_clipped = writeLatticeHex(polyin, hEdge)
%WRITELATTICEHEX Generate hexagonal lattice vertices to fill polygon
%   Detailed explanation goes here

%% initial calculations

% get bounding box
[lim_poly_x, lim_poly_y] = boundingbox(polyin);

% vertices for one cell
vertices_cell = [...
    0 0; ...
    hEdge 0; ...
    nan nan; ...
    0 0; ...
    hEdge * cos(2 * pi / 3) hEdge * sin(2 * pi / 3); ...
    nan nan; ...
    0 0; ...
    hEdge * cos(2  *pi / 3) -hEdge * sin(2 * pi / 3); ...
    nan nan; ...
    ];

lim_cell_x = [min(vertices_cell(:, 1)) max(vertices_cell(:, 1))];
lim_cell_y = [min(vertices_cell(:, 2)) max(vertices_cell(:, 2))];

% center
center_x = sum(lim_poly_x) / 2;
center_y = sum(lim_poly_y) / 2;

% figure out lattice offset to produce nice edges (useful for top panel)
best_offset_x = 0;
best_offset_x_score = inf; % lower is better
for offset_x = linspace(lim_cell_x(1), lim_cell_x(2), 60)
    % figure out left and right ediges
    edge1 = mod(lim_poly_x(1) - (center_x + offset_x), lim_cell_x(2) - lim_cell_x(1));
    edge2 = mod(lim_poly_x(2) - (center_x + offset_x), lim_cell_x(2) - lim_cell_x(1));
    
%     % bad match
%     if edge1 > lim_cell_x(2) || edge2 > lim_cell_x(2)
%         continue;
%     end
    
    % edges that fall between 0 and 10 (lim_cell_x(2)) will intersect
    % horizontal range, so score based on square distance to ideal (5)
    score = (edge1 - lim_cell_x(2) / 2) ^ 2 + (edge2 - lim_cell_x(2) / 2) ^ 2;
    if score < best_offset_x_score
        best_offset_x = offset_x;
        best_offset_x_score = score;
    end
end

center_x = center_x + best_offset_x;

%% generate lattice

% multiply by two since per row, cell alternates (i.e., [0, 0], [2, 0], [4,
% 0], etc; alternate rows have alternates, [1, 1], [3, 1], [5, 1], etc)
repeat_every_x = lim_cell_x(2) - lim_cell_x(1);

% generate x multipliers
% figure out how many cells to span the space to the left of the center
% coordinate; and how many cells spane the space to the right of the center
% (use conservative estimates to ensure fully spanned, hence floor and
% ceil)
% also: add one extra cell to ensure cell fully spans bounding box; add
% another extra cell to the right side to ensure that the alternate
% (shifted) rows still fully span bounding box
left = -(1 + floor((center_x - lim_poly_x(1)) / (2 * repeat_every_x)));
right = (2 +  ceil((lim_poly_x(2) - center_x) / (2 * repeat_every_x)));

idx_x = 2 * (left:right);
pos_x = center_x + idx_x * repeat_every_x;

% shift alternate rows
pos_x_alt = center_x + (idx_x - 1) * repeat_every_x;

% similar logic for y
repeat_every_y = lim_cell_y(2) - lim_cell_y(1);

% generate y multipliers
% figure out how many cells to span the space below the center coordinate;
% and how many cells to span the space above the center coordinate (use
% conservative estimates to ensure fully spanned)
below = -(1 + floor((center_y - lim_poly_y(1)) / repeat_every_y));
above =  (2 +  ceil((lim_poly_y(2) - center_y) / repeat_every_y));

idx_y = below:above;
pos_y = center_y + idx_y * repeat_every_y;

% shift alternate rows
pos_y_alt = center_y + idx_y * repeat_every_y - (repeat_every_y / 2);

% generate all permutations of pos_x and pos_y
pos = [repelem(pos_x', length(pos_y)) repmat(pos_y', length(pos_x), 1)];
pos_alt = [repelem(pos_x_alt', length(pos_y_alt)) repmat(pos_y_alt', length(pos_x_alt), 1)];

% combine
pos = [pos; pos_alt];

%% make vertices

% number of cells
cells = size(pos, 1);

% make vertices by repeating pos (x/y offsets) for each cell and adding it
% to the cell coordinates (repeated for the number of cells)
vertices = repelem(pos, size(vertices_cell, 1), 1) + repmat(vertices_cell, cells, 1);

%% clip to polyin
%inside_point_1 = isinterior(polyin, vertices(1:3:end, :));
%inside_point_2 = isinterior(polyin, vertices(2:3:end, :));

vertices_clipped = [];
for i = 1:3:size(vertices, 1)
    [v_in, ~] = intersect(polyin, vertices(i:(i+1), :));
    if ~isempty(v_in)
        vertices_clipped = [vertices_clipped; v_in; nan nan]; %#ok<AGROW>
    end
end

end
