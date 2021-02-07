function polys = cut(poly, max_size, buffer)
%CUT Summary of this function goes here
%   Detailed explanation goes here

if ~exist('buffer', 'var')
    buffer = 0;
end

% get bounding box
[bb_x, bb_y] = boundingbox(poly);

% make return variable
polys = poly;

% fits?
if (bb_x(2) - bb_x(1)) <= max_size(1) && (bb_y(2) - bb_y(1)) <= max_size(2)
    return;
end

% remove buffers
max_size = max_size - 2 * buffer;

% number of pieces
pieces_x = ceil((bb_x(2) - bb_x(1)) / max_size(1));
pieces_y = ceil((bb_y(2) - bb_y(1)) / max_size(2));

% figure out origin
origin_x = mean(bb_x) - max_size(1) * pieces_x / 2;
origin_y = mean(bb_y) - max_size(2) * pieces_y / 2;

% for each piece, perform cuts
for i = 1:pieces_x
    for j = 1:pieces_y
        x1 = origin_x + max_size(1) * (i - 1) - buffer;
        x2 = x1 + max_size(1) + buffer * 2;
        
        y1 = origin_y + max_size(2) * (j - 1) - buffer;
        y2 = y1 + max_size(2) + buffer * 2;

        % create mask
        mask = polyshape([x1 y1; x2 y1; x2 y2; x1 y2]);
        
        % cut polygon
        polys((i - 1) * pieces_y + j, :) = intersect(mask, poly);
    end
end

end
