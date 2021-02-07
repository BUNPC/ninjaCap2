function panel = createChinStrapSlot(panel, outline, idx, width, slot_dim)
%CREATECHINSTRAPSLOT Summary of this function goes here
%   Detailed explanation goes here

% idx should be in order (uses -1 and +1 to access adjacenet points, which
% requires idx to be in order)
idx = sort(idx);

% mid point
m = mean(outline(idx, :), 1);

% add outline on either side of slot
slot_width = slot_dim(1) + width * 2;

% OUTSIDE

% move away from outline by desired slot width
edge_outside = [lineMoveAlong(m, outline(idx(1), :), slot_width / 2); ...
    lineMoveAlong(m, outline(idx(2), :), slot_width / 2)];

% top edge will overlap with the existing outside outline, some move in by
% that amount
slot_outside = [...
    lineMoveOrthogonal(edge_outside(1, :), m, width, outline(idx(1) - 1, :)); ... % top: point 1
    lineMoveOrthogonal(edge_outside(2, :), m, width, outline(idx(2) + 1, :)); ... % top: point 2
    lineMoveOrthogonal(edge_outside(2, :), m, - slot_dim(2) - width, outline(idx(2) + 1, :)); ... % bottom: point 2
    lineMoveOrthogonal(edge_outside(1, :), m, - slot_dim(2) - width, outline(idx(1) - 1, :))]; ... % bottom point 1

% INSIDE

% move away from outline by desired slot width
edge_inside = [lineMoveAlong(m, outline(idx(1), :), slot_dim(1) / 2); ...
    lineMoveAlong(m, outline(idx(2), :), slot_dim(1) / 2)];

% top edge will overlap with the existing outside outline, some move in by
% that amount
slot_inside = [...
    edge_inside(1, :); ... % top: point 1
    edge_inside(2, :); ... % top: point 2
    lineMoveOrthogonal(edge_inside(2, :), m, - slot_dim(2), outline(idx(2) + 1, :)); ... % bottom: point 2
    lineMoveOrthogonal(edge_inside(1, :), m, - slot_dim(2), outline(idx(1) - 1, :))]; ... % bottom point 1

% create slot
slot = subtract(polyshape(slot_outside), polyshape(slot_inside));

% add slot to panel
panel(:, 1) = union(panel(:, 1), slot);

% % plot outline
% figure;
% plot(outline(:, 1), outline(:, 2));
% hold on;
% plot(slot_outside(:, 1), slot_outside(:, 2), 'rx');
% plot(slot_inside(:, 1), slot_inside(:, 2), 'rx');
% hold off;

% figure;
% plot(panel);

end

