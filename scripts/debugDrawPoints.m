function debugDrawPoints(pos, labels, filled)
%DEBUGDRAWPOINTS Summary of this function goes here
%   Detailed explanation goes here

if ~exist('filled', 'var') || isempty(filled)
    filled = false;
end
if ~exist('labels', 'var')
    labels = {};
end

% draw points
if filled
    scatter(pos(:, 1), pos(:, 2), 'filled');
else
    scatter(pos(:, 1), pos(:, 2));
end

% draw labels
for i = 1:length(labels)
    if filled
        text(pos(i, 1), pos(i, 2), labels(i), 'VerticalAlignment', 'bottom');
    else
        text(pos(i, 1), pos(i, 2), labels(i), 'VerticalAlignment', 'top');
    end
end

end

