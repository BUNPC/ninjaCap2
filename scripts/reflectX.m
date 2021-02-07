function [poly, grommets] = reflectX(poly, grommets)
%REFLECTX Summary of this function goes here
%   Detailed explanation goes here

% find centroid
[x, y] = boundingbox(poly);
c = [mean(x) mean(y)];

for i = 1:numel(poly)
    v = poly(i).Vertices;
    
    % move to origin
    v = bsxfun(@minus, v, c);
    
    % reflect
    v = bsxfun(@times, v, [-1 1]);
    
    % move back to center
    v = bsxfun(@plus, v, c);
    
    poly(i).Vertices = v;
end

if exist('grommets', 'var')
    for i = 1:length(grommets)
        v = grommets(i).posPanel;
        v = c + ((v - c) .* [-1 1]);
        grommets(i).posPanel = v;
    end
end

end

