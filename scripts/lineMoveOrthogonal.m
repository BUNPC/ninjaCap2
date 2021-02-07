function p = lineMoveOrthogonal(p0, p1, distance, towards)
%LINEMOVEORTHOGONAL Summary of this function goes here
%   Move orthogonally from line at p0 by distance. If towards is specified,
%   will move in whichever direction gets closer to towards (if distance is
%   positive) or whichever direction gets further from towards (if distance
%   is negative).

% get line normal
d = p1 - p0;
line_normal = d([2 1]) .* [-1 1];
line_normal = line_normal ./ norm(line_normal);


% move the lines along the line normal, ensure that it is moving outward
p = p0 + distance .* line_normal;

% if towards is specified, figure out which direction gets us closer or
% further
if exist('towards', 'var') && ~isempty(towards)
    p_alt = p0 - distance .* line_normal;
    
    % measure which point is closer
    rel_distance = norm(p_alt - towards) - norm(p - towards);
    
    % compare relative distance to sign of the distance
    if sign(rel_distance) ~= sign(distance)
        p = p_alt;
    end
end


end

