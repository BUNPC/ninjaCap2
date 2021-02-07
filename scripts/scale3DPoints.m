function ptsScaled = scale3DPoints(vHead, standard, pts)
%SCALE3DPOINTS Summary of this function goes here

% make sure head mapping matches
if size(standard.vHead, 1) ~= size(vHead, 1)
    ptsScaled = [];
    return;
end

% extra
x = interp1(standard.vHead(:, 1), vHead(:, 1), pts(:, 1), 'linear', 'extrap');
y = interp1(standard.vHead(:, 2), vHead(:, 2), pts(:, 2), 'linear', 'extrap');
z = interp1(standard.vHead(:, 3), vHead(:, 3), pts(:, 3), 'linear', 'extrap');

% return new points
ptsScaled = [x y z];

end