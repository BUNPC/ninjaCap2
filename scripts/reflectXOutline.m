function [outline, varargout] = reflectXOutline(outline, varargin)
%REFLECTXOUTLINE Summary of this function goes here
%   Detailed explanation goes here

assert(nargin == nargout);

% find centroid
mn = min(outline, [], 1);
mx = max(outline, [], 1);
x = [mn(1) mx(1)];
y = [mn(2) mx(2)];
c = [mean(x) mean(y)];
    
% move to origin
outline = bsxfun(@minus, outline, c);

% reflect
outline = bsxfun(@times, outline, [-1 1]);

% move back to center
outline = bsxfun(@plus, outline, c);

% other arguments
for i = 1:(nargin - 1)
    disp(i)
    v = bsxfun(@minus, varargin{i}, c);
    v = bsxfun(@times, v, [-1 1]);
    v = bsxfun(@plus, v, c);
    varargout{i} = v;
end

end
