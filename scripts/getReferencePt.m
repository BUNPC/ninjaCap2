function pt = getReferencePt(refpts, name)
%GETREFERENCEPT Summary of this function goes here
%   Detailed explanation goes here

idx = find(strcmp(refpts.labels, name), 1);
pt = refpts.pos(idx, :);

end
