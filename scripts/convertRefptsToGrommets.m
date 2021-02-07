function grommets = convertRefptsToGrommets(refpts, subset)
%CONVERTREFPTSTOGROMMETS Summary of this function goes here
%   Detailed explanation goes here

grommets = struct('type', {}, 'flags', {}, 'posHead', {}, 'panel', {}, 'posPanel', {});

for i = 1:size(refpts.pos, 1)
    % allow selecting a subset
    if exist('subset', 'var') && ~isempty(subset) && ~ismember(refpts.labels(i), subset)
        continue;
    end
    
    grommet = struct('type', 'source', 'flags', {refpts.labels(i)}, 'posHead', refpts.pos(i, :), ...
        'panel', '', 'panelIndex', [], 'posPanel', [nan nan]);
    grommets = [grommets grommet]; %#ok<AGROW> % append
end

end
