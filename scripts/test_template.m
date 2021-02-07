% reset grommets
grommets = convertRefptsToGrommets(refpts);

% grommets
grommets = mapGrommets(grommets, vHead, fHead, refpts, 'top', topOutline, 0, true);
grommets = mapGrommets(grommets, vHead, fHead, refpts, 'sideLeft', sideLeftOutline, sideOutsidePoints, true);
grommets = mapGrommets(grommets, vHead, fHead, refpts, 'sideRight', sideRightOutline, sideOutsidePoints, true);

% load reference points
ref = load('mapGrommets.mat');

% draw: top
figure;
plot(topOutline(:, 1), topOutline(:, 2))
hold on;
idx = arrayfun(@(x) strcmp(x.panel, 'top'), grommets);
if any(idx)
    debugDrawPoints(cat(1, grommets(idx).posPanel), cat(1, grommets(idx).flags), true);
end
debugDrawPoints(ref.topGrommetPositions, ref.topGrommetLabels);
hold off;
axis equal;

% draw: left side
figure;
plot(sideLeftOutline(:, 1), sideLeftOutline(:, 2))
hold on;
idx = arrayfun(@(x) strcmp(x.panel, 'sideLeft'), grommets);
if any(idx)
    debugDrawPoints(cat(1, grommets(idx).posPanel), cat(1, grommets(idx).flags), true);
end
debugDrawPoints(ref.sideLeftGrommetPositions, ref.sideLeftGrommetLabels);
hold off;

% draw: right side
figure;
plot(sideRightOutline(:, 1), sideRightOutline(:, 2))
hold on;
idx = arrayfun(@(x) strcmp(x.panel, 'sideRight'), grommets);
if any(idx)
    debugDrawPoints(cat(1, grommets(idx).posPanel), cat(1, grommets(idx).flags), true);
end
debugDrawPoints(ref.sideRightGrommetPositions, ref.sideRightGrommetLabels);
hold off;
