function [v, f, refpts] = loadHead(cir)
%LOADHEAD Load AtlasViewer atlas for head model
%   This function loads a surface for a generic head from AtlasViewer. For
%   production caps, a saved viewer state (atlasViewer.mat) should be used 
%   instead.

% get atlas diredtory
atlas = getAtlasDir();

%% load surface frmo AtlasViewer
headsurf = initHeadsurf();
headsurf = getHeadsurf(headsurf, atlas);
v = headsurf.mesh.vertices;
f = headsurf.mesh.faces;

%% load referene points from AtlasViewer
refpts = initRefpts();
refpts = getRefpts(refpts, atlas);

%% scale if needed
% calculate scaling constant based on target circumference
if exist('cir', 'var') && ~isempty(cir)
    scale = cir / refpts.eeg_system.lengths.circumference;
    
    % apply scaling
    v = v * scale;
    refpts.pos = refpts.pos * scale;
end

%% plot head surface
% markers for reference points on head
[xs, ys, zs] = sphere(20);
[fs, vs] = surf2patch(xs, ys, zs,' triangles');

% reference point positions
pos = refpts.pos;

figure;
h = trisurf(f, v(:, 1), v(:, 2), v(:, 3));
alpha(0.5);
set(h, 'linestyle', 'none');
set(h, 'facealpha', 1);
light;
axis image;

hold on;
for ii = 1:size(pos, 1)
    % alternative style
    %hp = plot3(pos(ii, 1), pos(ii, 2), pos(ii, 3), 'k.');
    %ht = text(pos(ii, 1), pos(ii, 2) - 100, pos(ii, 3), refpts.labels{ii});
    %set(hp, 'markersize', 16);
    
    h2 = trisurf(fs, pos(ii, 1) + vs(:, 1), pos(ii,2) + vs(:, 2), pos(ii, 3) + vs(:, 3));
    alpha(0.5);
    set(h2, 'linestyle', 'none');
end
hold off;

end