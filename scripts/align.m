%% open atlas
s = load('template/state.mat');

figure(3);
h = trisurf(s.fHead, s.vHead(:, 1), s.vHead(:, 2), s.vHead(:, 3));
set(h, 'LineStyle', 'none');
set(h, 'FaceAlpha', 0.7);
xlabel('x'); ylabel('y'); zlabel('z');
axis image;

hold on;
plot3(s.sideLeft3DVertices(:, 1), s.sideLeft3DVertices(:, 2), s.sideLeft3DVertices(:, 3), 'r', 'LineWidth', 3);
plot3(s.sideRight3DVertices(:, 1), s.sideRight3DVertices(:, 2), s.sideRight3DVertices(:, 3), 'r', 'LineWidth', 3);

for i = 1:size(s.refpts.pos, 1)
    plot3(s.refpts.pos(i, 1), s.refpts.pos(i, 2), s.refpts.pos(i, 3), 'k.');
    text(s.refpts.pos(i, 1), s.refpts.pos(i, 2), s.refpts.pos(i, 3), s.refpts.labels{i});
end

hold off;

%% align top
alignPanel('top', {'CPz' 'P1' 'Pz' 'P2' 'FOz' 'O1' 'Oz' 'O2' 'I1' 'Iz' 'I2' 'FP1' 'FPz' 'FP2' 'AF3' 'AFz' 'AF4' 'F1' 'Fz' 'F2' 'FC1' 'FCz' 'FC2' 'Cz' 'C1' 'C2' 'CP1' 'CP2' 'PO3' 'PO4'});

%% align side left
alignPanel('sideLeft', {'AF7' 'F3' 'F5' 'F7' 'F9' 'FC3' 'FC5' 'FT7' 'FT9' 'C3' 'C5' 'T7' 'LPA' 'CP3' 'CP5' 'TP7' 'TP9' 'P3' 'P5' 'P7' 'P9' 'PO3' 'PO7' 'PO9'});

%% align side right
alignPanel('sideRight', {'AF8' 'F4' 'F6' 'F8' 'F10' 'FC4' 'FC6' 'FT8' 'FT10' 'C4' 'C6' 'T8' 'RPA' 'CP4' 'CP6' 'TP8' 'TP10' 'P4' 'P6' 'P8' 'P10' 'PO4' 'PO8' 'PO10'});
