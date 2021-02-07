% This script was designed to take a measured head circumference, and
% produce a fitted headcap filled with a hexagonal pattern.

% Dependencies are in the private folder

%% Parameters

% whether or not to show useful debugging visualizations
debug = true;

%% Load head model
% Use the view state (atlasViewer.mat).
% fpath = 'D:\Office\Research\Software - Scripts\Matlab\ninjaCap\atlasViewer.mat';
% load the AtlasViewer state file, including a head model and grommets
[vHead, fHead, refpts, grommets] = loadAtlasViewer('atlasViewer.mat');

% no state file? you can alternatively use the following to produce a
% generic head, but this code requires AtlasViewer installed:
%[vHead, fHead, refpts] = loadHead();
%grommets = [];

%% Load side piece
% New function reads in outline from a SVG file
fname = 'svg/side_piece2.svg';
pathid = 'side_piece';
% side panel outline and number of side panel points that define the outside vertices
[side_piece, sideOutsidePoints] = svgImport(fname, pathid);

%% Segment the head
% Divides the head into three sections to create outlines for the top and
% side panels
[sideLeftOutline, sideRightOutline, topOutline, sideLeft3DVertices, sideRight3DVertices, top3DVertices] = segmentHead(vHead, fHead, refpts, side_piece);
sideChinStrapLine = [44 45]; % indices of the two points that define the chin strap line

%% Map grommet positions

if ~isempty(grommets)
    grommets = mapGrommets(grommets, vHead, fHead, refpts, 'top', topOutline, 0, top3DVertices, debug);
    grommets = mapGrommets(grommets, vHead, fHead, refpts, 'sideLeft', sideLeftOutline, sideOutsidePoints, sideLeft3DVertices, debug);
    grommets = mapGrommets(grommets, vHead, fHead, refpts, 'sideRight', sideRightOutline, sideOutsidePoints, sideRight3DVertices, debug);
end

%% Build basic panels

% create polygons
topPoly = polyshape(topOutline);
sideLeftPoly = polyshape(sideLeftOutline);
sideRightPoly = polyshape(sideRightOutline);

% hexagon edge for lattice
hEdge = 10;

% write the hexagonal pattern into the outlines
topLattice = writeLatticeHex(topPoly, hEdge);
sideLeftLattice = writeLatticeHex(sideLeftPoly, hEdge);
sideRightLattice = writeLatticeHex(sideRightPoly, hEdge);

% debug
if debug
    h = figure; h.Position = h.Position .* [1 1 1 3];
    subplot(3, 1, 1); debugLattice(sideLeftPoly, sideLeftLattice); axis equal;
    subplot(3, 1, 2); debugLattice(topPoly, topLattice); axis equal;
    subplot(3, 1, 3); debugLattice(sideRightPoly, sideRightLattice); axis equal;
end

%% Align panel connections

% SETTINGS
% if there are short lattice edges at the top of the side panel, just move
% the nearest vertex to the edge and eliminate the short edge (threshold in
% mm; empty to disable)
align_eliminate_short_edges = hEdge * 0.2;
% if there are nearby vertices on the top panel, combine them (threshold in
% mm; empty to disable)
align_combine_top_vertices = hEdge * 0.2;
% shift vertices at the top of the side panels to align with the top panel
% (true/false)
align_shift_side_vertices = true;
% if a vertex is moved substantially on the side panel, add another edge to
% reinforce it (threshold in mm, empty to disable)
align_reinforce_side_vertices = [];
% whether or not all edge vertices of side and top panels should be
% connected; recommended if no stitching outline (true/false)
align_fully_connect = true;

% get edges
edgeTopLeft = topOutline([1 2], :);
edgeTopRight = topOutline([4 3], :);
edgeSideLeft = sideLeftOutline((sideOutsidePoints+1):end, :);
edgeSideRight = sideRightOutline((sideOutsidePoints+1):end, :);

% eliminate short edges on side panels
if ~isempty(align_eliminate_short_edges)
    edgeSideLeftVerticeIdx = getVertNearLine(sideLeftLattice, edgeSideLeft);
    sideLeftLattice = eliminateShortEdges(sideLeftLattice, edgeSideLeftVerticeIdx, align_eliminate_short_edges);
    
    edgeSideRightVerticeIdx = getVertNearLine(sideRightLattice, edgeSideRight);
    sideRightLattice = eliminateShortEdges(sideRightLattice, edgeSideRightVerticeIdx, align_eliminate_short_edges);
end

% debug
if debug
    h = figure; h.Position = h.Position .* [1 1 1 3];
    subplot(3, 1, 1); debugLattice(sideLeftPoly, sideLeftLattice); axis equal;
    subplot(3, 1, 2); debugLattice(topPoly, topLattice); axis equal;
    subplot(3, 1, 3); debugLattice(sideRightPoly, sideRightLattice); axis equal;
end

% get edge vertices
[edgeTopLeftVerticeIdx, distTopLeftVertice] = getVertNearLine(topLattice, edgeTopLeft);
[edgeTopRightVerticeIdx, distTopRightVertice] = getVertNearLine(topLattice, edgeTopRight);
[edgeSideLeftVerticeIdx, distSideLeftVertice] = getVertNearLine(sideLeftLattice, edgeSideLeft);
[edgeSideRightVerticeIdx, distSideRightVertice] = getVertNearLine(sideRightLattice, edgeSideRight);

% plot showing edge vertices for side panel
if debug
    h = figure; h.Position = h.Position .* [1 1 1 3];
    subplot(3, 2, [1 2]); debugEdge(sideLeftPoly, sideLeftLattice(edgeSideLeftVerticeIdx, :), distSideLeftVertice(edgeSideLeftVerticeIdx)); axis equal;
    subplot(3, 2, 3); debugEdge(topPoly, topLattice(edgeTopLeftVerticeIdx, :), distTopLeftVertice(edgeTopLeftVerticeIdx)); axis equal;
    subplot(3, 2, 4); debugEdge(topPoly, topLattice(edgeTopRightVerticeIdx, :), distTopRightVertice(edgeTopRightVerticeIdx)); axis equal;
    subplot(3, 2, [5 6]); debugEdge(sideRightPoly, sideRightLattice(edgeSideRightVerticeIdx, :), distSideRightVertice(edgeSideRightVerticeIdx)); axis equal;
end

% combine vertices on top if they are close
if ~isempty(align_combine_top_vertices)
    [topLattice, distTopLeftVertice] = combineVertices(topLattice, distTopLeftVertice, align_combine_top_vertices);
    [topLattice, distTopRightVertice] = combineVertices(topLattice, distTopRightVertice, align_combine_top_vertices);
end

% adjust vertices on side panels
if align_shift_side_vertices
    sideLeftLattice = shiftVertices(sideLeftLattice, distSideLeftVertice, distTopLeftVertice(edgeTopLeftVerticeIdx), edgeSideLeft, align_reinforce_side_vertices, align_fully_connect);
    sideRightLattice = shiftVertices(sideRightLattice, distSideRightVertice, distTopRightVertice(edgeTopRightVerticeIdx), edgeSideRight, align_reinforce_side_vertices, align_fully_connect);
end

% debug
if debug
    h = figure; h.Position = h.Position .* [1 1 3 1];
    subplot(1, 5, [1 2]); debugLattice(sideLeftPoly, sideLeftLattice); axis equal;
    subplot(1, 5, 3); debugLattice(topPoly, topLattice); axis equal;
    subplot(1, 5, [4 5]); debugLattice(sideRightPoly, sideRightLattice); axis equal;
end

%% Extend edge lattice pieces
% uses edge vertices entries from last section

width_extensions = 2;

% extend edge vertices (adding overlap that may help in bonding pieces)
extend_edge_vertices = 5;

% outline extensions (need to move in by half the extension width, to line
% up correctly with the outline)

vertices = [...
    lineMoveAlong(edgeTopLeft(1, :), edgeTopLeft(2, :), width_extensions / 2); ...
    lineMoveAlong(edgeTopLeft(2, :), edgeTopLeft(1, :), width_extensions / 2); ...
    lineMoveAlong(edgeTopRight(1, :), edgeTopRight(2, :), width_extensions / 2); ...
    lineMoveAlong(edgeTopRight(2, :), edgeTopRight(1, :),  width_extensions / 2)];
topExtensionsOutline = makeExtensions(topPoly, vertices, extend_edge_vertices);

vertices = [...
    lineMoveAlong(edgeSideLeft(1, :), edgeSideLeft(2, :), width_extensions / 2); ...
    lineMoveAlong(edgeSideLeft(end, :), edgeSideLeft(end - 1, :), width_extensions / 2)];
sideLeftExtensionsOutline = makeExtensions(sideLeftPoly, vertices, extend_edge_vertices);

vertices = [...
    lineMoveAlong(edgeSideRight(1, :), edgeSideRight(2, :), width_extensions / 2); ...
    lineMoveAlong(edgeSideRight(end, :), edgeSideRight(end - 1, :), width_extensions / 2)];
sideRightExtensionsOutline = makeExtensions(sideRightPoly, vertices, extend_edge_vertices);

% recalculate side points near edge, since new lines may have been drawn
% if `align_fully_connect` is true
[edgeSideLeftVerticeIdx, ~] = getVertNearLine(sideLeftLattice, edgeSideLeft);
[edgeSideRightVerticeIdx, ~] = getVertNearLine(sideRightLattice, edgeSideRight);

% lattice extensions
topExtensions = makeExtensions(topPoly, topLattice(edgeTopLeftVerticeIdx | edgeTopRightVerticeIdx, :), extend_edge_vertices);
sideLeftExtensions = makeExtensions(sideLeftPoly, sideLeftLattice(edgeSideLeftVerticeIdx, :), extend_edge_vertices);
sideRightExtensions = makeExtensions(sideRightPoly, sideRightLattice(edgeSideRightVerticeIdx, :), extend_edge_vertices);

%% Get stitching nubs

% get edges
edgeTopLeft = topOutline([1 2], :);
edgeTopRight = topOutline([4 3], :);
edgeSideLeft = sideLeftOutline((sideOutsidePoints+1):end, :);
edgeSideRight = sideRightOutline((sideOutsidePoints+1):end, :);

topStitchNubs = [];
idx = getVertNearLine(topLattice, edgeTopLeft);
topStitchNubs = [topStitchNubs; topLattice(idx, :)];
idx = getVertNearLine(topLattice, edgeTopRight);
topStitchNubs = [topStitchNubs; topLattice(idx, :)];

idx = getVertNearLine(sideLeftLattice, edgeSideLeft);
sideLeftStitchNubs = sideLeftLattice(idx, :);

idx = getVertNearLine(sideRightLattice, edgeSideRight);
sideRightStitchNubs = sideRightLattice(idx, :);

%% Make 2D panels

% SETTINGS
width_lattice = 2;
width_outline_out = 2;
width_outline_stitch = 0;
width_nub_stitch = 2.5; % diameter
% width_extensions = 2; DEFINED ABOVE!

% each panel is a row of polyshapes comprised of:
% column 1: outside border, non-stitching side
% column 2: outside border, stitching side
% column 3: inside lattice
% column 4: stitch nubs
% column 5: grommet borders

% outlines are drawn at double width, then clipped, so that final
% dimensions exactly match

% make top panel
outline_out = [topOutline([2 3], :); ... % top
    nan nan; ...
    topOutline([4 1], :)]; % bottom

outline_stitch = [topOutline([1 2], :); ... % left
    nan nan; ...
    topOutline([3 4], :)]; % right

topPanel = capConstructPanel(topPoly, ...
    outline_out, width_outline_out, ... % top / bottom
    outline_stitch, width_outline_stitch, ... % left / right
    topLattice, width_lattice, ...
    topExtensionsOutline, topExtensions, width_extensions, ...
    topStitchNubs, width_nub_stitch);

% make left side panel
sideLeftPanel = capConstructPanel(sideLeftPoly, ...
    sideLeftOutline(1:sideOutsidePoints, :), width_outline_out, ... % bottom edge
    sideLeftOutline((sideOutsidePoints+1):end, :), width_outline_stitch, ... % top edge
    sideLeftLattice, width_lattice, ...
    sideLeftExtensionsOutline, sideLeftExtensions, width_extensions, ...
    sideLeftStitchNubs, width_nub_stitch);

% make right side panel
sideRightPanel = capConstructPanel(sideRightPoly, ...
    sideRightOutline(1:sideOutsidePoints, :), width_outline_out, ... % bottom edge
    sideRightOutline((sideOutsidePoints+1):end, :), width_outline_stitch, ... % top edge
    sideRightLattice, width_lattice, ...
    sideRightExtensionsOutline, sideRightExtensions, width_extensions, ...
    sideRightStitchNubs, width_nub_stitch);


%% Cut panels if needed

% SETTINGS
max_size = [270 270]; % mm
buffer = 0; % mm in overlap along boundaries

topPanel = cut(topPanel, max_size, buffer);
sideLeftPanel = cut(sideLeftPanel, max_size, buffer);
sideRightPanel = cut(sideRightPanel, max_size, buffer);

grommets = cutGrommets(grommets, topPanel, sideLeftPanel, sideRightPanel);

if debug
    figure;
    plot(topPanel);
    axis equal;
end

%% Make 2D panels

%save folder
stlfpath = 'stl/cap/';

% make height vector
names = {'outline-out' 'outline-stitch' 'lattice' 'stitch-nub' 'grommet'};

% construct STL files
for i = 1:size(topPanel, 1) % for each top piece
    for j = 1:min(length(names), size(topPanel, 2)) % for each part
        nm = sprintf('top%d-%s.stl', i, names{j});
        
        % is empty? skip
        if numboundaries(topPanel(i, j)) == 0
            if exist(nm, 'file')
                delete(nm);
            end
            continue;
        end
        
        % extrude the part
        [f, v] = extrude(topPanel(i, j), 0);
        
        % save the STL file
        saveModel([stlfpath nm], f, v);
        
    end
end

for i = 1:size(sideLeftPanel, 1) % for each left side piece
    for j = 1:min(length(names), size(sideLeftPanel, 2)) % for each part
        nm = sprintf('sideLeft-%s.stl', names{j});
        
        % is empty? skip
        if numboundaries(sideLeftPanel(i, j)) == 0
            if exist(nm, 'file')
                delete(nm);
            end
            continue;
        end
        
        % extrude the part
        [f, v] = extrude(sideLeftPanel(i, j), 0);
        
        % save the STL file
        saveModel([stlfpath nm], f, v);
    end
end

for i = 1:size(sideRightPanel, 1) % for each right side piece
    for j = 1:min(length(names), size(sideRightPanel, 2)) % for each part
        nm = sprintf('sideRight-%s.stl', names{j});
        
        % is empty? skip
        if numboundaries(sideRightPanel(i, j)) == 0
            if exist(nm, 'file')
                delete(nm);
            end
            continue;
        end
        
        % extrude the part
        [f, v] = extrude(sideRightPanel(i, j), 0);
        
        % save the STL file
        saveModel([stlfpath nm], f, v);
    end
end


%% get save positions of grommets and other elements in STLs for python script
IDXchStrapPoints = [11 12];
holders.sideLeftHolders=[];
holders.sideRightHolders=[];
holders.topHolders=[];
aux.sideLeftAUX=[];
aux.sideRightAUX=[];
aux.topAUX=[];

[STLcoords] = getStlCoordinates(grommets, holders, aux, sideLeftOutline, sideRightOutline, sideLeftPanel, topOutline, sideOutsidePoints, IDXchStrapPoints);
save('RefPts', '-struct', 'STLcoords')

