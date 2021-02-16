function panel(atlasViewerFile, circumferenceDesired)

%% Parameters
debug = true;
sWidth = 0.375;
sExtensions = 1;
%% Load head model
% Use the view state (atlasViewer.mat).
% fpath = 'D:\Office\Research\Software - Scripts\Matlab\ninjaCap\atlasViewer.mat';
% load the AtlasViewer state file, including a head model and grommets
[vHead, fHead, refpts, grommets, springList] = loadAtlasViewer(atlasViewerFile);

%% Get 3D seam points
% segmentHead function generates 3D and 2D outline points for left, top and right panels.
% But we use only 3D outline points since we are going to generate 2D outlines in a different way
% using spring relaxation method. Also, top3DVertices given by segmentHead
% function is not accurate. Correct one is generated below using sideLeft3DVertices and
% sideRight3DVertices. In the future we might want to change segmentHead
% function. But for now it works.

fname = 'svg/side_piece2.svg';
pathid = 'side_piece';
% side panel outline and number of side panel points that define the outside vertices
[side_piece, sideOutsidePoints] = svgImport(fname, pathid);
[sideLeftOutline, sideRightOutline, topOutline, sideLeft3DVertices, sideRight3DVertices, top3DVertices] = segmentHead(vHead, fHead, refpts, side_piece);
ind = find(isnan(sideLeft3DVertices(:,1)) == 0);
sideLeft3DVertices = sideLeft3DVertices(ind,:);
ind = find(isnan(sideRight3DVertices(:,1)) == 0);
sideRight3DVertices = sideRight3DVertices(ind,:);
top3DVertices = [flip(sideLeft3DVertices,1); sideRight3DVertices; sideLeft3DVertices(end,:)];
topOutline_corner_idx = [1;size(sideLeft3DVertices,1);size(sideLeft3DVertices,1)+1;size(sideLeft3DVertices,1)+size(sideRight3DVertices,1)];
%% Assign grommets to panels

left_cutoff = max(max(sideLeft3DVertices(:,1)));
right_cutoff = min(sideRight3DVertices(:,1));
left_idx = [];
top_idx = [];
right_idx = [];
for j = 1:length(grommets)
    if grommets(j).posHead(1) <= left_cutoff
        grommets(j).panel = 'sideLeft';
        left_idx = [left_idx; j];
    elseif grommets(j).posHead(1) < right_cutoff
        grommets(j).panel = 'top';
        top_idx = [top_idx; j];
    else
        grommets(j).panel = 'sideRight';
        right_idx = [right_idx; j];
    end
end
%% get 2D panel seams and grommet positions on 2D panels

% get refPts connected neighbours information
[rightSide_refPts_neighbours, leftSide_refPts_neighbours, top_refPts_neighbours] =  get_refpts_neighbours();

right_cutoff = min(sideRight3DVertices(:,1));
refpts_sideRight_idx = find(refpts.pos(:,1) >= right_cutoff);
[sideRightOutline, grommets, ~, rightEarSlit] = mapPanelandGrommets(grommets, springList, refpts, refpts_sideRight_idx, sideRight3DVertices, rightSide_refPts_neighbours, 'sideRight');

left_cutoff = min(sideLeft3DVertices(:,1));
refpts_sideLeft_idx = find(refpts.pos(:,1) <= left_cutoff);
[sideLeftOutline, grommets, ~, leftEarSlit] = mapPanelandGrommets(grommets, springList, refpts, refpts_sideLeft_idx, sideLeft3DVertices, leftSide_refPts_neighbours, 'sideLeft');

refpts_top_idx = find(refpts.pos(:,1) > left_cutoff & refpts.pos(:,1) < right_cutoff);
[topOutline, grommets, topIdx, Cz_pos] = mapPanelandGrommets(grommets, springList, refpts, refpts_top_idx, top3DVertices, top_refPts_neighbours, 'top',topOutline_corner_idx);

%% Add the side piece to side panels

fixedPoints = [sideRightOutline(1,1:2); sideRightOutline(end,1:2)];
movingPoints = [side_piece(end, :); side_piece(1, :)];

tform = fitgeotrans(movingPoints, fixedPoints, 'nonreflectivesimilarity');

sideRight_side_piece = transformPointsForward(tform, side_piece);
sideRight_side_piece = (sideRight_side_piece);
sideRightIdx = [1;size(sideRightOutline,1)];
sideRightOutline = [sideRightOutline(:,1:2); sideRight_side_piece];

side_piece(:,1) = -1*side_piece(:,1);
fixedPoints = [sideLeftOutline(1,1:2); sideLeftOutline(end,1:2)];
movingPoints = [side_piece(end, :); side_piece(1, :)];
tform = fitgeotrans(movingPoints, fixedPoints, 'nonreflectivesimilarity');

sideLeft_side_piece = transformPointsForward(tform, side_piece);
sideLeft_side_piece = (sideLeft_side_piece);
sideLeftIdx = [1;size(sideLeftOutline,1)];
sideLeftOutline = [sideLeftOutline(:,1:2); sideLeft_side_piece];
topOutline = topOutline(:,1:2);
%%
% plot 2D panels with grommets positions
figure;
plot(sideRightOutline(:,1), sideRightOutline(:,2));
hold on;
a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#NONE'), grommets)).posPanel);
if ~isempty(a)
    scatter(a(:, 1), a(:, 2));
end
hold off;
axis equal

figure;
plot(sideLeftOutline(:,1), sideLeftOutline(:,2));
hold on;
a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#NONE'), grommets)).posPanel);
if ~isempty(a)
    scatter(a(:, 1), a(:, 2));
end
hold off;
axis equal

figure;
plot(topOutline(:,1), topOutline(:,2));
hold on;
a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#NONE'), grommets)).posPanel);
if ~isempty(a)
    scatter(a(:, 1), a(:, 2));
end
hold off;
axis equal
%%
% fill cap with hexagons
C = fillCapWithHexagons(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx);

% get panel (poly shpaes) from outlinea and hexagonal connections
[sideLeftPanel, sideRightPanel, topPanel] = getPanels(C, sWidth);

%% add marker grommet automatically. It is decided that it is best to add marker grommet while designing probe instead of adding automatically.
% nGrommets = length(grommets);
% grommets(nGrommets+1).type = '#PMRK1';
% grommets(nGrommets+1).rot = 0;
% grommets(nGrommets+1).panel = 'top';
% grommets(nGrommets+1).posPanel= Cz_pos(1:2);
% grommets(nGrommets+1).flags = {};
% grommets(nGrommets+1).optType = 3;
%%
% move grommet positions as panels position is moved in filling cap with
% hexagons
outlines_min = [min(topOutline(:,1)) min(topOutline(:,2)); min(sideLeftOutline(:,1)) min(sideLeftOutline(:,2)); min(sideRightOutline(:,1)) min(sideRightOutline(:,2))];
for u = 1:length(grommets)
    if strcmp(grommets(u).panel,'sideRight')
        grommets(u).posPanel = grommets(u).posPanel-[outlines_min(3,1) outlines_min(3,2)];
    end
    if strcmp(grommets(u).panel,'sideLeft')
        grommets(u).posPanel =grommets(u).posPanel-[outlines_min(2,1) outlines_min(2,2)];
    end
    if strcmp(grommets(u).panel,'top')
        grommets(u).posPanel = grommets(u).posPanel-[outlines_min(1,1) outlines_min(1,2)];
    end
end

%% Add extensions if needed
if sExtensions == 1
   [LE,RE,TLE,TRE] =  secondOrderExtensions(grommets,C);
   
   ext_lattice_left = [];
   ext_lattice_right = [];
   ext_lattice_top = [];
   
   % add left extensions and also make corresponidng top extensions wider
   for u = 1:length(LE)
       if LE(u) == 1
           for iL = 1:4:length(C(2).eExtensionsTopRot(u).pts)
               ext_lattice_left = [ext_lattice_left; [C(2).eExtensionsTopRot(u).pts(iL) C(2).eExtensionsTopRot(u).pts(iL+1)];...
                       [C(2).eExtensionsTopRot(u).pts(iL+2) C(2).eExtensionsTopRot(u).pts(iL+3)]; [NaN NaN]];
               ext_lattice_top = [ext_lattice_top; [C(2).eExtensionsTop(u).pts(iL) C(2).eExtensionsTop(u).pts(iL+1)];...
                       [C(2).eExtensionsTop(u).pts(iL+2) C(2).eExtensionsTop(u).pts(iL+3)]; [NaN NaN]];
           end
       end
   end
   
   % add right extensions and also make corresponidng top extensions wider
    for u = 1:length(RE)
       if RE(u) == 1
           for iL = 1:4:length(C(3).eExtensionsTopRot(u).pts)
                ext_lattice_right = [ ext_lattice_right; [C(3).eExtensionsTopRot(u).pts(iL) C(3).eExtensionsTopRot(u).pts(iL+1)];...
                       [C(3).eExtensionsTopRot(u).pts(iL+2) C(3).eExtensionsTopRot(u).pts(iL+3)]; [NaN NaN]];
               ext_lattice_top = [ext_lattice_top; [C(3).eExtensionsTop(u).pts(iL) C(3).eExtensionsTop(u).pts(iL+1)];...
                       [C(3).eExtensionsTop(u).pts(iL+2) C(3).eExtensionsTop(u).pts(iL+3)]; [NaN NaN]];
           end
       end
    end
   
    % add topLeft extensions and also make corresponidng left extensions wider
    for u = 1:length(TLE)
       if TLE(u) == 1
           for iL = 1:4:length(C(2).eExtensionsRot(u).pts)
                ext_lattice_top = [ ext_lattice_top; [C(2).eExtensionsRot(u).pts(iL) C(2).eExtensionsRot(u).pts(iL+1)];...
                       [C(2).eExtensionsRot(u).pts(iL+2) C(2).eExtensionsRot(u).pts(iL+3)]; [NaN NaN]];
               ext_lattice_left = [ext_lattice_left; [C(2).eExtensions(u).pts(iL) C(2).eExtensions(u).pts(iL+1)];...
                       [C(2).eExtensions(u).pts(iL+2) C(2).eExtensions(u).pts(iL+3)]; [NaN NaN]];
           end
       end
    end
    
     % add topRight extensions and also make corresponidng right extensions wider
    for u = 1:length(TRE)
       if TRE(u) == 1
           for iL = 1:4:length(C(3).eExtensionsRot(u).pts)
                ext_lattice_top = [ ext_lattice_top; [C(3).eExtensionsRot(u).pts(iL) C(3).eExtensionsRot(u).pts(iL+1)];...
                       [C(3).eExtensionsRot(u).pts(iL+2) C(3).eExtensionsRot(u).pts(iL+3)]; [NaN NaN]];
               ext_lattice_right = [ext_lattice_right; [C(3).eExtensions(u).pts(iL) C(3).eExtensions(u).pts(iL+1)];...
                       [C(3).eExtensions(u).pts(iL+2) C(3).eExtensions(u).pts(iL+3)]; [NaN NaN]];
           end
       end
    end
    
    if ~isempty(ext_lattice_left)
       ext_left_poly = polybuffer(ext_lattice_left, 'lines', 1, 'JointType', 'square');
       sideLeftPanel = union(sideLeftPanel,ext_left_poly );       
    end
    
    if ~isempty(ext_lattice_right)
       ext_right_poly = polybuffer(ext_lattice_right, 'lines', 1, 'JointType', 'square');
       sideRightPanel = union(sideRightPanel,ext_right_poly );       
    end
    
    if ~isempty(ext_lattice_top)
       ext_top_poly = polybuffer(ext_lattice_top, 'lines', 1, 'JointType', 'square');
       topPanel = union(topPanel,ext_top_poly );       
    end
end

%% Make overlapping struts wider
left_overlap_lattice = [];
for u = 1:size(C(2).eSeamExt,1)
    left_overlap_lattice = [left_overlap_lattice; [C(2).eSeamExt(u,1:2);C(2).eSeamExt(u,3:4)]; [NaN NaN]];
end
left_overlap_poly = polybuffer(left_overlap_lattice, 'lines', 1, 'JointType', 'square');
sideLeftPanel  = union(sideLeftPanel,left_overlap_poly);

right_overlap_lattice = [];
for u = 1:size(C(3).eSeamExt,1)
    right_overlap_lattice = [right_overlap_lattice; [C(3).eSeamExt(u,1:2);C(3).eSeamExt(u,3:4)]; [NaN NaN]];
end
right_overlap_poly = polybuffer(right_overlap_lattice, 'lines', 1, 'JointType', 'square');
sideRightPanel  = union(sideRightPanel,right_overlap_poly);

top_overlap_lattice = [];
for u = 1:size(C(2).eSeamExtRot,1)
    top_overlap_lattice = [top_overlap_lattice; [C(2).eSeamExtRot(u,1:2);C(2).eSeamExtRot(u,3:4)]; [NaN NaN]];
end

for u = 1:size(C(3).eSeamExtRot,1)
    top_overlap_lattice = [top_overlap_lattice; [C(3).eSeamExtRot(u,1:2);C(3).eSeamExtRot(u,3:4)]; [NaN NaN]];
end
top_overlap_poly = polybuffer(top_overlap_lattice, 'lines', 1, 'JointType', 'square');
topPanel  = union(topPanel,top_overlap_poly);

%% 

% add outlines (side_piece) to polygons also extend top outline by 5mm to
% overlap with side outlines
sideLeftOutline_poly = polybuffer(sideLeftOutline(sideLeftIdx(2)+1:end,:)-min(sideLeftOutline,[],1), 'lines', 1, 'JointType', 'square');
sideRightOutline_poly = polybuffer(sideRightOutline(sideRightIdx(2)+1:end,:)-min(sideRightOutline,[],1), 'lines', 1, 'JointType', 'square');
topOutline_side1 = topOutline(topIdx(2):topIdx(3),:);
topOutline_side2 = [topOutline(topIdx(4):end,:); topOutline(1,:)];

temp_OL = sideLeftOutline(sideLeftIdx(2)-1:sideLeftIdx(2),:);
theta = asin( (temp_OL(2,1)- temp_OL(1,1)) / norm(temp_OL(2,:)-temp_OL(1,:)));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
side1_ext1 = (R * (sideLeftOutline(sideLeftIdx(2)+2,:)-sideLeftOutline(sideLeftIdx(2)+1,:))' + topOutline_side1(1,:)')';

temp_OL = sideRightOutline(sideRightIdx(2)-1:sideRightIdx(2),:);
theta = asin( (temp_OL(2,1)- temp_OL(1,1)) / norm(temp_OL(2,:)-temp_OL(1,:)));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
side1_ext2 = (R * (sideRightOutline(sideRightIdx(2)+2,:)-sideRightOutline(sideRightIdx(2)+1,:))' + topOutline_side1(end,:)')';

temp_OL = sideLeftOutline(1:2,:);
theta = asin( (temp_OL(2,1)- temp_OL(1,1)) / norm(temp_OL(2,:)-temp_OL(1,:)));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
side2_ext2 = (R * (sideLeftOutline(end-1,:)-sideLeftOutline(end,:))' + topOutline_side2(end,:)')';

temp_OL = sideRightOutline(1:2,:);
theta = asin( (temp_OL(2,1)- temp_OL(1,1)) / norm(temp_OL(2,:)-temp_OL(1,:)));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
side2_ext1 = (R * (sideRightOutline(end-1,:)-sideRightOutline(end,:))' + topOutline_side2(1,:)')';

topOutline_side1 = [side1_ext1; topOutline_side1; side1_ext2];
topOutline_side2 = [side2_ext1; topOutline_side2; side2_ext2];
topOutline_lattice = [topOutline_side1; [NaN NaN]; topOutline_side2];
topOutline_poly =  polybuffer(topOutline_lattice-min(topOutline,[],1), 'lines', 1, 'JointType', 'square');

%%

% Make data format as previous ninjaCap so that from here on we can use
% prevoius ninjacap code
sideLeftPanel(3) = sideLeftPanel;
sideLeftPanel(1) = sideLeftOutline_poly;
sideLeftPanel(4) = sideLeftPanel(2);
sideLeftPanel(5) = sideLeftPanel(2);

sideRightPanel(3) = sideRightPanel;
sideRightPanel(1) = sideRightOutline_poly;
sideRightPanel(4) = sideRightPanel(2);
sideRightPanel(5) = sideRightPanel(2);

topPanel(3) = topPanel;
topPanel(1) = topOutline_poly;
topPanel(4) = topPanel(2);
topPanel(5) = topPanel(2);


%% Cut panels if needed

% SETTINGS
max_size = [270 270]; % mm
buffer = 5; % mm in overlap along boundaries

topPanel = cut(topPanel, max_size, buffer);
sideLeftPanel = cut(sideLeftPanel, max_size, buffer);
sideRightPanel = cut(sideRightPanel, max_size, buffer);

grommets = cutGrommets(grommets, topPanel, sideLeftPanel, sideRightPanel);


%% Debug grommets

if debug
    figure;
    plot(topPanel);
    hold on;
    a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'top') && ~strcmp(x.type, '#NONE') && ~ismember('short-separation', x.flags), grommets)).posPanel);
    if ~isempty(a)
        scatter(a(:, 1), a(:, 2));
    end
    hold off;
    axis equal;
    
    figure;
    plot(sideLeftPanel);
    hold on;
    a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideLeft') && ~strcmp(x.type, '#NONE') && ~ismember('short-separation', x.flags), grommets)).posPanel);
    if ~isempty(a)
        scatter(a(:, 1), a(:, 2));
    end
    hold off;
    axis equal;
    
    figure;
    plot(sideRightPanel);
    hold on;
    a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel, 'sideRight') && ~strcmp(x.type, '#NONE') && ~ismember('short-separation', x.flags), grommets)).posPanel);
    if ~isempty(a)
        scatter(a(:, 1), a(:, 2));
    end
    hold off;
    axis equal;
end

%% Make 2D panels

%save folder
stlfpath = 'stl/cap/';
if ~isfolder('stl/cap')
    mkdir stl cap
end

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
        
        
        TR = triangulation(f,v);
        % save the STL file
        stlwrite( TR,[stlfpath nm],'text');
        
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
        
        TR = triangulation(f,v);
        % save the STL file
        stlwrite( TR,[stlfpath nm],'text');
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
        
        TR = triangulation(f,v);
        % save the STL file
        stlwrite( TR,[stlfpath nm],'text');
    end
end



%% Save positions of grommets and other elements in STLs for python script
IDXchStrapPoints = [10 11];
holders.sideLeftHolders=[];
holders.sideRightHolders=[];
holders.topHolders=[];
aux.sideLeftAUX=[];
aux.sideRightAUX=[];
aux.topAUX=[];

shifted_sideLeftOutline = sideLeftOutline-min(sideLeftOutline,[],1);
shifted_sideRightOutline = sideRightOutline-min(sideRightOutline,[],1);
shifted_topOutline = topOutline-min(topOutline,[],1);
[STLcoords] = getStlCoordinates(grommets, holders, aux, shifted_sideLeftOutline, shifted_sideRightOutline, sideLeftPanel, shifted_topOutline, sideOutsidePoints, IDXchStrapPoints);
% add ear slit pos
STLcoords.sideLeftEar = leftEarSlit(1:2) - min(sideLeftOutline,[],1) -[10 10];
STLcoords.sideRightEar = rightEarSlit(1:2) - min(sideRightOutline,[],1) -[10 10];

% add marker optode at Cz position
% STLcoords.
save('stl\cap\RefPts', '-struct', 'STLcoords')

%% RUN SD SEPARATION VALIDATION
sdseps= sdsepValid(grommets);

end