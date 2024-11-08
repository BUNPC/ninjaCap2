function [v, f, refpts, grommets, sprintList] = loadAtlasViewer_AVfile(SD)
%LOADATLASVIEWER Summary of this function goes here
%   Detailed explanation goes here

% load('probe.SD','-mat')
% short separation threshold
short_separation_threshold = 12;

% load file
% av = load(atlasviewer_file);

%% get surface
v = SD.mesh.vertices;
f = SD.mesh.faces;

%% get referene points
refpts = SD.Landmarks;

%%
% add grommet rot information as zeros if they are not available. (Check
% this if this needs to be here or somewhere else)
%  if ~isfield(av.probe,'SrcGrommetRot')
%  end

%% get grommets
grommets = struct('type', {}, 'rot',[],'flags', {}, 'posHead', {}, 'panel', {}, 'panelIndex', {}, 'posPanel', {}, 'MList', {}, 'optType', {});

%sources
for i = 1:SD.nSrcs
    grommet = struct('type', SD.SrcGrommetType{i}, 'rot', SD.SrcGrommetRot(i),'flags', {{}}, 'posHead', SD.optpos_reg(i, :), ...
        'panel', '', 'panelIndex', [], 'posPanel', [nan nan], 'MList', [], 'optType', 1);
    grommets = [grommets grommet]; %#ok<AGROW> % append
end
%detectors
for i = 1:SD.nDets
    % is short separation?
    p = SD.optpos_reg(SD.nSrcs + i, :);
    
    % make grommet
    grommet = struct('type', SD.DetGrommetType{i},'rot', SD.DetGrommetRot(i), 'flags', {{}}, 'posHead', p, ...
        'panel', '', 'panelIndex', [], 'posPanel', [nan nan], 'MList', [], 'optType', 2);
    
    % relevant sources
    rs = SD.MeasList(SD.MeasList(:, 2) == i, 1);
    
    % save relevant sources for measurement list
    grommet.MList = rs;
    
    dist = sqrt(sum(bsxfun(@minus, SD.optpos_reg(rs, :), p) .^ 2, 2));
    if ~isempty(rs) && all(dist < short_separation_threshold)
        % flag it as a short separation detector
        grommet.flags{end + 1} = 'short-separation';
        
        % edit all sources
        for j = rs'
            if ~ismember('has-short-separation', grommets(j).flags)
                grommets(j).flags{end + 1} = 'has-short-separation';
            end
        end
    end
    
    grommets = [grommets grommet]; %#ok<AGROW> % append
end
% dummy
% for i = 1:size(av.probe.al, 1)
nOptodes = SD.nSrcs + SD.nDets;
% for i = 1:size(av.probe.al,1)
%     % is short separation?
%     anchorOptIdx = av.probe.al{i,1};
%     if anchorOptIdx>nOptodes
%         iDummy = anchorOptIdx - nOptodes;
%         p = av.probe.optpos_reg(av.probe.al{i, 1}, :);
%         grommet = struct('type', av.probe.DummyGrommetType{iDummy},'rot', av.probe.DummyGrommetRot(iDummy), 'flags', {{av.probe.al{i, 2}}}, 'posHead', p, ...
%             'panel', '', 'panelIndex', [], 'posPanel', [nan nan], 'MList', [], 'optType', 3);
%         grommets = [grommets grommet]; %#ok<AGROW> % append
%     end
% end

for i = nOptodes+1:size(SD.optpos_reg,1)
     p = SD.optpos_reg(i,:);
     iDummy = i-nOptodes;
     grommet = struct('type', SD.DummyGrommetType{iDummy},'rot', SD.DummyGrommetRot(iDummy), 'flags', {{}}, 'posHead', p, ...
            'panel', '', 'panelIndex', [], 'posPanel', [nan nan], 'MList', [], 'optType', 3);
     grommets = [grommets grommet]; %#ok<AGROW> % append
end

% nOptodes = av.probe.nsrc + av.probe.ndet;
% for i = 1:size(av.probe.al, 1)
% 
%     % is short separation?
% 
%     anchorOptIdx = av.probe.al{i,1};
% 
%     if anchorOptIdx>nOptodes
% 
%         iDummy = anchorOptIdx - nOptodes;
% 
%         p = av.probe.optpos_reg(av.probe.al{i, 1}, :);
% 
%         grommet = struct('type', av.probe.DummyGrommetType{iDummy}, 'flags', {{av.probe.al{i, 2}}}, 'posHead', p, ...
% 
%         'panel', '', 'panelIndex', [], 'posPanel', [nan nan], 'MList', [], 'optType', 3);
% 
%         grommets = [grommets grommet]; %#ok<AGROW> % append
% 
%     end
% 
% end

%% check and unify different "Dummy" types:
% Set #DUMMY, 'none', 'None', 'NONE', and any other ID without # to '#NONE'
for gg = 1:numel(grommets)
   switch grommets(gg).type
       case {'none','None','NONE','#DUMMY','#NONE'}
          grommets(gg).type = '#DUMMY';
   end
   if ~strfind(grommets(gg).type, '#')
       grommets(gg).type = '#DUMMY';
   end
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
    ht = text(pos(ii, 1), pos(ii, 2), pos(ii, 3), refpts.labels{ii});
    %set(hp, 'markersize', 16);
    
    h2 = trisurf(fs, pos(ii, 1) + vs(:, 1), pos(ii,2) + vs(:, 2), pos(ii, 3) + vs(:, 3));
    alpha(0.5);
    set(h2, 'linestyle', 'none');
end
hold off;

%% get spring list
sprintList = SD.SpringList;

end

