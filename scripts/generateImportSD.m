function generateImportSD(filename)
%GENERATEIMPORTSD Summary of this function goes here
%   Detailed explanation goes here

% MODIFIED VERSION OF: menuItemImportProbe_Callback (no longer a UI prompt)

global atlasViewer

probe        = atlasViewer.probe;
refpts       = atlasViewer.refpts;
headsurf     = atlasViewer.headsurf;
labelssurf   = atlasViewer.labelssurf;
fwmodel      = atlasViewer.fwmodel;
imgrecon     = atlasViewer.imgrecon;
digpts       = atlasViewer.digpts;

probe = getProbe(probe, filename, headsurf, digpts, refpts);
probe = viewProbe(probe, 'unregistered');

% This is done to not display dummy points by default. It does nothing 
% if the method isn't spring registration.
probe = setProbeDisplay(probe,headsurf);

% get pathname
[pathname, nm, ~] = fileparts(which(filename));

% add grommet information to probe (new SD structure)
sdbuf = load(fullfile(pathname, [nm '.SD']), '-mat');
probe.SrcGrommetType = sdbuf.SD.SrcGrommetType;
probe.DetGrommetType = sdbuf.SD.DetGrommetType;
probe.DummyGrommetType = sdbuf.SD.DummyGrommetType;

% add grommet rotation information to probe. If it is not available in the 
% loaded SD file then make rotation of the grommets to zero.
if isfield(sdbuf.SD,'SrcGrommetRot')
    probe.SrcGrommetRot = sdbuf.SD.SrcGrommetRot;
else
    probe.SrcGrommetRot = zeros(size(sdbuf.SD.SrcGrommetType));
end
if isfield(sdbuf.SD,'DetGrommetRot')
    probe.DetGrommetRot = sdbuf.SD.DetGrommetRot;
else
    probe.DetGrommetRot = zeros(size(sdbuf.SD.DetGrommetType));
end
if isfield(sdbuf.SD,'DummyGrommetRot')
    probe.DummyGrommetType = sdbuf.SD.DummyGrommetRot;
else
    probe.DummyGrommetRot = zeros(size(sdbuf.SD.DummyGrommetType));
end

atlasViewer.probe        = probe;
atlasViewer.dirnameProbe = pathname;
atlasViewer.labelssurf   = labelssurf;
atlasViewer.digpts       = digpts;
atlasViewer.fwmodel      = fwmodel;
atlasViewer.imgrecon     = imgrecon;

