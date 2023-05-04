function generateScaleCircumference(target_circumference)
%GENERATESCALECIRCUMFERENCE Summary of this function goes here
%   Detailed explanation goes here

% correction factor based on weird inaccuracy in AtlasViewerScaling
c = 1.025;
i = 1.0321;
l = 0.9721;

% linear scaling factor based on ration between target circumference and 
% standard circumference
r = target_circumference / 56;
scaling = 1; % it was 1.09 but that was making cap of size 61cm.
circumference   = 56 * c * r * scaling;
NzCzIz          = 36.43 * i * r * scaling;
LPACzRPA        = 35.18 * l * r * scaling;

fprintf('Circumference: %.2f\n' ,circumference);
fprintf('Nz-Cz-Iz: %.2f\n', NzCzIz);
fprintf('LPA-Cz-RPA: %.2f\n', LPACzRPA);

% MODIFIED VERSION OF: menuItemRegisterAtlasToHeadSize_Callback (no longer prompt)

global atlasViewer;

digpts       = atlasViewer.digpts;
refpts       = atlasViewer.refpts;

headsizeParams = {num2str(circumference), num2str(NzCzIz), num2str(LPACzRPA)};
digpts.headsize = setHeadsize(digpts.headsize, headsizeParams);
digpts = calcDigptsFromHeadsize(digpts, refpts);

atlasViewer.digpts = digpts;

% little hack to allow script atlasviewer
atlasViewer.handles.figure.HandleVisibility = 'on';
figure(atlasViewer.handles.figure);

% get handles
handles = guidata(atlasViewer.handles.figure);
hObject = handles.menuItemRegisterAtlasToHeadSize;

% register atlas to head size
AtlasViewerGUI('menuItemRegisterAtlasToDigpts_Callback', hObject, struct(), handles);

atlasViewer.probe.optpos_reg = atlasViewer.probe.optpos_reg*target_circumference/61.3;
atlasViewer.probe.srcpos = atlasViewer.probe.srcpos*target_circumference/61.3;
atlasViewer.probe.detpos = atlasViewer.probe.detpos*target_circumference/61.3;
if isfield(atlasViewer.probe,'dummypos')
    atlasViewer.probe.dummypos = atlasViewer.probe.dummypos*target_circumference/61.3;
end
atlasViewer.probe.optpos = atlasViewer.probe.optpos*target_circumference/61.3;
atlasViewer.probe.registration.sl(:,3) = atlasViewer.probe.registration.sl(:,3)*target_circumference/61.3;
end

