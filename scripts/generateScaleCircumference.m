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

circumference   = 56 * c * r * 1.09;
NzCzIz          = 36.43 * i * r * 1.09;
LPACzRPA        = 35.18 * l * r * 1.09;

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

end

