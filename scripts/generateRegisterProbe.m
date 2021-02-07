function generateRegisterProbe()
%GENERATEREGISTERPROBE Summary of this function goes here
%   Detailed explanation goes here

global atlasViewer;

% little hack to allow script atlasviewer
atlasViewer.handles.figure.HandleVisibility = 'on';
figure(atlasViewer.handles.figure);

% get handles
handles = guidata(atlasViewer.handles.figure);
hObject = handles.pushbuttonRegisterProbeToSurface;

% register atlas to head size
AtlasViewerGUI('pushbuttonRegisterProbeToSurface_Callback', hObject, struct(), handles);

end

