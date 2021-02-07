function generateLoad1010()
%GENERATELOAD1010 Summary of this function goes here
%   Detailed explanation goes here

global atlasViewer;

% little hack to allow script atlasviewer
atlasViewer.handles.figure.HandleVisibility = 'on';
figure(atlasViewer.handles.figure);

% get handles
handles = guidata(atlasViewer.handles.figure);
hObject = handles.menuItemShow10_5;

% register atlas to head size
AtlasViewerGUI('menuItemShowRefpts_Callback', hObject, struct(), handles);

end

