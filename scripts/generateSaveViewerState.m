function generateSaveViewerState()
%GENERATESAVEVIEWERSTATE Summary of this function goes here
%   Detailed explanation goes here

global atlasViewer;

% little hack to allow script atlasviewer
atlasViewer.handles.figure.HandleVisibility = 'on';
figure(atlasViewer.handles.figure);

% get handles
handles = guidata(atlasViewer.handles.figure);
hObject = handles.menuItemSaveViewerState;

% register atlas to head size
AtlasViewerGUI('menuItemSaveViewerState_Callback', hObject, struct(), handles);

end

