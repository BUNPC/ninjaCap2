%% CODE FOR NINJACAP GENERATION
% Instructions to the user:
% 0. set up and install everything as described in the github HELP file
% here: https://github.com/neuluce/ninjaCap
% 1. provide probe file with name "probe.SD" in dir \userinput\probe.SD
% 2. run generate.m with head circumference(HC) in cm as input argument: e.g. generate(56)
% 3. when the files are generated, a blender workspace will be opened. Run
% the default python script (click <Text> -> <Run Script>, or press Alt+P.
% 4. your generated and assembled cap is saved in the folder \print\. Here
% you also find a default 3D printing profile for Cura LULZBOT TAZ 6 using
% ninjaflex.
function [] = generateNinjaCap_AVfile(HC)

load('probe.SD','-mat')
dirSave = pwd;
rootDir = fileparts(which(mfilename));

%% Change the current folder to the folder of this m-file.
if(~isdeployed)
    rootDir = fileparts(which(mfilename));
    cd(rootDir);
end

panel_AVfile(SD,HC)
%% add this folder and subfolders to paths
% addpath(genpath(pwd))

%% AUTOMATED PANEL CREATION FROM SD FILES
% % if you renamed your files, update the path and names here
% if ~exist('HC','var')
%     HC = 56; % default cap size
% end
% circumference = HC; % cm (actual circumference, correction applied automatically)

%% CHECK FOR EXISTING AltasViewer FILE?
% if exist([pwd filesep 'atlasViewer.mat'], 'file')
%     disp('An "atlasViewer.mat" file already exists.');
%     prompt = ['To overwrite file press Y/y otherwise N/n and confirm with Enter:'];
%     x = input(prompt,'s');
%     if (x == 'y') || (x=='Y')
%         delete([pwd filesep 'atlasViewer.mat'])
%     else
%         error('Stopped')
%     end
% end

%% PREARE ATLASVIEWER.MAT

% % launch AtlasViewer
% % DEPENDENCY HOMER 3 >> getAppDir()
% %AtlasViewerGUI(pwd, [getAppDir() 'PACKAGES' filesep 'AtlasViewerGUI' filesep 'Data' filesep 'Colin'], 'userargs');     
% % appdir hardcoded
% avMainFilepath = which('AtlasViewerGUI.m');
% avAppDir = fileparts(avMainFilepath);
% avDataDir = filesepStandard([avAppDir, '/Data/Colin']);
% AtlasViewerGUI(pwd, avDataDir, 'userargs');
% 
% % scale head to corrected circumference
% generateScaleCircumference(circumference);
% 
% % load 10-10 points into AtlasViewer
% generateLoad10_5();
% 
% % align to head surface
% generateRegisterProbe();
% 
% % save state file
% generateSaveViewerState();

%% RUN PANEL.m
% cd(rootDir);
% copyfile( [dirSave filesep 'atlasViewer.mat'],'.' );
% panel('atlasViewer.mat', circumference);
drawnow

%% DELETE FILES
% delete([pwd filesep 'atlasViewer.mat']);
% delete([pwd filesep 'digpts.txt']);
% delete([pwd filesep 'digpts2mc.txt']);

%% RUN blender
system('wrkspace.blend')

%% COPY BLENDER OUTPUT TO WORKING DIRECTORY
if ~isdir([dirSave filesep 'print'])
    mkdir([dirSave filesep 'print']);
end
copyfile( ['print' filesep '*.stl'], [dirSave filesep 'print'] );
copyfile( ['print' filesep '*.curaprofile'], [dirSave filesep 'print'] );

cd(dirSave)

end