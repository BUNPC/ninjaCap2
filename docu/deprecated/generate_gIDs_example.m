% load relevant SD file and generate gIDs.mat
% MAY 2020
dir = 'C:\Users\avolu\Google Drive\NinjaCap_Printing\LindsayButler';
sdname = 'AutismLanguage.sd';
load(fullfile(dir,sdname), '-mat');

% Currently available grommets (see grommet_lookup.pdf)
%   #NFLPS  - NIRS: Fiber-based optode holder (CW6 fibers) with Low
%           Profile and single Short-separation hole
%   #NFHPS  - NIRS: Fiber-based optode holder (CW6 fibers) with High 
%           Profile and single Short-separation hole
%   #NFWMS  - NIRS: Fiber-based Wide optode holder (CW6 fibers) 
%           with Multiple Short-separation holes
%   #NFDSO  - NIRS: Fiber-based Dual Ss Optode holder
%   #NOND1  - NIRS: Open fNIRS Dual optode holder for optodes in printed case
%   #PMRK1  - Position MaRKer in small cross shape, symmetric
%   #EBPAS  - EEG electrode holder from Brain Products: ActiCap Snap
%   #EECEH  - EEG Easy Cap Electrode Holder: Ring for commercial electrode holders
%   #DUMMY  - Dummy, nothing placed (use this for e.g. Anchor Dummy Optodes)

% assign (special) optodes here, e.g. a high profile optode to reduce pressure over frontal regions
optode_frontal_list = [1 2 3 9 10 11 17 18 19 21 23 24 33 34 35 37 39 40]; 
optode_frontal_type = '#NFHPS';

% default optode holder
optode_default = '#NFLPS';%'#NOND1'; 

rhoSD_ssThresh = 15;
ml = SD.MeasList;
lst = find(ml(:,4)==1);
rhoSD = zeros(length(lst),1);
posM = zeros(length(lst),3);
for iML = 1:length(lst)
    rhoSD(iML) = sum((SD.SrcPos(ml(lst(iML),1),:) - SD.DetPos(ml(lst(iML),2),:)).^2).^0.5;
    posM(iML,:) = (SD.SrcPos(ml(lst(iML),1),:) + SD.DetPos(ml(lst(iML),2),:)) / 2;
end
lstSS = lst(find(rhoSD<=rhoSD_ssThresh));

for i = 1:size(lstSS,1)
    foo = ml(lstSS(i),:);
    ss_det(i) = foo(2);
end

for i = 1:(size(SD.SrcPos,1)+size(SD.DetPos,1)); gIDs{i}=optode_default; end % assign all sources and detectors to default optode profile

% provide the list of optodes that you want to be printed in a different
% profile. Eg. for frontal high profile optodes.
if exist('optode_frontal_list') == 1
    for i = 1:size(optode_frontal_list,2); gIDs{optode_frontal_list(i)}=optode_frontal_type; end  % set frontal optodes to high profile
end

for i = 1:size(ss_det,2); gIDs{size(SD.SrcPos,1)+ss_det(i)}='#DUMMY'; end  % set ss detectors as dummy (no print)

for i = size(gIDs,2)+1 : size(gIDs,2)+ size(SD.DummyPos,1); gIDs{i}='#DUMMY'; end  % set dummy optodes as dummy (no print)




save(fullfile(dir, 'gIDs.mat'), 'gIDs')
