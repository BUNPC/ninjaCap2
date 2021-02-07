% correct template by enforcing symmetry

%% TOP PANEL
tpanel = load('template/map-top.mat');
tpanelNew = tpanel;
% points on center line
[~, cntrlineIDX] = ismember({'CPz','Pz','FOz','Oz','Iz','FPz','AFz','Fz','FCz','Cz'}, tpanel.labels);
% Left/right pairs
[~, lIDX] = ismember({'P1','O1','FP1','AF3','F1','FC1','C1','CP1','PO3','I1'}, tpanel.labels);
[~, rIDX] = ismember({'P2','O2','FP2','AF4','F2','FC2','C2','CP2','PO4','I2'}, tpanel.labels);
lrIDX = [lIDX; rIDX];
      % corresponding center points
[~, cIDX] = ismember({'Pz','Oz','FPz','AFz','Fz','FCz','Cz','CPz','FOz','Iz'}, tpanel.labels);

%% enforce symmetry of left/right points
% average horizontal distance to corresponding midpoint
hdist = (abs(tpanel.points(lrIDX(1,:),1)-tpanel.points(cIDX,1)) + ...
    abs(tpanel.points(lrIDX(2,:),1)-tpanel.points(cIDX,1)))/2;
      
%% correct center line (X coordinate):
cx = tpanel.outline(3,1)/2;
tpanelNew.points(cntrlineIDX,1)=cx;

%% horizontal correction: re-set L/R points around center line 
tpanelNew.points(lrIDX(1,:),1) = cx - hdist;
tpanelNew.points(lrIDX(2,:),1) = cx + hdist;

%% vertical correction: re-set y coordinate to average y of L/R points
tpanelNew.points(lrIDX(1,:),2) = (tpanelNew.points(lrIDX(1,:),2)+tpanelNew.points(lrIDX(2,:),2))/2;
tpanelNew.points(lrIDX(2,:),2) = (tpanelNew.points(lrIDX(1,:),2)+tpanelNew.points(lrIDX(2,:),2))/2;

figure
plot(tpanel.panel)
hold on 
plot(tpanel.points(:,1), tpanel.points(:,2), 'xk')
plot(tpanelNew.points(:,1), tpanelNew.points(:,2), 'or') 
title('top panel')

%% save
save('template/map-top.mat', '-struct', 'tpanelNew')

%% SIDE PANELS
sRpanel = load('template/map-sideRight.mat');
sRpanelNew = sRpanel;
sLpanel = load('template/map-sideLeft.mat');
sLpanelNew = sLpanel;
%% assign channels (all identical but last two)
[~, lidx] = ismember({'AF7','F3','F5','F7','F9','FC3','FC5','FT7','FT9','C3','C5','T7','LPA','CP3','CP5','P3','P5','P9','PO7','TP7','P7','PO9'}, sLpanel.labels);
[~, ridx] = ismember({'AF8','F4','F6','F8','F10','FC4','FC6','FT8','FT10','C4','C6','T8','RPA','CP4','CP6','P4','P6','P10','PO8','TP8','P8','PO10'}, sRpanel.labels);
chidx = [ridx; lidx];

%% find center of bounding boxes
[xL, yL] = boundingbox(sLpanel.panel);
cSL = [mean(xL) mean(yL)];
[xR, yR] = boundingbox(sRpanel.panel);
cSR = [mean(xR) mean(yR)];

%% calculate average position for each L-R point using center of bounding
% box for reflecting: mirror points
Rpoints(:,2) = sRpanel.points(:,2);
Rpoints(:,1) = sRpanel.points(:,1)-cSR(1);
Lpoints(:,2) = sLpanel.points(:,2);
Lpoints(:,1) = sLpanel.points(:,1)-cSL(1);

%% Save average positions: 
sRpanelNew.points(chidx(1,:),1) = (Rpoints(chidx(1,:),1)-Lpoints(chidx(2,:),1))/2 + cSR(1);
sRpanelNew.points(chidx(1,:),2) = (Rpoints(chidx(1,:),2)+Lpoints(chidx(2,:),2))/2;
sLpanelNew.points(chidx(2,:),1) = (Lpoints(chidx(2,:),1)-Rpoints(chidx(1,:),1))/2 + cSL(1);
sLpanelNew.points(chidx(2,:),2) = (Lpoints(chidx(2,:),2)+Rpoints(chidx(1,:),2))/2;


figure
plot(sRpanel.panel)
hold on
plot(sRpanel.points(:,1), sRpanel.points(:,2), 'kx')
plot(sRpanelNew.points(:,1), sRpanelNew.points(:,2), 'or')
title('right side panel')
figure
plot(sLpanel.panel)
hold on
plot(sLpanel.points(:,1), sLpanel.points(:,2), 'kx')
plot(sLpanelNew.points(:,1), sLpanelNew.points(:,2), 'or')
title('left side panel')

%% save
save('template/map-sideLeft.mat', '-struct', 'sLpanelNew')
save('template/map-sideRight.mat', '-struct', 'sRpanelNew')
