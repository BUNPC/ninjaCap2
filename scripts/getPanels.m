function [sideLeftPanel, sideRightPanel, topPanel] = getPanels(o,sWidth, outlines_min)



%% create panels

% sidePanels
for u = 2:3
    v = o(u).v;
    e = o(u).e;
    eSeamExt = o(u).eSeamExt;
    Outline = [o(u).xoutline o(u).youtline];
    v(:,1) = v(:,1) + outlines_min(u,1);
    v(:,2) = v(:,2) + outlines_min(u,2);
    eSeamExt(:,[1 3]) = eSeamExt(:,[1 3]) + outlines_min(u,1);
    eSeamExt(:,[2 4]) = eSeamExt(:,[2 4]) + outlines_min(u,2);

    lattice = [];
    for iE=1:size(e,1)
        lattice = [lattice; v(e(iE,1),:); v(e(iE,2),:); [NaN NaN] ];
    end
    for iS=1:size(eSeamExt,1)
        lattice = [lattice; eSeamExt(iS,1:2); eSeamExt(iS,3:4); [NaN NaN] ];
    end
    
    if u == 2
        sideLeftPanel = polybuffer(lattice, 'lines', sWidth, 'JointType', 'square');
    end
     if u == 3
        sideRightPanel = polybuffer(lattice, 'lines', sWidth, 'JointType', 'square');
    end
end

%topPanel
v = o(1).v;
e = o(1).e;
Outline = [o(u).xoutline o(u).youtline];
v(:,1) = v(:,1) + outlines_min(1,1);
v(:,2) = v(:,2) + outlines_min(1,2);
lattice = [];
for iE=1:size(e,1)
    lattice = [lattice; v(e(iE,1),:); v(e(iE,2),:); [NaN NaN] ];
end
for u = 2:3
    eSeamExt = o(u).eSeamExtRot;
    eSeamExt(:,[1 3]) = eSeamExt(:,[1 3]) + outlines_min(1,1);
    eSeamExt(:,[2 4]) = eSeamExt(:,[2 4]) + outlines_min(1,2);
    for iS=1:size(eSeamExt,1)
        lattice = [lattice; eSeamExt(iS,1:2); eSeamExt(iS,3:4); [NaN NaN] ];
    end
end

topPanel = polybuffer(lattice, 'lines', sWidth, 'JointType', 'square');
