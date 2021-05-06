function grommets = cutGrommets(grommets, topPanel, sideLeftPanel, sideRightPanel)

grommets = cutPanelGrommets(grommets, 'top', topPanel);
grommets = cutPanelGrommets(grommets, 'sideLeft', sideLeftPanel);
grommets = cutPanelGrommets(grommets, 'sideRight', sideRightPanel);

for u = 1:length(grommets)
    if isempty(grommets(u).panelIndex)
        grommets(u).panelIndex = 0;
    end
end

function grommets = cutPanelGrommets(grommets, panel_name, panel)
    if size(panel, 1) == 1
        for j = 1:length(grommets)
            if strcmp(grommets(j).panel, panel_name)
                grommets(j).panelIndex = 1;
            end
        end
        return;
    end
    
    for i = 1:size(panel, 1)
        % get bounding box
        p = union(panel(i, :));
        [x, y] = boundingbox(p);
        p = polyshape([x(1) y(1); x(2) y(1); x(2) y(2); x(1) y(2)]);
        
        % add buffer
        p = polybuffer(p, 0);
        
        for j = 1:length(grommets)
            if strcmp(grommets(j).panel, panel_name) && isempty(grommets(j).panelIndex) && isinterior(p, grommets(j).posPanel)
                grommets(j).panelIndex = i;
            end
        end
    end
end

end
