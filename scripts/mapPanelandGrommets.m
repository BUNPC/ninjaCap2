function [sideOutline, grommets, topIdx, ear_slit_pos] = mapPanelandGrommets(grommets, springList, refpts, refpts_idx, panel3DVertices, refPts_neighbours, panel_name,outline_idx_to_include)

panel3DVertices_01mm = getEquiDistPtsAlongOutline(panel3DVertices,0.1);
refpts_outline = [];
refpts_outline_label = {};
outline_refpts_idx = [];
iii = [];
for u =1:20:size(panel3DVertices_01mm,1)
    [min_dist, idx1] = min(sqrt(sum((refpts.pos(refpts_idx,:)-panel3DVertices_01mm(u,:)).^2,2)));
    if min_dist < 20
        closest_pt = refpts.pos(refpts_idx(idx1),:);
        [min_dist, idx2] = min(sqrt(sum((panel3DVertices_01mm-closest_pt).^2,2)));
        refpts_outline = [refpts_outline; panel3DVertices_01mm(idx2,:)];
        outline_refpts_idx = [outline_refpts_idx; idx2];
        refpts_outline_label{end+1} = refpts.labels{refpts_idx(idx1)};
        iii = [iii; idx1];
    end
end
[refpts_outline,ia,ic] = unique(refpts_outline,'rows','stable');
refpts_outline_label = refpts_outline_label(ia);
outline_refpts_idx = outline_refpts_idx(ia);
if nargin > 7
    [outline_refpts_idx, sort_idx] = sort(outline_refpts_idx);
    refpts_outline = refpts_outline(sort_idx,:);
    refpts_outline_label = refpts_outline_label(sort_idx);
    idx_to_add = [];
    for  u = 1:length(outline_idx_to_include)
        pt = panel3DVertices(outline_idx_to_include(u),:);
        dist = sqrt(sum((panel3DVertices_01mm-pt).^2,2));
        [~,idx] = min(dist);
%         pt = panel3DVertices_01mm(idx,:);
        idx_to_add = [idx_to_add; idx];
    end
    idx_to_add  = sort(idx_to_add );
    topIdx = [];
    for u = 1:length(idx_to_add)
        idxs = find(outline_refpts_idx < idx_to_add(u));
        idx_pos = length(idxs)+1;
        topIdx = [topIdx; idx_pos];
        outline_refpts_idx = [outline_refpts_idx(1:length(outline_refpts_idx) < idx_pos); idx_to_add(u); outline_refpts_idx(1:length(outline_refpts_idx) >= idx_pos)];
        refpts_outline = [refpts_outline(1:length(refpts_outline) < idx_pos,:); panel3DVertices_01mm(idx_to_add(u),:); refpts_outline(1:length(refpts_outline) >= idx_pos,:)];
        [min_dist, idx1] = min(sqrt(sum((refpts.pos(refpts_idx,:)-panel3DVertices_01mm(idx_to_add(u),:)).^2,2))); 
        refpts_outline_label = [refpts_outline_label(1:length(refpts_outline_label) < idx_pos) refpts.labels{refpts_idx(idx1)} refpts_outline_label(1:length(refpts_outline_label) >= idx_pos)];
    end
else
    topIdx = [];
end
%%
vHex = [];
eHex= [];
hHex = [];

for u = 1:size(refPts_neighbours,1)
    for v = 1:9
        if ~isempty(refPts_neighbours{u,v})
            ref_idx = find(ismember(refpts.labels, refPts_neighbours{u,v})==1);
            if isempty(ref_idx)
                [u v]
            end
            ref_pos = refpts.pos(ref_idx,:);
            if isempty(vHex)
                vHex = ref_pos;
                idx = 1;
            else
                idx = find(ismember(vHex,ref_pos,'rows')==1);
                if isempty(idx)
                    vHex = [vHex; ref_pos];
                    idx = size(vHex,1);
                end
            end
            if v == 1
                v_idx = idx;
            else
                if isempty(eHex)
                    eHex = [v_idx idx];
                else
                    if (sum(ismember(eHex,[v_idx idx],'rows')) == 0) && (sum(ismember(eHex,[idx v_idx],'rows')) == 0)
                        eHex = [eHex; [v_idx idx]];
                    end
                end
                
            end
        end
    end
end
vHex_refpts_length = size(vHex,1);

%%
for u = 1:length(outline_refpts_idx)
    vHex = [vHex; panel3DVertices_01mm(outline_refpts_idx(u),:)];
    idx = find(ismember(refpts.labels,refpts_outline_label(u))==1);
    pos = refpts.pos(idx,:);
    idx = find(ismember(vHex,pos,'rows')==1);
    eHex = [eHex; [idx size(vHex,1)]];
    if u > 1
        eHex = [eHex; [size(vHex,1)-1 size(vHex,1)]];
    end
end
% for top panel adding edge between last outline point and first outline
% point
if strcmp(panel_name,'top')
    idx = find(ismember(vHex,panel3DVertices_01mm(outline_refpts_idx(1),:),'rows')==1);
    eHex = [eHex; [size(vHex,1) idx]];
end
vHex_outline_idx = vHex_refpts_length+1:size(vHex,1);
for  u = 1:length(outline_refpts_idx)-1
    out_idx1 = find(ismember(vHex,refpts_outline(u,:),'rows')==1);
    out_idx2 = find(ismember(vHex,refpts_outline(u+1,:),'rows')==1);
    ref_idx1 = find(ismember(refpts.labels,refpts_outline_label(u))==1);
    pos = refpts.pos(ref_idx1,:);
    ref_idx1 = find(ismember(vHex,pos,'rows')==1);
    ref_idx2 = find(ismember(refpts.labels,refpts_outline_label(u+1))==1); 
    pos = refpts.pos(ref_idx2,:);
    ref_idx2 = find(ismember(vHex,pos,'rows')==1);

    if  sqrt(sum((vHex(out_idx1,:)-vHex(ref_idx1,:)).^2,2)) >5
        eHex = [eHex; [ref_idx1 out_idx2]];
    end
    if  sqrt(sum((vHex(out_idx2,:)-vHex(ref_idx2,:)).^2,2)) >5
        eHex = [eHex; [ref_idx2 out_idx1]];
    end
end

%%
ref_length = size(vHex,1);
gidx = [];
for u = 1:length(grommets)
 if strcmp(grommets(u).panel,panel_name)
     gidx = [gidx; u];
 end
end

gvHex = [];
geHex = [];
gvOrder = [];
for u = 1:length(gidx)
    ii = find(springList(:,1) == gidx(u) | springList(:,2) == gidx(u));
    connections = springList(ii,1:2);
    connections = [gidx(u); setdiff(connections(:),gidx(u))];
    for v = 1:length(connections)
        if strcmp(grommets(connections(v)).panel,panel_name)
            if isempty(gvHex)
                gvHex = grommets(connections(v)).posHead;
                idx = 1;
                gvOrder = [gvOrder; connections(v)];
            else
                idx = find(ismember(gvHex,grommets(connections(v)).posHead,'rows')==1);
                if isempty(idx)
                    gvHex = [gvHex; grommets(connections(v)).posHead];
                    idx = size(gvHex,1);
                    gvOrder = [gvOrder; connections(v)];
                end
            end
            if v == 1
                v_idx = idx;
            else
                if isempty(geHex)
                    geHex = [v_idx idx];
                else
                    if (sum(ismember(geHex,[v_idx idx],'rows')) == 0) && (sum(ismember(geHex,[idx v_idx],'rows')) == 0)
                        geHex = [geHex; [v_idx idx]];
                    end
                end
            end
        end
    end
end

vHex = [vHex; gvHex];
geHex = geHex+ref_length;
eHex = [eHex; geHex];
vHex_grommets_idx = vHex_refpts_length+length(vHex_outline_idx)+1:size(vHex,1);
ref_elength = size(eHex,1)-size(geHex,1);
for u = ref_length+1:size(vHex,1)
    dist = sqrt(sum((vHex(1:ref_length,:)-vHex(u,:)).^2,2));
%     [~,min_idx] = min(dist);
    [~,sort_idx] = sort(dist);
    for v = 1:2
        eHex = [eHex; [u sort_idx(v)]];
    end
end

%%
hHex  = sqrt(sum((vHex(eHex(:,1),:)-vHex(eHex(:,2),:)).^2,2));
 
pTopInt = [160 140];
vHex2 = SpringRelax_func( vHex, eHex, hHex );

for ii=1:2400
    vHex2 = SpringRelax_func( vHex2, eHex, hHex );
    
    lst = find(abs(vHex2(:,3))>0.1);
    vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;
    
    lst = find(vHex2(:,3)<0.5);
    vec = vHex2(lst,1:2)-ones(length(lst),1)*pTopInt;
    vec = vec ./ sum( vec.^2, 2).^0.5;
    vHex2(lst,1:2) = vHex2(lst,1:2) + 0.1 * vec;
end

for ii=1:1800
    vHex2 = SpringRelax_func( vHex2, eHex, hHex );
    
    lst = find(abs(vHex2(:,3))>0.1);
    vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;
end

eHexLen = sum((vHex2(eHex(:,1),:)-vHex2(eHex(:,2),:)).^2,2).^0.5;
% mean( (eHexLen-hHex).^2 ).^0.5;
% dLen = eHexLen-hHex;
% figure
% hist( dLen, [-10:1:10] )

% figure
% plot3( vHex(:,1), vHex(:,2), vHex(:,3), 'k.');
% hold on
% plot3( vHex2(:,1), vHex2(:,2), vHex2(:,3), 'r.');
% for u = 1:ref_elength
%     line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','b','LineWidth',0.5)
% end
% for u = ref_elength+1:size(eHex,1)
%     line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','g','LineWidth',0.5)
% end
% for u = 1:ref_elength
%     hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color','k','LineWidth',0.5);
%     if dLen(u)>3
%         set(hl,'color','r');
%     elseif dLen(u)<-3
%         set(hl,'color','c');
%     end
% end
% for u =  ref_elength+1:size(eHex,1)
%     hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color','g','LineWidth',0.5);
%     if dLen(u)>3
%         set(hl,'color','r');
%     elseif dLen(u)<-3
%         set(hl,'color','c');
%     end
% end
% hold off
% axis equal

sideOutline = vHex2(vHex_outline_idx,:);
if strcmp(panel_name,'sideLeft') || strcmp(panel_name,'sideRight')
    sideOutline  = flipud(sideOutline);
end
if strcmp(panel_name,'sideLeft')
    idx = find(ismember(refpts.labels,'LPA')==1);
    LPA_pos = refpts.pos(idx,:);
    vHex_idx = find(ismember(vHex, LPA_pos ,'rows')==1);
    ear_slit_pos = vHex2(vHex_idx ,:);
elseif strcmp(panel_name,'sideRight')
    idx = find(ismember(refpts.labels,'RPA')==1);
    RPA_pos = refpts.pos(idx,:);
    vHex_idx = find(ismember(vHex, RPA_pos ,'rows')==1);
    ear_slit_pos = vHex2(vHex_idx ,:);
else
    ear_slit_pos = [];
end

for u = 1:length(gvOrder)
    grommets(gvOrder(u)).posPanel = vHex2(vHex_grommets_idx(u),1:2);
end
