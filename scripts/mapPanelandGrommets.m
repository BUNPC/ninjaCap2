function [sideOutline, grommets, topIdx, ear_slit_or_Cz_pos] = mapPanelandGrommets(vHead, fHead,grommets, springList, refpts, refpts_idx, panel3DVertices, refPts_neighbours, panel_name,outline_idx_to_include)

make_video =false;

% get points on outline seam that are closer to refPts
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

% If outline_idx_to_include parameter is passed then add those points to seam outline
% points. Here we used it to include corner points on top seam outline
if nargin > 9
    [outline_refpts_idx, sort_idx] = sort(outline_refpts_idx);
    refpts_outline = refpts_outline(sort_idx,:);
    refpts_outline_label = refpts_outline_label(sort_idx);
    idx_to_add = [];
    for  u = 1:length(outline_idx_to_include)
        pt = panel3DVertices(outline_idx_to_include(u),:);
        dist = sqrt(sum((panel3DVertices_01mm-pt).^2,2));
        [~,idx] = min(dist);
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

% create vertices and edges list for spring relaxation
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
% add connections along outline seam points and between seam and refPts
for u = 1:length(outline_refpts_idx)
    vHex = [vHex; panel3DVertices_01mm(outline_refpts_idx(u),:)];
    idx = find(ismember(refpts.labels,refpts_outline_label(u))==1);
    pos = refpts.pos(idx,:);
    idx = find(ismember(vHex,pos,'rows')==1);
    if ~isempty(idx)
        eHex = [eHex; [idx size(vHex,1)]];
        if u > 1
            eHex = [eHex; [size(vHex,1)-1 size(vHex,1)]];
        end
    end
end

% for top panel adding edge between last and first point
if strcmp(panel_name,'top')
    idx = find(ismember(vHex,panel3DVertices_01mm(outline_refpts_idx(1),:),'rows')==1);
    eHex = [eHex; [size(vHex,1) idx]];
end
vHex_outline_idx = vHex_refpts_length+1:size(vHex,1);

% Add cross connections between seam outline refPts. Currently it cross
% connections are if distance between refPt and outline is greater than
% 5mm. But I think it also should work if ignore distance condition.
for  u = 1:length(outline_refpts_idx)-1
    out_idx1 = find(ismember(vHex,refpts_outline(u,:),'rows')==1);
    out_idx2 = find(ismember(vHex,refpts_outline(u+1,:),'rows')==1);
    ref_idx1 = find(ismember(refpts.labels,refpts_outline_label(u))==1);
    pos = refpts.pos(ref_idx1,:);
    ref_idx1 = find(ismember(vHex,pos,'rows')==1);
    ref_idx2 = find(ismember(refpts.labels,refpts_outline_label(u+1))==1); 
    pos = refpts.pos(ref_idx2,:);
    ref_idx2 = find(ismember(vHex,pos,'rows')==1);

%     if  sqrt(sum((vHex(out_idx1,:)-vHex(ref_idx1,:)).^2,2)) >5
    if ~isempty(out_idx2) & ~isempty(ref_idx1)
        eHex = [eHex; [ref_idx1 out_idx2]];
    end
%     end
%     if  sqrt(sum((vHex(out_idx2,:)-vHex(ref_idx2,:)).^2,2)) >5
    if ~isempty(out_idx1) & ~isempty(ref_idx2)
        eHex = [eHex; [ref_idx2 out_idx1]];
    end
%     end
end

%%

% add connection between grommets based on springList and also connect
% grommet to nearest two refPts.
ref_length = size(vHex,1);
gidx = [];
for u = 1:length(grommets)
 if strcmp(grommets(u).panel,panel_name) %&& grommets(u).optType ~= 3
     gidx = [gidx; u];
 end
end

gvHex = [];
optType = [];
geHex = [];
gvOrder = [];
for u = 1:length(gidx)
    ii = find((springList(:,1) == gidx(u) | springList(:,2) == gidx(u)));
    connections = springList(ii,1:2);
    connections = [gidx(u); setdiff(connections(:),gidx(u))];
    for v = 1:length(connections)
        if strcmp(grommets(connections(v)).panel,panel_name)
            if isempty(gvHex)
                gvHex = grommets(connections(v)).posHead;
                optType = grommets(connections(v)).optType;
                idx = 1;
                gvOrder = [gvOrder; connections(v)];
            else
                idx = find(ismember(gvHex,grommets(connections(v)).posHead,'rows')==1);
                if isempty(idx)
                    gvHex = [gvHex; grommets(connections(v)).posHead];
                    optType = [optType; grommets(connections(v)).optType];
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

gIdx = size(vHex,1);
sIdx = find(optType ==1)+gIdx;
dIdx = find(optType ==2)+gIdx;
g_idx = size(vHex,1)+1;
vHex = [vHex; gvHex];
geHex = geHex+ref_length;
eHex = [eHex; geHex];
vHex_grommets_idx = vHex_refpts_length+length(vHex_outline_idx)+1:size(vHex,1);
ref_elength = size(eHex,1)-size(geHex,1);
for u = ref_length+1:size(vHex,1)
    dist = sqrt(sum((vHex(1:ref_length,:)-vHex(u,:)).^2,2));
    [~,sort_idx] = sort(dist);
    for v = 1:4
        eHex = [eHex; [u sort_idx(v)]];
    end
end

%%
max_h = max(vHex(:,2));


% vHex(:,2) = -(vHex(:,2)-max_h);
% temp = vHex(:,2);
% vHex(:,2) = vHex(:,3);
% vHex(:,3) = temp;


% if strcmp(panel_name,'top')
%     max_x = max(vHex(:,1));
%     vHex(:,1) = -(vHex(:,1)-max_x);
% end

% get euclidean distnace between connections.
hHex_temp  = sqrt(sum((vHex(eHex(:,1),:)-vHex(eHex(:,2),:)).^2,2));
hHex = hHex_temp;
%% Length Deformation to fix posterior issue
n_firstorder_egdes = size(eHex,1);

% Find bottom row edges
[min_z ] = min(vHex(:,3));
vidx_bottom_row = find(vHex(:,3) < min_z+10);
edges_last_row = ismember(eHex(:,1),vidx_bottom_row) & ismember(eHex(:,2),vidx_bottom_row);
edges_last_row = find(edges_last_row ==1);

% second order sprinng connections
for u = 1:size(vHex,1)
   vert_distnaces = sqrt(sum((vHex-vHex(u,:)).^2,2));
   id = find(vert_distnaces <= 40);
   for v = 1:length(id)
       if (sum(ismember(eHex,[u id(v)],'rows')) == 0) & sum(ismember(eHex,[id(v) u],'rows')) == 0
           if u ~= id(v)
            eHex = [eHex; [u id(v)]];
           end
       end
   end
end

hHex = sqrt(sum((vHex(eHex(:,1),:)-vHex(eHex(:,2),:)).^2,2));
scale_factor = ones(size(eHex,1),1);
scale_factor(n_firstorder_egdes+1:end) = 0.5;
pTopInt = [160 140];

hHex_deform = hHex;
count = 1;
rms_error = [];
% max_dLen = [];
% max_idx = [];

while count <= 200
    %count
    vHex2 = SpringRelax_func_scaled( vHex, eHex, hHex_deform, scale_factor );
    for ii=1:2400
        vHex2 = SpringRelax_func_scaled( vHex2, eHex, hHex_deform,scale_factor );

        lst = find(abs(vHex2(:,3))>0.1);
        vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;

        lst = find(vHex2(:,3)<0.5);
        vec = vHex2(lst,1:2)-ones(length(lst),1)*pTopInt;
        vec = vec ./ sum( vec.^2, 2).^0.5;
        vHex2(lst,1:2) = vHex2(lst,1:2) + 0.1 * vec;    
    end

    for ii=1:1800
        vHex2 = SpringRelax_func_scaled( vHex2, eHex, hHex_deform,scale_factor);

        lst = find(abs(vHex2(:,3))>0.1);
        vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;
    end
    
    eHexLen = sum((vHex2(eHex(:,1),:)-vHex2(eHex(:,2),:)).^2,2).^0.5;
    dLen = eHexLen-hHex;
    dLen_last_row =  dLen(edges_last_row);
    dLen = dLen(1:n_firstorder_egdes);
    rms_error(count) = sqrt(sum(dLen.^2)/length(dLen));
    
    [max_dLen, max_idx] = max(dLen_last_row);
    max_idx = edges_last_row(max_idx);
    display([count rms_error(count) max_dLen max_idx])
    hHex_deform( max_idx) = hHex_deform(max_idx)-0.4*dLen(max_idx);
%     display([count rms_error(count) max(abs(dLen))])
%     hHex_deform(1:n_firstorder_egdes) = hHex_deform(1:n_firstorder_egdes)-0.1*dLen;
%     if count==1
%         dLen_1 = dLen;
%     end
    count = count+1;
%     save([panel_name filesep 'iter_' num2str(count) '.mat'],'vHex2','rms_error','dLen')
%     figure(10); plot_panel(vHex2,eHex(1:n_firstorder_egdes,:), dLen)
%     figure(100);hist([dLen dLen_1],[-1:0.025:1]*1.5)
%     drawnow
end

%%
% % get deodesic distnaces between connections
% hHex = getGeodesicdist(vHead, fHead, vHex, eHex);
% 
% hHex_diff = hHex-hHex_temp;
% nidx = find(hHex_diff < 0);
% hHex(nidx) = hHex_temp(nidx);
% gidx = find(hHex_diff >= 4);
% hHex(gidx) = hHex_temp(gidx);
% zero_dist_idx = find(hHex == 0);
% if ~isempty(zero_dist_idx)
%     hHex(zero_dist_idx) = 0.5;
% end

% save(['deform_correction_results_sDist' filesep 'sDist_' panel_name '.mat'],'hHex');

% load(['deform_correction_results_sDist' filesep 'sDist_' panel_name '.mat'])
 
% run spring relaxation and push cap to 2D plane

% save(['panel_' panel_name '.mat'],'vHex','eHex','hHex','g_idx');
%% old code
%         vHex2 = SpringRelax_func( vHex, eHex, hHex );
% 
%         pTopInt = [160 140];
%         if 1
%             ind = 1;
%             movieFileName = ['capMapping' panel_name '.mp4'];
%             figure(101)
%             vidfile = VideoWriter(movieFileName,'MPEG-4');
%             open(vidfile);
% 
%             custom_colormap = zeros(31,3);
%             temp = linspace(1,0,16);
%             custom_colormap(1:15,1) = temp(1:15); 
%             custom_colormap(2:16,2) = flip(temp(1:15)); 
%             custom_colormap(17:end,3) = flip(temp(1:15)); 
%             custom_colormap(17:end,2) = temp(2:16); 
%             color_var = linspace(-0.5,0.5,31);
%         end
% 
%         for ii=1:2400
%             vHex2 = SpringRelax_func( vHex2, eHex, hHex );
%         %     figure; plot3(vHex2(:,1),vHex2(:,2),vHex2(:,3),'.'); axis equal
% 
%             lst = find(abs(vHex2(:,3))>0.1);
%             vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;
% 
% 
%             lst = find(vHex2(:,3)<0.5);
%             vec = vHex2(lst,1:2)-ones(length(lst),1)*pTopInt;
%             vec = vec ./ sum( vec.^2, 2).^0.5;
%             vHex2(lst,1:2) = vHex2(lst,1:2) + 0.1 * vec;
%         %     figure; plot3(vHex2(:,1),vHex2(:,2),vHex2(:,3),'.'); axis equal
%             if make_video
%                 if ii == 1 || rem(ii,10) == 0
%                     %% make movie
%                     eHexLen = sum((vHex2(eHex(:,1),:)-vHex2(eHex(:,2),:)).^2,2).^0.5;
%                     dLen = eHexLen-hHex;
%                 %     plot3( vHex(:,1), vHex(:,2), vHex(:,3), 'k.');
% 
%                     plot3( vHex2(1:gIdx,1), vHex2(1:gIdx,2), vHex2(1:gIdx,3), 'k.');
%                     hold on
%                     plot3( vHex2(sIdx,1), vHex2(sIdx,2), vHex2(sIdx,3), 'r.', 'MarkerSize',12);
%                     hold on
%                     plot3( vHex2(dIdx,1), vHex2(dIdx,2), vHex2(dIdx,3), 'b.', 'MarkerSize',12);
%                 %     for u = 1:ref_elength
%                 %         line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','b','LineWidth',0.5)
%                 %     end
%                 %     for u = ref_elength+1:size(eHex,1)
%                 %         line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','g','LineWidth',0.5)
%                 %     end
%                     for u = 1:ref_elength
%                         hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color',[255 255 102]/255.0,'LineWidth',0.5);
%                         [~,color_idx] = min(abs(color_var-dLen(u)));
%                         set(hl,'color',custom_colormap(color_idx,:));
%         %                 if dLen(u)>3
%         %                     set(hl,'color','r');
%         %                 elseif dLen(u)<-3
%         %                     set(hl,'color','c');
%         %                 end
%                     end
%                     for u =  ref_elength+1:size(eHex,1)
%                         hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color',[102 255 102]/255.0,'LineWidth',0.5);
%                         [~,color_idx] = min(abs(color_var-dLen(u)));
%                         set(hl,'color',custom_colormap(color_idx,:));
%         %                 if dLen(u)>3
%         %                     set(hl,'color','r');
%         %                 elseif dLen(u)<-3
%         %                     set(hl,'color','c');
%         %                 end
%                     end
%                     hold off
%                     axis equal
%         %             colorbar;
%                     drawnow
%                     F = getframe(gcf);
%                     writeVideo(vidfile,F);
%                     ind = ind+1;
%                 end
%             end
%         end
%         % if ~strcmp(panel_name,'top')
%         loop_count  = 1;
%         continue_loop = true;
%         while (continue_loop)
%             if loop_count == 1
%                 hHex_outline_adjust = hHex;
%             end
%             for ii=1:1800
%                 vHex2 = SpringRelax_func( vHex2, eHex, hHex_outline_adjust );
% 
%                 lst = find(abs(vHex2(:,3))>0.1);
%                 vHex2(lst,3) = vHex2(lst,3) - sign(vHex2(lst,3))*0.1;
%                 if ii == 1800 
%                     if ii == 1 || rem(ii,10) == 0
%                         %% make movie
%                         eHexLen = sum((vHex2(eHex(:,1),:)-vHex2(eHex(:,2),:)).^2,2).^0.5;
%                         dLen = eHexLen-hHex;
%                     %     plot3( vHex(:,1), vHex(:,2), vHex(:,3), 'k.');
% 
%                         plot3( vHex2(1:gIdx,1), vHex2(1:gIdx,2), vHex2(1:gIdx,3), 'k.');
%                         hold on
%                         plot3( vHex2(sIdx,1), vHex2(sIdx,2), vHex2(sIdx,3), 'r.', 'MarkerSize',12);
%                         hold on
%                         plot3( vHex2(dIdx,1), vHex2(dIdx,2), vHex2(dIdx,3), 'b.', 'MarkerSize',12);
%                     %     for u = 1:ref_elength
%                     %         line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','b','LineWidth',0.5)
%                     %     end
%                     %     for u = ref_elength+1:size(eHex,1)
%                     %         line([vHex(eHex(u,1),1) vHex(eHex(u,2),1)], [vHex(eHex(u,1),2) vHex(eHex(u,2),2)], [vHex(eHex(u,1),3) vHex(eHex(u,2),3)],'color','g','LineWidth',0.5)
%                     %     end
%                         for u = 1:ref_elength
%                             hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color',[255 255 102]/255.0,'LineWidth',0.5);
%                             [~,color_idx] = min(abs(color_var-dLen(u)));
%                             set(hl,'color',custom_colormap(color_idx,:));
%             %                 if dLen(u)>3
%             %                     set(hl,'color','r');
%             %                 elseif dLen(u)<-3
%             %                     set(hl,'color','c');
%             %                 end
%                         end
%                         for u =  ref_elength+1:size(eHex,1)
%                             hl=line([vHex2(eHex(u,1),1) vHex2(eHex(u,2),1)], [vHex2(eHex(u,1),2) vHex2(eHex(u,2),2)], [vHex2(eHex(u,1),3) vHex2(eHex(u,2),3)],'color',[102 255 102]/255.0,'LineWidth',0.5);
%                             [~,color_idx] = min(abs(color_var-dLen(u)));
%                             set(hl,'color',custom_colormap(color_idx,:));
%             %                 if dLen(u)>3
%             %                     set(hl,'color','r');
%             %                 elseif dLen(u)<-3
%             %                     set(hl,'color','c');
%             %                 end
%                         end
%                         hold off
%                         axis equal
%             %             colorbar;
%                         drawnow
%         %                 saveas(gcf,['deform_correction_results_sDist' filesep 'final_frame' panel_name '_' num2str(loop_count) '.fig'])
%                         F = getframe(gcf);
%                         writeVideo(vidfile,F);
%                         ind = ind+1;
%         %                 save(['deform_correction_results_sDist' filesep 'error_info_' panel_name '_' num2str(loop_count) '.mat'],'vHex2','eHex','hHex','eHexLen','vHex_outline_idx');
%                         if 1 %strcmp(panel_name,'top')
%                             continue_loop = false;
%                         else
%         %                 [ans,hHex_outline_adjust] = repeat_loop(vHex2,eHex,hHex,hHex_outline_adjust,eHexLen,vHex_outline_idx);
%         %                 if ~ans || loop_count > 5
%         %                     continue_loop = false; 
%         %                 end
%                         loop_count = loop_count+1;
%                         end
%                     end
%                 end
%             end
%         end
%         % end
% 
%         if make_video
%             close(vidfile)
%         end
%%
eHexLen = sum((vHex2(eHex(:,1),:)-vHex2(eHex(:,2),:)).^2,2).^0.5;
mean( (eHexLen-hHex).^2 ).^0.5;
dLen = eHexLen-hHex;
figure
hist( dLen, [-10:1:10] )

% temp = vHex2(:,2);
% vHex2(:,2) = vHex2(:,3);
% vHex2(:,3) = temp;

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

% retrieve 2D seam outline points and also LPA and RPA positions for side
% panels for positioning ear slit.
sideOutline = vHex2(vHex_outline_idx,:);
if strcmp(panel_name,'sideLeft') || strcmp(panel_name,'sideRight')
    sideOutline  = flipud(sideOutline);
end
if strcmp(panel_name,'sideLeft')
    idx = find(ismember(refpts.labels,'LPA')==1);
    LPA_pos = refpts.pos(idx,:);
%     LPA_pos = [LPA_pos(1) LPA_pos(3) -(LPA_pos(2)-max_h)];
    vHex_idx = find(ismember(vHex, LPA_pos ,'rows')==1);
    ear_slit_or_Cz_pos = vHex2(vHex_idx ,:);
elseif strcmp(panel_name,'sideRight')
    idx = find(ismember(refpts.labels,'RPA')==1);
    RPA_pos = refpts.pos(idx,:);
%     RPA_pos = [RPA_pos(1) RPA_pos(3) -(RPA_pos(2)-max_h)];
    vHex_idx = find(ismember(vHex, RPA_pos ,'rows')==1);
    ear_slit_or_Cz_pos = vHex2(vHex_idx ,:);
else
    idx = find(ismember(refpts.labels,'Cz')==1);
    LPA_pos = refpts.pos(idx,:);
%     LPA_pos = [-(LPA_pos(1)-max_x) LPA_pos(3) -(LPA_pos(2)-max_h)];
    vHex_idx = find(ismember(vHex, LPA_pos ,'rows')==1);
    ear_slit_or_Cz_pos = vHex2(vHex_idx ,:);
end

side3Doutline = vHex(vHex_outline_idx,:);
figure; plot3(side3Doutline(:,1), side3Doutline(:,2),side3Doutline(:,3));
hold on;
a = cat(1, grommets(arrayfun(@(x) strcmp(x.panel,  panel_name) && ~strcmp(x.type, '#NONE'), grommets)).posHead);
if ~isempty(a)
    scatter3(a(:, 1), a(:, 2), a(:, 3));
end
hold off;
axis equal
% assign grommet 2 panel position
for u = 1:length(gvOrder)
    grommets(gvOrder(u)).posPanel = vHex2(vHex_grommets_idx(u),1:2);
end





