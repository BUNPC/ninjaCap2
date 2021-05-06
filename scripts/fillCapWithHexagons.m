function o = fillCapWithHexagons(sideLeftOutline, sideLeftIdx, sideRightOutline, sideRightIdx, topOutline, topIdx)

o(1).v = [50 50]; % first vertex location
o(1).xoutline = topOutline(:,1) - min(topOutline(:,1));
o(1).youtline = topOutline(:,2) - min(topOutline(:,2));

o(2).v = [100 100]; % first vertex location
o(2).xoutline = sideLeftOutline(:,1) - min(sideLeftOutline(:,1));
o(2).youtline = sideLeftOutline(:,2) - min(sideLeftOutline(:,2));

o(3).v = [100 100]; % first vertex location
o(3).xoutline = sideRightOutline(:,1) - min(sideRightOutline(:,1));
o(3).youtline = sideRightOutline(:,2) - min(sideRightOutline(:,2));


hEdge = 10; % hexagon edge length


nOutlines = length(o);

%%
% set up the seams and get the length along the left and right seams

iSeam = 1;
% o(iSeam).seamIs = [1 0 1 0 0];
xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamLength = zeros(length(topOutline),1);
o(iSeam).seamIs = zeros(length(topOutline),1);
o(iSeam).seamIs(topIdx(1)) = 1;
for ii = topIdx(1)+1:topIdx(2)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end
o(iSeam).seamIs(topIdx(3)) = 1;
for ii = topIdx(3)+1:topIdx(4)
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end


iSeam = 2;
o(iSeam).seamStartIdx = sideLeftIdx(1); % I am pretty sure that start at back of head and move along the top of the head to the front
o(iSeam).seamEndIdx = sideLeftIdx(2); %19
o(iSeam).seamLength = zeros(length(sideLeftOutline),1);
o(iSeam).seamIs = zeros(length(sideLeftOutline),1);

xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end

iSeam = 3;
o(iSeam).seamStartIdx = sideRightIdx(1); % I am pretty sure that start at back of head and move along the top of the head to the front
o(iSeam).seamEndIdx = sideRightIdx(2);%20;
o(iSeam).seamLength = zeros(length(sideRightOutline),1);
o(iSeam).seamIs = zeros(length(sideRightOutline),1);

xx = o(iSeam).xoutline;
yy = o(iSeam).youtline;
o(iSeam).seamIs(o(iSeam).seamStartIdx) = 1;
for ii = o(iSeam).seamStartIdx+1 : o(iSeam).seamEndIdx
    o(iSeam).seamLength(ii) = o(iSeam).seamLength(ii-1) + norm( [xx(ii)-xx(ii-1) yy(ii)-yy(ii-1)] );
    o(iSeam).seamIs(ii) = 1;
end


%%
% create the mask from the outline of the cap
% figure(1)

for iO = 1:nOutlines
    
    [xgrid,ygrid] = meshgrid(1:max(o(iO).xoutline)*1.05,1:max(o(iO).youtline)*1.05);

    o(iO).Imask = inpolygon(xgrid,ygrid,o(iO).xoutline,o(iO).youtline);

%     subplot(1,nOutlines,iO)
%     imagesc(o(iO).Imask)
%     axis image
%     set(gca,'ydir',"reverse")
    
end


%%
% fill mask with hexagons
% figure(1)

for iO = 1:nOutlines

%     subplot(1,nOutlines,iO)
    [o(iO).v, o(iO).e, o(iO).vOut, o(iO).eOut] = fillCapWithHexagons_func( o(iO).v, hEdge, o(iO).Imask );

end

% save temp.mat o

%%

% for outline edge (not outline seam)
% pull outside vertices to the outline

for iO = 1:length(o)
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    for iE = 1:size(o(iO).eOut,1)
        iE
        xy2 = [o(iO).vOut(o(iO).eOut(iE,1),:) o(iO).vOut(o(iO).eOut(iE,2),:)];
        out = lineSegmentIntersect(xy1, xy2);
        idx = find(out.intAdjacencyMatrix==1);
        if o(iO).seamIs(idx)==0
            o(iO).v(end+1,:) = o(iO).vOut(o(iO).eOut(iE,1),:);
            o(iO).v(end+1,:) = [out.intMatrixX(idx(1)) out.intMatrixY(idx(1))]; %vOut(eOut(iE,2),:);
            o(iO).e(end+1,:) = [size(o(iO).v,1)-1 size(o(iO).v,1)];
        end
    end
end

%%
% Work on struts crossing the seam
for iO = 2:3
    iO
    v = o(iO).v;
    vOut = o(iO).vOut;
    e = o(iO).e;
    eOut = o(iO).eOut;
    
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    eSeamX = zeros(size(eOut,1),1);
    eSeamXidx = zeros(size(eOut,1),1);
    eSeamXpts = zeros(size(eOut,1),6);
    eSeamXlen = zeros(size(eOut,1),1);
    eSeamXtheta = zeros(size(eOut,1),1);
    eSeamXrot = zeros(size(eOut,1),2,2);
    eSeamXptsRot = zeros(size(eOut,1),6);
    
    for iE = 1:size(eOut,1)
        iE
        idx = []; foo=0;
        while isempty(idx) % this was thrown in to deal with an outward edge that didn't cross the outline. This happens because of the discrete mask I use prior to this
            xy2 = [vOut(eOut(iE,1),:)+foo*(vOut(eOut(iE,1),:)-vOut(eOut(iE,2),:)) vOut(eOut(iE,2),:)+foo*(vOut(eOut(iE,2),:)-vOut(eOut(iE,1),:))];
            foo = foo + 0.1;
            out = lineSegmentIntersect(xy1, xy2);
            idx = find(out.intAdjacencyMatrix==1);
        end
        idx = idx(1)
        eSeamXidx(iE) = idx;
        if o(iO).seamIs(idx)==1
            eSeamX(iE) = 1;
            eSeamXpts(iE,1:4) = [vOut(eOut(iE,1),:) out.intMatrixX(idx) out.intMatrixY(idx)];
            eSeamXlen(iE) = o(iO).seamLength(idx)*(1-out.intNormalizedDistance1To2(idx)) + o(iO).seamLength(idx+1)*out.intNormalizedDistance1To2(idx);
            
            % get rotation to top panel
            theta = asin( (xy1(idx,3)-xy1(idx,1)) / norm(xy1(idx,3:4) - xy1(idx,1:2)) );
            eSeamXtheta(iE) = theta;
            if iO==2 % left side
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                pt = get_point_on_outline([(o(1).xoutline(1:topIdx(2))) (o(1).youtline(1:topIdx(2)))], eSeamXlen(iE));
                eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + pt';
                eSeamXptsRot(iE,3:4) = pt';
%                 eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + [0 eSeamXlen(iE)]';
%                 eSeamXptsRot(iE,3:4) = [0 eSeamXlen(iE)]';
            else % right side
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                pt = get_point_on_outline([flipud(o(1).xoutline(topIdx(3):topIdx(4))) flipud(o(1).youtline(topIdx(3):topIdx(4)))], eSeamXlen(iE));
                eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + pt';
                eSeamXptsRot(iE,3:4) = pt';
%                 eSeamXptsRot(iE,1:2) = R * (eSeamXpts(iE,1:2)-eSeamXpts(iE,3:4))' + [75 eSeamXlen(iE)]';
%                 eSeamXptsRot(iE,3:4) = [max(o(1).xoutline) eSeamXlen(iE)]';
            end
            eSeamXrot(iE,:,:) = R;
        end
    end

    o(iO).eSeamX = eSeamX;
    o(iO).eSeamXidx = eSeamXidx;
    o(iO).eSeamXpts = eSeamXpts;
    o(iO).eSeamXlen = eSeamXlen;
    o(iO).eSeamXtheta = eSeamXtheta;
    o(iO).eSeamXrot = eSeamXrot;
    o(iO).eSeamXptsRot = eSeamXptsRot;
end

% move side panel extended edges to nearest top panel vertex

for iO = 2:3
    for iE = 1:length(o(iO).eSeamX)
        if o(iO).eSeamX(iE)==1
            
            % find nearest vertex inside the top panel
            p1 = o(iO).eSeamXptsRot(iE,[3:4]);
            p2lst = o(1).vOut(o(1).eOut(:,1),:);
            rho = sum((ones(size(p2lst,1),1)*p1 - p2lst).^2,2).^0.5;
            [foo,idx] = min(rho);
            
            % move side panel extended edge to the top panel vertex
            o(iO).eSeamXptsRot(iE,[5:6]) = p2lst(idx,:);
            
            % rotate back into side panel space
            p0 = p2lst(idx,:) - p1;
            theta = o(iO).eSeamXtheta(iE);
            R = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            o(iO).eSeamXpts(iE,[5:6]) = R * p0' + o(iO).eSeamXpts(iE,3:4)';
        end
    end
    
    lst = find(o(iO).eSeamX==1);
    o(iO).eSeamExt = o(iO).eSeamXpts(lst,[1 2 5 6]);
    o(iO).eSeamExtRot = o(iO).eSeamXptsRot(lst,[1 2 5 6]);

end

%%
% Remove duplicates nodes in v, e, vout and eout
for iO = 1:3
        v = o(iO).v;
        e = o(iO).e;
        % clean vertices and edges
        [u,I,J] = unique(v, 'rows', 'first');
        duplicated_rows = setdiff((1:size(v,1))',I);
        for uu = 1:length(duplicated_rows)
            idxs = find(v(:,1) == v(duplicated_rows(uu),1) & v(:,2) == v(duplicated_rows(uu),2));
            idx = setdiff(idxs,duplicated_rows);
            e(e==duplicated_rows(uu)) = idx;
        end

        nNodes = size(v,1);
        map = (1:nNodes)';
        map(duplicated_rows) = [];
        mapTemp = (1:length(map))';
        nodeMap = zeros(nNodes,1);
        nodeMap(map) = mapTemp;

        edgesNew = nodeMap(e);
        nodesNew = v;
        nodesNew(duplicated_rows,:) = [];

        zero_idx = find(edgesNew(:,1) == 0 | edgesNew(:,2)==0);
        edgesNew(zero_idx,:) = [];
 
        o(iO).v= nodesNew;
        edgesNew = sort(edgesNew,2);
        [u,I,J] = unique(edgesNew, 'rows', 'first');
        duplicated_rows = setdiff((1:size(edgesNew,1))',I);
        edgesNew(duplicated_rows,:)= [];
        o(iO).e = edgesNew;
end
% for iO = 1:3
%         v = o(iO).vOut;
%         e = o(iO).eOut;
%         [u,I,J] = unique(e(:,1), 'last');
%         duplicated_rows = setdiff((1:size(e(:,1),1))',I);
%         e(duplicated_rows,:) = [];
%         o(iO).eOut = e;
% end
%%
for iO = 2:3
    v = o(iO).v;
    e = o(iO).e;
%     left_extensions = [];
%     leftExtension_lattice = [];
%     TleftExtension_lattice = [];
    for pt = 1:size(o(iO).eSeamExt,1) 
    %     if o(iO).eSeamXpts(pt,1) ~= 0 && o(iO).eSeamXpts(pt,2) ~= 0
    %         if o(iO).seamIs(pt)==1
                vidx = find(v(:,1) == o(iO).eSeamExt(pt,1) & v(:,2) == o(iO).eSeamExt(pt,2));
                eidx = find(e(:,1) == vidx | e(:,2) == vidx);
                
                 pts = [];
                 vidxs = e(eidx,:);
                 other_vidx = setdiff(vidxs(:),vidx);
                 for u =1:length(other_vidx)
                     pts = [pts v(vidx,:) v(other_vidx(u),:)];
                 end
                 o(iO).eExtensions(pt).pts = pts;
%                 if length(eidx) == 1
%                     other_vidx = setdiff(e(eidx,:),vidx);
%                     other_eidx = setdiff(find(e(:,1) == other_vidx | e(:,2) == other_vidx),eidx);
%                     other_vidxs = setdiff([e(other_eidx(1),:) e(other_eidx(2),:)], other_vidx);
%                     o(iO).eExtensions(pt).pts = [v(vidx,:) v(other_vidx,:) v(other_vidx,:) v(other_vidxs(1),:) v(other_vidx,:) v(other_vidxs(2),:)];
%                 elseif length(eidx) == 2
%                     other_vidxs = setdiff([e(eidx(1),:) e(eidx(2),:)], vidx);
%                     o(iO).eExtensions(pt).pts = [v(vidx,:) v(other_vidxs(1),:) v(vidx,:) v(other_vidxs(2),:)];
%                 end
    %             for iL = 1:4:length(left_extensions(pt).pts)
    %                 leftExtension_lattice = [leftExtension_lattice; [left_extensions(pt).pts(iL)+min(sideLeftOutline(:,1)) left_extensions(pt).pts(iL+1)+min(sideLeftOutline(:,2))];...
    %                     [left_extensions(pt).pts(iL+2)+min(sideLeftOutline(:,1)) left_extensions(pt).pts(iL+3)+min(sideLeftOutline(:,2))]; [NaN NaN]];
    %             end
                Movingpts = [o(iO).eSeamExt(pt,1) o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3) o(iO).eSeamExt(pt,4)];
                FixedPts = [o(iO).eSeamExtRot(pt,1) o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3) o(iO).eSeamExtRot(pt,4)];
                tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
                Xpts = o(iO).eExtensions(pt).pts(1:2:end);
                Ypts = o(iO).eExtensions(pt).pts(2:2:end);
                [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
                temppts = [Txpts;Typts];
                temppts = temppts;%+[min(topOutline(:,1)); min(topOutline(:,2))];
                o(iO).eExtensionsRot(pt).pts = temppts(:)';
    
    %             for iL = 1:4:length(Tleft_extensions(pt).pts)
    %                 TleftExtension_lattice = [TleftExtension_lattice; [Tleft_extensions(pt).pts(iL) Tleft_extensions(pt).pts(iL+1)];...
    %                     [Tleft_extensions(pt).pts(iL+2) Tleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
    %             end
    %         end
    %     end
    end
end

%%

iO = 3;
v = o(1).v;
e = o(1).e;
for iO = 2:3
% topright_extensions = [];
% toprightExtension_lattice = [];
% TtoprightExtension_lattice = [];
for pt = 1:size(o(iO).eSeamExtRot,1) 
%     if o(iO).eSeamX(pt)==1
%     if o(iO).eSeamXptsRot(pt,5) ~= 0 && o(iO).eSeamXptsRot(pt,6) ~= 0
         vidx = find(v(:,1) == o(iO).eSeamExtRot(pt,3) & v(:,2) == o(iO).eSeamExtRot(pt,4));
         eidx = find(e(:,1) == vidx | e(:,2) == vidx);
         
         pts = [];
         vidxs = e(eidx,:);
         other_vidx = setdiff(vidxs(:),vidx);
         for u =1:length(other_vidx)
             pts = [pts v(vidx,:) v(other_vidx(u),:)];
         end
         
%          if length(eidx) == 2
%              vidxs = e(eidx,:);
%              other_vidx = setdiff(vidxs(:),vidx);
%             pts = [v(vidx,:) v(other_vidx(1),:) v(vidx,:) v(other_vidx(2),:)];
%          elseif length(eidx) == 1
%              pts = [];
%              for ii = 1:length(eidx)
%                  other_vidx = setdiff(e(eidx(ii),:),vidx);
%     %              pts = [pts; v(vidx,:)  
%                  for jj = 1:length(other_vidx)
%                     pts = [pts v(vidx,:)  v(other_vidx(jj),:)];
%                     other_eidx = setdiff(find(e(:,1) == other_vidx(jj) | e(:,2) == other_vidx(jj)),eidx(ii));
%                     for kk = 1:length(other_eidx)
%                         other_vidx_2 = setdiff(e(other_eidx(kk),:),other_vidx(jj));
%                         for ll = 1:length(other_vidx_2 )
%                              pts = [pts v(other_vidx(jj),:)  v(other_vidx_2(ll),:)];
%                         end
%                     end
%                  end
%              end
%          end
         o(iO).eExtensionsTop(pt).pts  = pts;
%          for iL = 1:4:length(topright_extensions(pt).pts)
%             toprightExtension_lattice = [toprightExtension_lattice; [topright_extensions(pt).pts(iL) topright_extensions(pt).pts(iL+1)];...
%                 [topright_extensions(pt).pts(iL+2) topright_extensions(pt).pts(iL+3)]; [NaN NaN]];
%          end
         
        FixedPts = [o(iO).eSeamExt(pt,1) o(iO).eSeamExt(pt,2); o(iO).eSeamExt(pt,3) o(iO).eSeamExt(pt,4)];
        Movingpts = [o(iO).eSeamExtRot(pt,1) o(iO).eSeamExtRot(pt,2); o(iO).eSeamExtRot(pt,3) o(iO).eSeamExtRot(pt,4)];
        tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
        Xpts = o(iO).eExtensionsTop(pt).pts(1:2:end);
        Ypts = o(iO).eExtensionsTop(pt).pts(2:2:end);
        [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
        temppts = [Txpts;Typts];
        if iO == 2
            temppts = temppts;%+[min(sideLeftOutline(:,1)); min(sideLeftOutline(:,2))];
        end
        if iO == 3
            temppts = temppts;%+[min(sideRightOutline(:,1)); min(sideRightOutline(:,2))];
        end
        o(iO).eExtensionsTopRot(pt).pts = temppts(:)';
%         for iL = 1:4:length(Ttopright_extensions(pt).pts)
%             TtoprightExtension_lattice = [TtoprightExtension_lattice; [Ttopright_extensions(pt).pts(iL) Ttopright_extensions(pt).pts(iL+1)];...
%                 [Ttopright_extensions(pt).pts(iL+2) Ttopright_extensions(pt).pts(iL+3)]; [NaN NaN]];
%         end
%     end
%     end
end
end



