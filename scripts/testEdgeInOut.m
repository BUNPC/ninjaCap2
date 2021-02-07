% load template
% load temp
% load panels.mat
%%
testPlotPanels( o );

%%
% for outline edge (not outline seam)
% pull outside vertices to the outline

for iO = 1:length(o)
    xx = o(iO).xoutline;
    yy = o(iO).youtline;
    
    xy1 = [xx yy [xx(2:end); xx(1)] [yy(2:end); yy(1)] ];
    
    for iE = 1:size(o(iO).eOut,1)
        xy2 = [o(iO).vOut(o(iO).eOut(iE,1),:) o(iO).vOut(o(iO).eOut(iE,2),:)];
        out = lineSegmentIntersect(xy1, xy2);
        idx = find(out.intAdjacencyMatrix==1);
        if o(iO).seamIs(idx)==0
            o(iO).v(end+1,:) = o(iO).vOut(o(iO).eOut(iE,1),:);
            o(iO).v(end+1,:) = [out.intMatrixX(idx) out.intMatrixY(idx)]; %vOut(eOut(iE,2),:);
            o(iO).e(end+1,:) = [size(o(iO).v,1)-1 size(o(iO).v,1)];
        end
    end
end

%%
% replot
testPlotPanels( o );

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
% Work on struts crossing the seam
for iO = 2:3
    
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
        idx = []; foo=0;
        while isempty(idx) % this was thrown in to deal with an outward edge that didn't cross the outline. This happens because of the discrete mask I use prior to this
            xy2 = [vOut(eOut(iE,1),:) vOut(eOut(iE,2),:)+foo*(vOut(eOut(iE,2),:)-vOut(eOut(iE,1),:))];
            foo = foo + 0.1;
            out = lineSegmentIntersect(xy1, xy2);
            idx = find(out.intAdjacencyMatrix==1);
        end
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


%%
% replot
testPlotPanels( o );

%%
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
% replot
testPlotPanels( o );

%%
% Get second order extensions for left side panel
iO = 2;
v = o(iO).v;
e = o(iO).e;
left_extensions = [];
leftExtension_lattice = [];
TleftExtension_lattice = [];
for pt = 1:size(o(iO).eSeamXpts,1) 
    pt
    if o(iO).eSeamXpts(pt,1) ~= 0 && o(iO).eSeamXpts(pt,2) ~= 0
%         if o(iO).seamIs(pt)==1
            vidx = find(v(:,1) == o(iO).eSeamXpts(pt,1) & v(:,2) == o(iO).eSeamXpts(pt,2));
            eidx = find(e(:,1) == vidx | e(:,2) == vidx);
            if length(eidx) == 1
                other_vidx = setdiff(e(eidx,:),vidx);
                other_eidx = setdiff(find(e(:,1) == other_vidx | e(:,2) == other_vidx),eidx);
                other_vidxs = setdiff([e(other_eidx(1),:) e(other_eidx(2),:)], other_vidx);
                left_extensions(pt).pts = [v(vidx,:) v(other_vidx,:) v(other_vidx,:) v(other_vidxs(1),:) v(other_vidx,:) v(other_vidxs(2),:)];
            elseif length(eidx) == 2
                other_vidxs = setdiff([e(eidx(1),:) e(eidx(2),:)], vidx);
                left_extensions(pt).pts = [v(vidx,:) v(other_vidxs(1),:) v(vidx,:) v(other_vidxs(2),:)];
            end
            for iL = 1:4:length(left_extensions(pt).pts)
                leftExtension_lattice = [leftExtension_lattice; [left_extensions(pt).pts(iL)+min(sideLeftOutline(:,1)) left_extensions(pt).pts(iL+1)+min(sideLeftOutline(:,2))];...
                    [left_extensions(pt).pts(iL+2)+min(sideLeftOutline(:,1)) left_extensions(pt).pts(iL+3)+min(sideLeftOutline(:,2))]; [NaN NaN]];
            end
            Movingpts = [o(iO).eSeamXpts(pt,1) o(iO).eSeamXpts(pt,2); o(iO).eSeamXpts(pt,5) o(iO).eSeamXpts(pt,6)];
            FixedPts = [o(iO).eSeamXptsRot(pt,1) o(iO).eSeamXptsRot(pt,2); o(iO).eSeamXptsRot(pt,5) o(iO).eSeamXptsRot(pt,6)];
            tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
            Xpts = left_extensions(pt).pts(1:2:end);
            Ypts = left_extensions(pt).pts(2:2:end);
            [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
            temppts = [Txpts;Typts];
            temppts = temppts+[min(topOutline(:,1)); min(topOutline(:,2))];
            Tleft_extensions(pt).pts = temppts(:)';
            for iL = 1:4:length(Tleft_extensions(pt).pts)
                TleftExtension_lattice = [TleftExtension_lattice; [Tleft_extensions(pt).pts(iL) Tleft_extensions(pt).pts(iL+1)];...
                    [Tleft_extensions(pt).pts(iL+2) Tleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
%         end
    end
end
foo = polybuffer(leftExtension_lattice, 'lines', 1, 'JointType', 'square');
figure(44)
plot(foo)
axis image
% foo = polybuffer(TleftExtension_lattice, 'lines', 1, 'JointType', 'square');
% figure(55)
% plot(foo)
%%
% get second order extensions for right side panel
iO = 3;
v = o(iO).v;
e = o(iO).e;
right_extensions = [];
rightExtension_lattice = [];
TrightExtension_lattice = [];
for pt = 1:size(o(iO).eSeamXpts,1) 
    pt
    if o(iO).eSeamXpts(pt,1) ~= 0 && o(iO).eSeamXpts(pt,2) ~= 0
        vidx = find(v(:,1) == o(iO).eSeamXpts(pt,1) & v(:,2) == o(iO).eSeamXpts(pt,2));
        eidx = find(e(:,1) == vidx | e(:,2) == vidx);
        if length(eidx) == 1
            other_vidx = setdiff(e(eidx,:),vidx);
            other_eidx = setdiff(find(e(:,1) == other_vidx | e(:,2) == other_vidx),eidx);
            other_vidxs = setdiff([e(other_eidx(1),:) e(other_eidx(2),:)], other_vidx);
            right_extensions(pt).pts = [v(vidx,:) v(other_vidx,:) v(other_vidx,:) v(other_vidxs(1),:) v(other_vidx,:) v(other_vidxs(2),:)];
        elseif length(eidx) == 2
            other_vidxs = setdiff([e(eidx(1),:) e(eidx(2),:)], vidx);
            right_extensions(pt).pts = [v(vidx,:) v(other_vidxs(1),:) v(vidx,:) v(other_vidxs(2),:)];
        end
        for iL = 1:4:length(right_extensions(pt).pts)
            rightExtension_lattice = [rightExtension_lattice; [right_extensions(pt).pts(iL)+min(sideLeftOutline(:,1)) right_extensions(pt).pts(iL+1)+min(sideLeftOutline(:,2))];...
                [right_extensions(pt).pts(iL+2)+min(sideLeftOutline(:,1)) right_extensions(pt).pts(iL+3)+min(sideLeftOutline(:,2))]; [NaN NaN]];
        end
        Movingpts = [o(iO).eSeamXpts(pt,1) o(iO).eSeamXpts(pt,2); o(iO).eSeamXpts(pt,5) o(iO).eSeamXpts(pt,6)];
        FixedPts = [o(iO).eSeamXptsRot(pt,1) o(iO).eSeamXptsRot(pt,2); o(iO).eSeamXptsRot(pt,5) o(iO).eSeamXptsRot(pt,6)];
        tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
        Xpts = right_extensions(pt).pts(1:2:end);
        Ypts = right_extensions(pt).pts(2:2:end);
        [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
        temppts = [Txpts;Typts];
        temppts = temppts+[min(topOutline(:,1)); min(topOutline(:,2))];
        Tright_extensions(pt).pts = temppts(:)';
        for iL = 1:4:length(Tright_extensions(pt).pts)
            TrightExtension_lattice = [TrightExtension_lattice; [Tright_extensions(pt).pts(iL) Tright_extensions(pt).pts(iL+1)];...
                [Tright_extensions(pt).pts(iL+2) Tright_extensions(pt).pts(iL+3)]; [NaN NaN]];
        end
    end
end
foo = polybuffer(rightExtension_lattice, 'lines', 1, 'JointType', 'square');
figure(45)
plot(foo)
axis image
%%
% get second order extensions to top left side
iO = 2;
v = o(1).v;
e = o(1).e;

% % clean vertices and edges
% [u,I,J] = unique(v, 'rows', 'first');
% duplicated_rows = setdiff((1:size(v,1))',I);
% for uu = 1:length(duplicated_rows)
%     idxs = find(v(:,1) == v(duplicated_rows(uu),1) & v(:,2) == v(duplicated_rows(uu),2));
%     idx = setdiff(idxs,duplicated_rows);
%     e(e==duplicated_rows(uu)) = idx;
% end
% 
% nNodes = size(v,1);
% map = (1:nNodes)';
% map(duplicated_rows) = [];
% mapTemp = (1:length(map))';
% nodeMap = zeros(nNodes,1);
% nodeMap(map) = mapTemp;
% 
% edgesNew = nodeMap(e);
% nodesNew = v;
% nodesNew(duplicated_rows,:) = [];
% 
% zero_idx = find(edgesNew(:,1) == 0 | edgesNew(:,2)==0);
% edgesNew(zero_idx,:) = [];
% 
% v= nodesNew;
% e = edgesNew;


topleft_extensions = [];
topleftExtension_lattice = [];
TtopleftExtension_lattice = [];
for pt = 1:size(o(iO).eSeamXptsRot,1)
    if o(iO).eSeamX(pt)==1
        if o(iO).eSeamXptsRot(pt,5) ~= 0 && o(iO).eSeamXptsRot(pt,6) ~= 0
             vidx = find(v(:,1) == o(iO).eSeamXptsRot(pt,5) & v(:,2) == o(iO).eSeamXptsRot(pt,6));
             eidx = find(e(:,1) == vidx | e(:,2) == vidx);
             pts = [];
             for ii = 1:length(eidx)
                 other_vidx = setdiff(e(eidx(ii),:),vidx);
    %              pts = [pts; v(vidx,:)  
                 for jj = 1:length(other_vidx)
                    pts = [pts v(vidx,:)  v(other_vidx(jj),:)];
                    other_eidx = setdiff(find(e(:,1) == other_vidx(jj) | e(:,2) == other_vidx(jj)),eidx(ii));
                    for kk = 1:length(other_eidx)
                        other_vidx_2 = setdiff(e(other_eidx(kk),:),other_vidx(jj));
                        for ll = 1:length(other_vidx_2 )
                             pts = [pts v(other_vidx(jj),:)  v(other_vidx_2(ll),:)];
                        end
                    end
                 end
             end
             topleft_extensions(pt).pts  = pts;
             for iL = 1:4:length(topleft_extensions(pt).pts)
                topleftExtension_lattice = [topleftExtension_lattice; [topleft_extensions(pt).pts(iL) topleft_extensions(pt).pts(iL+1)];...
                    [topleft_extensions(pt).pts(iL+2) topleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
             end

            FixedPts = [o(iO).eSeamXpts(pt,1) o(iO).eSeamXpts(pt,2); o(iO).eSeamXpts(pt,5) o(iO).eSeamXpts(pt,6)];
            Movingpts = [o(iO).eSeamXptsRot(pt,1) o(iO).eSeamXptsRot(pt,2); o(iO).eSeamXptsRot(pt,5) o(iO).eSeamXptsRot(pt,6)];
            tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
            Xpts = topleft_extensions(pt).pts(1:2:end);
            Ypts = topleft_extensions(pt).pts(2:2:end);
            [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
            temppts = [Txpts;Typts];
            temppts = temppts+[min(sideLeftOutline(:,1)); min(sideLeftOutline(:,2))];
            Ttopleft_extensions(pt).pts = temppts(:)';
            for iL = 1:4:length(Ttopleft_extensions(pt).pts)
                TtopleftExtension_lattice = [TtopleftExtension_lattice; [Ttopleft_extensions(pt).pts(iL) Ttopleft_extensions(pt).pts(iL+1)];...
                    [Ttopleft_extensions(pt).pts(iL+2) Ttopleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
        end
    end
end

foo = polybuffer(TtopleftExtension_lattice, 'lines', 1, 'JointType', 'square');
figure(46)
plot(foo)
axis image

%%
%%
% get second order extensions to top right side
iO = 3;
v = o(1).v;
e = o(1).e;

% % clean vertices and edges
% [u,I,J] = unique(v, 'rows', 'first');
% duplicated_rows = setdiff((1:size(v,1))',I);
% for uu = 1:length(duplicated_rows)
%     idxs = find(v(:,1) == v(duplicated_rows(uu),1) & v(:,2) == v(duplicated_rows(uu),2));
%     idx = setdiff(idxs,duplicated_rows);
%     e(e==duplicated_rows(uu)) = idx;
% end
% 
% nNodes = size(v,1);
% map = (1:nNodes)';
% map(duplicated_rows) = [];
% mapTemp = (1:length(map))';
% nodeMap = zeros(nNodes,1);
% nodeMap(map) = mapTemp;
% 
% edgesNew = nodeMap(e);
% nodesNew = v;
% nodesNew(duplicated_rows,:) = [];
% 
% zero_idx = find(edgesNew(:,1) == 0 | edgesNew(:,2)==0);
% edgesNew(zero_idx,:) = [];
% 
% v= nodesNew;
% e = edgesNew;


topright_extensions = [];
toprightExtension_lattice = [];
TtoprightExtension_lattice = [];
for pt = 1:size(o(iO).eSeamXptsRot,1) 
    if o(iO).eSeamX(pt)==1
    if o(iO).eSeamXptsRot(pt,5) ~= 0 && o(iO).eSeamXptsRot(pt,6) ~= 0
         vidx = find(v(:,1) == o(iO).eSeamXptsRot(pt,5) & v(:,2) == o(iO).eSeamXptsRot(pt,6));
         eidx = find(e(:,1) == vidx | e(:,2) == vidx);
         if length(eidx) == 2
             vidxs = e(eidx,:);
             other_vidx = setdiff(vidxs(:),vidx);
            pts = [v(vidx,:) v(other_vidx(1),:) v(vidx,:) v(other_vidx(2),:)];
         elseif length(eidx) == 1
             pts = [];
             for ii = 1:length(eidx)
                 other_vidx = setdiff(e(eidx(ii),:),vidx);
    %              pts = [pts; v(vidx,:)  
                 for jj = 1:length(other_vidx)
                    pts = [pts v(vidx,:)  v(other_vidx(jj),:)];
                    other_eidx = setdiff(find(e(:,1) == other_vidx(jj) | e(:,2) == other_vidx(jj)),eidx(ii));
                    for kk = 1:length(other_eidx)
                        other_vidx_2 = setdiff(e(other_eidx(kk),:),other_vidx(jj));
                        for ll = 1:length(other_vidx_2 )
                             pts = [pts v(other_vidx(jj),:)  v(other_vidx_2(ll),:)];
                        end
                    end
                 end
             end
         end
         topright_extensions(pt).pts  = pts;
         for iL = 1:4:length(topright_extensions(pt).pts)
            toprightExtension_lattice = [toprightExtension_lattice; [topright_extensions(pt).pts(iL) topright_extensions(pt).pts(iL+1)];...
                [topright_extensions(pt).pts(iL+2) topright_extensions(pt).pts(iL+3)]; [NaN NaN]];
         end
         
        FixedPts = [o(iO).eSeamXpts(pt,1) o(iO).eSeamXpts(pt,2); o(iO).eSeamXpts(pt,5) o(iO).eSeamXpts(pt,6)];
        Movingpts = [o(iO).eSeamXptsRot(pt,1) o(iO).eSeamXptsRot(pt,2); o(iO).eSeamXptsRot(pt,5) o(iO).eSeamXptsRot(pt,6)];
        tform = fitgeotrans(Movingpts,FixedPts,'nonreflectivesimilarity');
        Xpts = topright_extensions(pt).pts(1:2:end);
        Ypts = topright_extensions(pt).pts(2:2:end);
        [Txpts,Typts] = transformPointsForward(tform, Xpts, Ypts);
        temppts = [Txpts;Typts];
        temppts = temppts+[min(sideRightOutline(:,1)); min(sideRightOutline(:,2))];
        Ttopright_extensions(pt).pts = temppts(:)';
        for iL = 1:4:length(Ttopright_extensions(pt).pts)
            TtoprightExtension_lattice = [TtoprightExtension_lattice; [Ttopright_extensions(pt).pts(iL) Ttopright_extensions(pt).pts(iL+1)];...
                [Ttopright_extensions(pt).pts(iL+2) Ttopright_extensions(pt).pts(iL+3)]; [NaN NaN]];
        end
    end
    end
end

foo = polybuffer(TtoprightExtension_lattice, 'lines', 1, 'JointType', 'square');
figure(48)
plot(foo)
axis image
%%
for pt = 15
%     pt = 10;
    TTtoprightExtension_lattice = [];
    for iL = 1:4:length(Ttopright_extensions(pt).pts)
        TTtoprightExtension_lattice = [TTtoprightExtension_lattice; [topright_extensions(pt).pts(iL) topright_extensions(pt).pts(iL+1)];...
        [topright_extensions(pt).pts(iL+2) topright_extensions(pt).pts(iL+3)]; [NaN NaN]];
    end
    foo = polybuffer(TTtoprightExtension_lattice, 'lines', 1, 'JointType', 'square');
    figure(45)
    plot(foo)
    axis image
    pause
end
%%
% create Side Left Panel

v = o(2).v;
e = o(2).e;
eSeamExt = o(2).eSeamExt;

v(:,1) = v(:,1) + min(sideLeftOutline(:,1));
v(:,2) = v(:,2) + min(sideLeftOutline(:,2));
eSeamExt(:,[1 3]) = eSeamExt(:,[1 3]) + min(sideLeftOutline(:,1));
eSeamExt(:,[2 4]) = eSeamExt(:,[2 4]) + min(sideLeftOutline(:,2));

lattice = [];
for iE=1:size(e,1)
    lattice = [lattice; v(e(iE,1),:); v(e(iE,2),:); [NaN NaN] ];
end
for iS=1:size(eSeamExt,1)
    lattice = [lattice; eSeamExt(iS,1:2); eSeamExt(iS,3:4); [NaN NaN] ];
end
stateLattice.sideLeftPanel = lattice;
for pt = 1:size(Ttopleft_extensions,2)
    if ~isempty(Ttopleft_extensions(pt).pts)
        pt
        for iL = 1:4:length(Ttopleft_extensions(pt).pts)
            if Ttopleft_extensions(pt).pts
            lattice = [lattice; [Ttopleft_extensions(pt).pts(iL) Ttopleft_extensions(pt).pts(iL+1)];...
                [Ttopleft_extensions(pt).pts(iL+2) Ttopleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
        end
    end
end


foo = polybuffer(lattice, 'lines', 1, 'JointType', 'square');

figure(22)
plot(foo)
hold on
plot(sideLeftPanel(1))
hold off
% foo = polybuffer(TtopleftExtension_lattice, 'lines', 1, 'JointType', 'square');
% hold on
% plot(foo)
% hold off
axis image
% save the STL file
% [f, v] = extrude(foo, 0.5);        
% saveModel('panelSideLeft.stl', f, v);


%%
% create Side Right Panel

v = o(3).v;
e = o(3).e;
eSeamExt = o(3).eSeamExt;

v(:,1) = v(:,1) + min(sideRightOutline(:,1));
v(:,2) = v(:,2) + min(sideRightOutline(:,2));
eSeamExt(:,[1 3]) = eSeamExt(:,[1 3]) + min(sideRightOutline(:,1));
eSeamExt(:,[2 4]) = eSeamExt(:,[2 4]) + min(sideRightOutline(:,2));

lattice = [];
for iE=1:size(e,1)
    lattice = [lattice; v(e(iE,1),:); v(e(iE,2),:); [NaN NaN] ];
end
for iS=1:size(eSeamExt,1)
    lattice = [lattice; eSeamExt(iS,1:2); eSeamExt(iS,3:4); [NaN NaN] ];
end

stateLattice.sideRightPanel = lattice;

for pt = 1:size(Ttopright_extensions,2)
    if ~isempty(Ttopright_extensions(pt).pts)
        pt
        for iL = 1:4:length(Ttopright_extensions(pt).pts)
            if Ttopright_extensions(pt).pts
            lattice = [lattice; [Ttopright_extensions(pt).pts(iL) Ttopright_extensions(pt).pts(iL+1)];...
                [Ttopright_extensions(pt).pts(iL+2) Ttopright_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
        end
    end
end


foo = polybuffer(lattice, 'lines', 1, 'JointType', 'square');
figure(23)
plot(foo)
hold on
plot(sideRightPanel(1))
hold off
axis image


% save the STL file
[f, v] = extrude(foo, 0.5);        
% saveModel('panelSideRight.stl', f, v);

%%
% create Top Panel

v = o(1).v;
e = o(1).e;
eSeamExt = [o(2).eSeamExtRot; o(3).eSeamExtRot];

v(:,1) = v(:,1) + min(topOutline(:,1));
v(:,2) = v(:,2) + min(topOutline(:,2));
eSeamExt(:,[1 3]) = eSeamExt(:,[1 3]) + min(topOutline(:,1));
eSeamExt(:,[2 4]) = eSeamExt(:,[2 4]) + min(topOutline(:,2));

lattice = [];
for iE=1:size(e,1)
    lattice = [lattice; v(e(iE,1),:); v(e(iE,2),:); [NaN NaN] ];
end
for iS=1:size(eSeamExt,1)
    lattice = [lattice; eSeamExt(iS,1:2); eSeamExt(iS,3:4); [NaN NaN] ];
end

stateLattice.topPanel = lattice;

for pt = 1:size(Tleft_extensions,2)
    if ~isempty(Tleft_extensions(pt).pts)
        pt
        for iL = 1:4:length(Tleft_extensions(pt).pts)
            if Tleft_extensions(pt).pts
            lattice = [lattice; [Tleft_extensions(pt).pts(iL) Tleft_extensions(pt).pts(iL+1)];...
                [Tleft_extensions(pt).pts(iL+2) Tleft_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
        end
    end
end

for pt = 1:size(Tright_extensions,2)
    if ~isempty(Tright_extensions(pt).pts)
        pt
        for iL = 1:4:length(Tright_extensions(pt).pts)
            if Tright_extensions(pt).pts
            lattice = [lattice; [Tright_extensions(pt).pts(iL) Tright_extensions(pt).pts(iL+1)];...
                [Tright_extensions(pt).pts(iL+2) Tright_extensions(pt).pts(iL+3)]; [NaN NaN]];
            end
        end
    end
end

foo = polybuffer(lattice, 'lines', 1, 'JointType', 'square');
% stateLattice.topPanel = foo;
figure(25)
plot(foo)
hold on
plot(topPanel(1))
hold off
% foo = polybuffer(topleftExtension_lattice, 'lines', 1, 'JointType', 'square');
% hold on
% plot(foo)
% hold off
% foo = polybuffer(toprightExtension_lattice, 'lines', 1, 'JointType', 'square');
% hold on
% plot(foo)
% hold off
% foo = polybuffer(TleftExtension_lattice, 'lines', 1, 'JointType', 'square');
% hold on
% plot(foo)
% hold off
axis image
max_size = [270 270]; % mm
buffer = 2; % mm in overlap along boundaries
fooCut = cut(foo, max_size, buffer);

figure(31)
subplot(1,2,1)
plot(fooCut(1))
subplot(1,2,2)
plot(fooCut(2))

%%

% save the STL file
[f, v] = extrude(fooCut(1), 0.5);        
% saveModel('panelTop_1.stl', f, v);

[f, v] = extrude(fooCut(2), 0.5);        
% saveModel('panelTop_2.stl', f, v);
%%
latticeExtensions.leftExtensions = left_extensif o(iO).eSeamX(pt)==1
ions;
latticeExtensions.rightExtensions = right_extensions;
latticeExtensions.topLeftExtensions = topleft_extensions;
latticeExtensions.topRightExtensions = topright_extensions;
latticeExtensions.TleftExtensions = Tleft_extensions;
latticeExtensions.TrightExtensions = Tright_extensions;
latticeExtensions.TtopLeftExtensions = Ttopleft_extensions;
latticeExtensions.TtopRightExtensions = Ttopright_extensions;
latticeExtensions.lefteSeamXpts = o(2).eSeamXpts;
latticeExtensions.righteSeamXpts = o(3).eSeamXpts;
latticeExtensions.toplefteSeamXpts = o(2).eSeamXptsRot;
latticeExtensions.toprighteSeamXpts = o(3).eSeamXptsRot;

%%
function pt = get_point_on_outline(outline_points, dist)
    
    dist_along_outline = 0;
    for u = 1:size(outline_points,1)-1
        dist_along_outline = dist_along_outline+sqrt(sum((outline_points(u,:)-outline_points(u+1,:)).^2));
        if dist < dist_along_outline
            leftover_dist = dist_along_outline-dist;
            unit_vector = (outline_points(u,:)-outline_points(u+1,:))/norm(outline_points(u,:)-outline_points(u+1,:));
            pt = outline_points(u+1,:) +(unit_vector*leftover_dist);
            break;
        end
        pt = outline_points(u+1,:);
    end
  
end
