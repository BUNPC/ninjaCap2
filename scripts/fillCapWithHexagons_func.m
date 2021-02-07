function [v, e, vOut, eOut] = fillCapWithHexagons_func( v, hEdge, Imask )

%%
% fill with hexagons
% initialize the first vertex and edges

voe = 0; % 0-horizontal edge points right. 1-horizontal edge points left

% three edges from each vertex
dv = [hEdge 0; hEdge*cos(2*pi/3) hEdge*sin(2*pi/3); hEdge*cos(2*pi/3) -hEdge*sin(2*pi/3) ];

oe = mod(voe+1,2); % flip from odd to even to change direction of horizonal edge
[ny,nx] = size(Imask); % dimensions of the mask

% this is the code to add the 3 new vertices from the first vertex and
% remove the ones that already exist. The redundancy check is not needed
% for this first time, but it tests the code for the loop later. This also
% makes sure the added vertices are within the mask.
vtmp = ones(3,1)*v + dv;
vCen = []; vCen(1,1:2) = v - dv(1,:);
vtmp = setdiff(vtmp,v,'rows','stable'); % find vertices that already exist and remove them
lst = [];
for ii=1:size(vtmp,1) % make sure the vertices are within the mask
    if vtmp(ii,1)<1 || vtmp(ii,1)>nx || vtmp(ii,2)<1 || vtmp(ii,2)>ny
        lst(ii) = 0;
    elseif  Imask(round(vtmp(ii,2)),round(vtmp(ii,1)))==1
        lst(ii) = 1;
    else
        lst(ii) = 0;
    end
end
vtmp = vtmp(find(lst==1),:);

nv = size(v,1);
nvtmp = size(vtmp,1);

vlst = nv + [1:nvtmp]; % this is the list of new vertices that we iterate over to add the next round of vertices
v(end+[1:nvtmp],:) = vtmp; % append new vertices to the vertex list
voe(end+[1:nvtmp]) = oe; % store if vertex is even or odd, that is right or left horizontal edge
e = [1 2; 1 3; 1 4]; % initial edge list

%%
% iterate the hexagon filling process

eOut = [];
vOut = [];
while ~isempty(vlst) % loop until no more new vertices to expand from

    if voe(vlst(end))==0 % added new vertices either with right or left horizontal edge
        vtmp = ones(3,1)*v(vlst(end),:) + dv;
        vCen(end+1,1:2) = v(vlst(end),:) - dv(1,:);
    else
        vtmp = ones(3,1)*v(vlst(end),:) - dv;
        vCen(end+1,1:2) = v(vlst(end),:) + dv(1,:);
    end
    oe = mod(voe(vlst(end))+1,2); % swap the direction of the next horizontal edge

    % add edges to the new vertices and draw them on the mask
    % but only add it if it doesn't already exist in the edge list
    nv = size(v,1);
    nVexist = 0;
    for ii = 1:size(vtmp,1)
        
        % add edges
        v1 = vlst(end);        
        [v3,jj,kk] = intersect(vtmp(ii,:),v,'rows');
        flag = 1;
        if isempty(kk)
            
            flag = 0;
            if vtmp(ii,1)<1 || vtmp(ii,1)>nx || vtmp(ii,2)<1 || vtmp(ii,2)>ny
                flag = 0;
            elseif  Imask(round(vtmp(ii,2)),round(vtmp(ii,1)))==1
                flag = 1;
            else
                flag = 0;
            end
        
            if flag==1
                v2 = nv+ii-nVexist;
            else
                v2 = 0;
                nVexist = nVexist + 1;
            end
        else
            v2 = kk;
            nVexist = nVexist + 1;
        end
        if isempty( find( (e(:,1)==v1 & e(:,2)==v2) | (e(:,1)==v2 & e(:,2)==v1) ) ) && v2>0
            e(end+1,:) = [v1 v2];
        end
        
%         figure(1)
%         hold on
        p1 = v(vlst(end),:);
        p2 = vtmp(ii,:);
        if flag==1
%             hl=plot( [p1(1) p2(1)], [p1(2) p2(2)], 'k-');
        else
%             hl=plot( [p1(1) p2(1)], [p1(2) p2(2)], 'r-');
            
            if ~isempty(vOut)
                ii1 = find( sum((ones(size(vOut,1),1)*p1 - vOut).^2,2).^0.5<1e-4, 1 );
                ii2 = find( sum((ones(size(vOut,1),1)*p2 - vOut).^2,2).^0.5<1e-4, 1 );
                if isempty(ii1), ii1=0; end
                if isempty(ii2), ii2=0; end
            else
               ii1 = 0;
               ii2 = 0;
            end
            if ii1==0
                vOut(end+1,1:2) = p1; 
                iE1 = size(vOut,1);
            else
                iE1 = ii1;
            end
% I commented this because I want to double count the vertex outside the polygon 
% because when the edges are brought back to the outline I will need to
% vertices
%            if ii2==0 | isempty(ii2) 
                vOut(end+1,1:2) = p2; 
                iE2 = size(vOut,1);
 %           else
 %               iE2 = ii2;
 %           end
            if isempty(eOut)
                eOut(end+1,1:2) = [iE1 iE2];
            elseif ~ismember([iE1 iE2],eOut,'rows')
                eOut(end+1,1:2) = [iE1 iE2];
            end
        end
%         set(hl,'linewidth',2)
%         hold off
    end
        
    % remove redundant vertices and those that fall outside the mask
    % need to remove edge for those that fall outside the mask
    vtmp = setdiff(vtmp,v,'rows','stable');
    lst = [];
    for ii=1:size(vtmp,1)
        if vtmp(ii,1)<1 || vtmp(ii,1)>nx || vtmp(ii,2)<1 || vtmp(ii,2)>ny
            lst(ii) = 0;
        elseif  Imask(round(vtmp(ii,2)),round(vtmp(ii,1)))==1
            lst(ii) = 1;
        else
            lst(ii) = 0;
        end
        if lst(ii)==0 % remove the edge
            [ir,ic] = find(e==(nv+ii));
            e(ir,:) = [];
        end
    end
    vtmp = vtmp(find(lst==1),:);
    
    nvtmp = size(vtmp,1);
    
    vlst(end) = []; % remove last vertex from the new vertex list
    vlst(end+[1:nvtmp]) = nv + [1:nvtmp]; % add the new vertices to the new vertex list
    v(end+[1:nvtmp],:) = vtmp; % add new vertices to the vertex array
    voe(end+[1:nvtmp]) = oe; % store if vertex is even or odd, that is right or left horizontal edge
    
    % draw the vertices on the mask
%     if nvtmp>0
%         figure(1)
%         hold on
%         plot(v(:,1),v(:,2),'r*')
%         plot(vCen(:,1),vCen(:,2),'b*')
%         hold off
%         drawnow
%     end

end  