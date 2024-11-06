function hHex = getGeodesicdist(vertices, faces, vHex, eHex)

global geodesic_library;                
geodesic_library = 'geodesic_debug'; 

mesh = geodesic_new_mesh(vertices,faces);         %initilize new mesh
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm
hHex = zeros(size(eHex,1),1);
figure
for u = 1:size(eHex,1)
    [u size(eHex,1)]
    pos1 = vHex(eHex(u,1),:);
    pos2 = vHex(eHex(u,2),:);

    vert_dist = sqrt(sum((vertices-pos1).^2,2));
    [min_dist, vertex_id] = min(vert_dist);

    source_points = {geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:))};

    geodesic_propagate(algorithm, source_points,[], 40);   %propagation stage of the algorithm (the most time-consuming)
    
     vert_dist = sqrt(sum((vertices-pos2).^2,2));
    [min_dist, vertex_id] = min(vert_dist);

    destination = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));

%%
%     vert_dist = sqrt(sum((vertices-pos2).^2,2));
%     [min_dist, vertex_id] = min(vert_dist); 
%     destination = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));
%     
%     vert_dist = sqrt(sum((vertices-pos1).^2,2));
%     [min_dist, vertex_id] = mink(vert_dist,3);
%     
%     face_id = find((faces(:,1)==vertex_id(1) & faces(:,2)==vertex_id(2) & faces(:,3)==vertex_id(3)) | ...
%         (faces(:,1)==vertex_id(1) & faces(:,2)==vertex_id(3) & faces(:,3)==vertex_id(2)) |...
%         (faces(:,1)==vertex_id(2) & faces(:,2)==vertex_id(1) & faces(:,3)==vertex_id(3)) |...
%         (faces(:,1)==vertex_id(2) & faces(:,2)==vertex_id(3) & faces(:,3)==vertex_id(1)) |...
%         (faces(:,1)==vertex_id(3) & faces(:,2)==vertex_id(1) & faces(:,3)==vertex_id(2)) |...
%         (faces(:,1)==vertex_id(3) & faces(:,2)==vertex_id(2) & faces(:,3)==vertex_id(1)));
%     if isempty(face_id)
%             vert_dist = sqrt(sum((vertices-pos1).^2,2));
%             [min_dist, vertex_id] = min(vert_dist);
%             source_points = {geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:))};
%     else
%     source_points = {geodesic_create_surface_point('face',face_id,pos1)};
%     end
    %%
%     geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)
%     
%     vert_dist = sqrt(sum((vertices-pos2).^2,2));
%     [min_dist, vertex_id] = mink(vert_dist,3);
%     
%     face_id = find((faces(:,1)==vertex_id(1) & faces(:,2)==vertex_id(2) & faces(:,3)==vertex_id(3)) | ...
%         (faces(:,1)==vertex_id(1) & faces(:,2)==vertex_id(3) & faces(:,3)==vertex_id(2)) |...
%         (faces(:,1)==vertex_id(2) & faces(:,2)==vertex_id(1) & faces(:,3)==vertex_id(3)) |...
%         (faces(:,1)==vertex_id(2) & faces(:,2)==vertex_id(3) & faces(:,3)==vertex_id(1)) |...
%         (faces(:,1)==vertex_id(3) & faces(:,2)==vertex_id(1) & faces(:,3)==vertex_id(2)) |...
%         (faces(:,1)==vertex_id(3) & faces(:,2)==vertex_id(2) & faces(:,3)==vertex_id(1)));
%     
%       if isempty(face_id)
%         vert_dist = sqrt(sum((vertices-pos2).^2,2));
%         [min_dist, vertex_id] = min(vert_dist);
%         destination = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));
%      else
%         destination = geodesic_create_surface_point('face',face_id,pos2);
%      end
  %%  
    
    path = geodesic_trace_back(algorithm, destination);     %find a shortest path from source to destination

    [source_id, distances] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1
    
    dist = 0;
    start_point = [path{1}.x path{1}.y path{1}.z];
    for v = 2:length(path)
        current_point = [path{v}.x path{v}.y path{v}.z];
        dist = dist+sqrt(sum((start_point-current_point).^2));
        start_point = current_point;
    end
    hHex(u) = dist;
    
    if 1
        if  u == 1
            colormap('default');
            trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),distances, 'FaceColor', 'interp', 'EdgeColor', 'k');       %plot the mesh
            daspect([1 1 1]);
        end

        hold on;
        plot3(source_points{1}.x, source_points{1}.y, source_points{1}.z, 'or', 'MarkerSize',3);    %plot sources

        plot3(destination.x, destination.y, destination.z, 'ok', 'MarkerSize',3);       %plot destination 
        [x,y,z] = extract_coordinates_from_path(path);                                  %prepare path data for plotting
        h = plot3(x*1.001,y*1.001,z*1.001,'r-','LineWidth',2);    %plot path
        legend(h,'geodesic curve');
    end
end