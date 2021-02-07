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