function [outline_pts, outline_labels] = getEquiDistPtsAlongOutline(outlinePts,dist)

left_over_dist = 0;
current_pt = 1;
total_points = size(outlinePts,1);
outline_pts = outlinePts(1,:);
outline_labels = {'OL1'};
while current_pt < total_points
    next_pt_dist = sqrt(sum((outlinePts(current_pt,:)-outlinePts(current_pt+1,:)).^2,2));
    left_over_dist = left_over_dist+next_pt_dist;
    if left_over_dist >= dist
        npoint = 0;
        dist_along_line = dist-(left_over_dist-next_pt_dist);
        while left_over_dist >= dist
            normalised_vector = (outlinePts(current_pt+1,:)-outlinePts(current_pt,:))/next_pt_dist;
            outline_pts = [outline_pts; outlinePts(current_pt,:)+(dist_along_line+npoint*dist)*normalised_vector];
            outline_labels{end+1} = strcat('OL',num2str(size(outline_pts,1)));
            left_over_dist = left_over_dist-dist;
            npoint = npoint+1;
        end
    end
    current_pt = current_pt+1;
end