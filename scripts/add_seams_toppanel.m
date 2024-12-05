function [topPanel_left, topPanel_right,topLoopPos] = add_seams_toppanel(topPanel_left, topPanel_right,topPanel_midpoint_x)


ymin = false;
ymax = false;
seam_width = 1.5;
for u = 1:size(topPanel_left,2)
    if ~ isempty(topPanel_left(1,u).Vertices)
        if islogical(ymin)
            ymin = min(topPanel_left(1,u).Vertices(:,2));
            ymax = max(topPanel_left(1,u).Vertices(:,2));
        else
            ymin = min(ymin,min(topPanel_left(1,u).Vertices(:,2)));
            ymax = max(ymax,max(topPanel_left(1,u).Vertices(:,2)));
        end
    end
end
seam_top = polyshape([topPanel_midpoint_x-seam_width ymin; topPanel_midpoint_x+seam_width ymin; topPanel_midpoint_x+seam_width ymax+5; topPanel_midpoint_x-seam_width ymax+5]);
LoopPos = [];
for u = linspace(ymin+5, ymax-20,10)
    LoopPos = [LoopPos ; topPanel_midpoint_x u];
end
topLoopPos{1} = LoopPos;
topLoopPos{3} = LoopPos;


ymin = false;
ymax = false;
for u = 1:size(topPanel_left,2)
    if ~ isempty(topPanel_left(2,u).Vertices)
        if islogical(ymin)
            ymin = min(topPanel_left(2,u).Vertices(:,2));
            ymax = max(topPanel_left(2,u).Vertices(:,2));
        else
            ymin = min(ymin,min(topPanel_left(2,u).Vertices(:,2)));
            ymax = max(ymax,max(topPanel_left(2,u).Vertices(:,2)));
        end
    end
end
seam_bottom = polyshape([topPanel_midpoint_x-seam_width ymin-5; topPanel_midpoint_x+seam_width ymin-5; topPanel_midpoint_x+seam_width ymax; topPanel_midpoint_x-seam_width ymax]);
LoopPos = [];
for u = linspace(ymin+20, ymax-5,10)
    LoopPos = [LoopPos ; topPanel_midpoint_x u];
end
topLoopPos{2} = LoopPos;
topLoopPos{4} = LoopPos;

topPanel_left(:,2) = [seam_top; seam_bottom];
topPanel_right(:,2) = [seam_top; seam_bottom];

