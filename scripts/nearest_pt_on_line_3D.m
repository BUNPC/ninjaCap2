function nearest_pt = nearest_pt_on_line_3D(line_pt1,line_pt2,pt)

t =  ((pt(1)-line_pt1(1))*(line_pt2(1)-line_pt1(1))+(pt(2)-line_pt1(2))*(line_pt2(2)-line_pt1(2))+(pt(3)-line_pt1(3))*(line_pt2(3)-line_pt1(3)))/norm(line_pt2-line_pt1)^2;
nearest_pt = [t*(line_pt2(1)-line_pt1(1))+line_pt1(1) t*(line_pt2(2)-line_pt1(2))+line_pt1(2) t*(line_pt2(3)-line_pt1(3))+line_pt1(3)];