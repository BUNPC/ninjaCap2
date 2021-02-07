function svgExport(fname, pts, pathid)
% Writes a very basic svg file, drawing lines point by point
% for import and editing in SVG tools such as inkscape
%
% INPUT ARGUMENTS:
% fname     - target file name
% pts       - outline, vector of points [x, y]
% pathid    - id of the path (string)

%  if no points are given use original default outline
if isempty(pts)
    pts = [1.20086, 8.090358;...
    2.0043, 6.786918;...
    5.368705, 6.083908;...
    5.368705, 4.778318;...
    2.0043, 4.376598;...
    1, 3.053072;...
    1.954085, 2.269718;...
    5.167845, 2.553072;...
    3.26088016, 1];
end

% open file
fid=fopen(fname, 'w');

sz=max(pts);
fprintf(fid, '<svg height="%1.4f" width="%1.4f">\n', sz(1), sz(2));


%% draw circles for each datapoint as a reference

%circle radius:
rad=0.1;
lwidth=0.01;
for ii=1:size(pts,1)
    fprintf(fid, '<circle \n');
    
    fprintf(fid, 'cx="%1.4f" \n', pts(ii,1));
    fprintf(fid, 'cy="%1.4f" \n', pts(ii,2));
    fprintf(fid, 'r="%1.4f" \n', rad);
    fprintf(fid, 'stroke="black" \n');
    fprintf(fid, 'stroke-width="%1.4f" \n', lwidth);
    fprintf(fid, 'fill="black" />\n');
end


%% draw path from cap outline points

fprintf(fid, '\n<path\n');
fprintf(fid, 'd="');
%start point (first coordinate)
fprintf(fid, "M %1.4f %1.4f ", pts(1,1), pts(1,2));
for ii=2:size(pts,1)
    fprintf(fid, 'L%1.4f %1.4f ', pts(ii,1), pts(ii,2));
end
fprintf(fid, ' "\n');
fprintf(fid, 'style="fill:none;stroke:#000000;stroke-opacity:1;stroke-width:%1.2f" \n', lwidth);
fprintf(fid, ['id="' pathid '" /> \n']);

fprintf(fid,'</svg>\n');

%% close file
fclose(fid);


end

