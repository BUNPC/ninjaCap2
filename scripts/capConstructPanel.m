function panel = capConstructPanel(poly, outline_out, width_outline_out, outline_stitch, width_outline_stitch, lattice, width_lattice, outline_extensions, extensions, width_extensions, stitch_nubs, width_nub_stitch)
%CAPCONSTRUCTPANEL Summary of this function goes here
%   Detailed explanation goes here

panel = polyshape();

% add outline outside
panel(1, 1) = polybuffer(outline_out, 'lines', width_outline_out);

% add outline stitching
if width_outline_stitch > 0
    panel(1, 2) = polybuffer(outline_stitch, 'lines', width_outline_stitch);
else
    panel(1, 2) = polyshape();
end

% add online lattice
panel(1, 3) = polybuffer(lattice, 'lines', width_lattice / 2, 'JointType', 'square');

% add stitching nubs
if width_nub_stitch > 0
    panel(1, 4) = polybuffer(stitch_nubs, 'points', width_nub_stitch / 2);
else
    panel(1, 4) = polyshape();
end

% add placeholder for grommets
panel(1, 5) = polyshape();

% clip
panel = intersect(panel, poly); % clip

% add extensions
if width_extensions > 0
    if ~isempty(outline_extensions)
        panel(1, 1) = union(panel(1, 1), polybuffer(outline_extensions, 'lines', width_extensions / 2, 'JointType', 'round'));
    end
    if ~isempty(extensions)
        panel(1, 3) = union(panel(1, 3), polybuffer(extensions, 'lines', width_extensions / 2, 'JointType', 'round'));
    end
end

end
