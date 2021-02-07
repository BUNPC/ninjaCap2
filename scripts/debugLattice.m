function debugLattice(poly, lattice)
%DEBUGLATTICE Summary of this function goes here
%   Detailed explanation goes here

hold on;
plot(poly);
plot(lattice(:, 1), lattice(:, 2), 'k');
hold off;

end

