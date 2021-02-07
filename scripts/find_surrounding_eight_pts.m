vHead = headsurf.mesh.vertices;
fHead = headsurf.mesh.faces;
%%
load('vHead.mat')
load('fHead.mat')
%%
nV = size(vHead,1);

h = 20;

% position of the grommet
p = [  220.6706  161.7568   49.9768];

% find nearest vertex on the head
lst = find( abs(vHead(:,1)-p(1)*ones(nV,1))<h & abs(vHead(:,2)-p(2)*ones(nV,1))<h & abs(vHead(:,3)-p(3)*ones(nV,1))<h );
rho = sum( (ones(length(lst),1)*p - vHead(lst,:)).^2, 2 ).^0.5;
[foo,idx] = min(rho);

iV = lst(idx);

% get surface normal
[ir,ic]=find(fHead==iV);

iF = ir(1);
if ic(1)==1
    p1 = vHead(fHead(iF,2),:);
    p2 = vHead(fHead(iF,3),:);
elseif ic(1)==2
    p1 = vHead(fHead(iF,1),:);
    p2 = vHead(fHead(iF,3),:);
elseif ic(1)==3
    p1 = vHead(fHead(iF,1),:);
    p2 = vHead(fHead(iF,2),:);
end

n = cross( p1-p, p2-p );
n = n / norm(n);

% get tangential axes
a1 = p1-p; 
a1 = a1 / norm(a1);
a2 = cross( a1, n );
a2 = a2 / norm(a2);

%%
% find points along positive a1 from p
lstSub = find( abs(vHead(:,1)-p(1)*ones(nV,1))<h & abs(vHead(:,2)-p(2)*ones(nV,1))<h & abs(vHead(:,3)-p(3)*ones(nV,1))<h );
lst = find( vHead(lstSub,:)*a1' > 0 );
rho = sum( (ones(length(lst),1)*p - vHead(lstSub(lst),:)).^2, 2 ).^0.5;

