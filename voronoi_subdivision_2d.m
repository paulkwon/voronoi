function world = voronoi_subdivision_2d
close all
addpath(fullfile(getenv('ywtoolbox'),'ocTree'));

n = 7;
o = rand(2,n); r = 0.02+0.05*rand(1,n);
show_objects(o,r);
% show_cont_grad(o,r);

% keyboard
world = ywOcTree2d('maxDepth',6);
world.Objects = {o,r};
world.divide(1);
world.plot();

% visualize octree construction
cc = get(groot, 'DefaultAxesColorOrder');
binLeaf = find(world.getIdxIsLeaf());
for i = 1:numel(binLeaf)
    binNo = binLeaf(i);
    binMinMax = world.BinBoundaries(binNo,:);
    pts = binMinMax([1 2; 1 4; 3 4; 3 2])';
    objID = world.queryNN(pts);
    if all(objID==objID(1)), continue; end
    scatter(pts(1,:), pts(2,:), 5, cc(objID,:), 'fill');
    patch(pts(1,:), pts(2,:), reshape(cc(objID,:),4,1,3),'facealpha',0.3,'edgecolor','interp');
end

%% visualize topology
binLeaf = find(world.getIdxIsLeaf());
E=zeros(2,0); V=zeros(2,0);
for i = 1:numel(binLeaf)
    binNo = binLeaf(i);
    binMinMax = world.BinBoundaries(binNo,:);
    pts = binMinMax([1 2; 1 4; 3 4; 3 2])';
    objID = world.queryNN(pts);
    if all(objID==objID(1)), continue; end
    
    if numel(unique(objID))==1, continue; end;
    
    faces = [1 2; 2 3; 3 4; 4 1];
    cut = objID(faces(:,1))~=objID(faces(:,2));
    voronoi_pts = (pts(:,faces(cut,1)) + pts(:,faces(cut,2)))/2;
    E = [E voronoi_pts];
    if numel(unique(objID))>=3,          
        V = [V mean(pts,2)];
    end
end
scatter(E(1,:), E(2,:), 5, 'k', 'fill');
scatter(V(1,:), V(2,:), 8, '^k', 'fill');

function show_objects(o,r)

cc = get(groot, 'DefaultAxesColorOrder');
n = size(o,2);
figure; scatter(o(1,:),o(2,:),20, cc(1:n,:), 'fill'); axis([0 1 0 1]); axis equal
t=0:pi/100:2*pi; 
hold on; plot(repmat(o(1,:),numel(t),1)+repmat(r,numel(t),1).*repmat(cos(t)',1,n), repmat(o(2,:),numel(t),1)+repmat(r,numel(t),1).*repmat(sin(t)',1,n))

function [x,y,px,py]=show_cont_grad(o,r)
[x,y] = meshgrid(0:0.025:1, 0:0.025:1);
z=zeros(size(x));
for i=1:numel(x)
    p = [x(i);y(i)];
    z(i) = objective(p,o,r);
end
hold on;
contour(x,y,z);
[px,py] = gradient(z,0.025,0.025);    % numerical gradient
quiver(x,y,px,py);

function val = objective(p,o,r)
n = size(o,2);
delta = repmat(p,1,n)-o;
%     hold on; quiver(p(1)*ones(1,n), p(2)*ones(1,n), delta(1,:),delta(2,:),'r');
dist = sqrt(sum(delta.^2)) - r;
[~,ind] = sort(dist);           % indices nearest
nn = 3;
dist_nn = dist(ind(1:nn));        % distances to nearest 3
val = sum(dist_nn.^2)-1/nn*sum(dist_nn)^2;