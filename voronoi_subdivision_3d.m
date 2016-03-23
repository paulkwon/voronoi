function voronoi_subdivision_3d
close all
addpath(fullfile(getenv('ywtoolbox'),'ocTree'));

n = 3;
o = rand(3,n); r = 0.05+0.1*rand(1,n);
show_objects(o,r);
% show_cont_grad(o,r);

world = ywOcTree('maxDepth',7);
world.Objects = {o,r};
world.divide(1);
% world.plot();
view(45,45)

cc = get(groot, 'DefaultAxesColorOrder');
binLeaf = find(world.getIdxIsLeaf());
% for i = 1:numel(binLeaf)
%     binNo = binLeaf(i);
%     binMinMax = world.BinBoundaries(binNo,:);
%     vertices = binMinMax([1 2 3; 1 2 6; 1 5 6; 1 5 3; 4 2 3; 4 2 6; 4 5 6; 4 5 3])';
%     objID = world.queryNN(vertices);
%     if all(objID==objID(1)), continue; end
%     
%     %     patch([0 0 1 1],[0 2 3 0],[1 1 1 1], reshape([0 0 1;0 0 1;0 1 0;1 0 0],4,1,3),'facealpha',0.3,'edgecolor','interp')
%     %     patch(pts(1,:), pts(2,:), pts(3,:), reshape(cc(objID,:),4,1,3),'facealpha',0.3,'edgecolor','interp');
%     faces = [1 2 3 4;2 3 7 6; 1 4 8 5; 3 7 8 4; 5 6 7 8; 1 2 6 5];
%     for k=1:6
%         m = faces(k,:);
%         face_vertices = vertices(:,m);
% %         patch(face_vertices(1,:), face_vertices(2,:), face_vertices(3,:), reshape(cc(objID(m),:),4,1,3),'facealpha',0.3,'edgecolor','interp');
%         %         face_edge = [1 2; 2 3; 3 4; 4 1];
%         %         cut = objID(m) ~= objID([m(2:end) m(1)]);
%         %         cut_vtx = (face_vertices(:,face_edge(cut,1))+face_vertices(:,face_edge(cut,2)))/2;
%         %         plot3(cut_vtx(1,:),cut_vtx(2,:),cut_vtx(3,:),'k');
%     end
%     
%     centeroid = (binMinMax(1:3)+binMinMax(4:6))'/2;
%     edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
%     cut = objID(edges(:,1))~=objID(edges(:,2));
%     voronoi_pts = (vertices(:,edges(cut,1)) + vertices(:,edges(cut,2)))/2;
%     for b=1:size(voronoi_pts,2)
%         draw_barrier(centeroid, voronoi_pts(:,b));
%     end 
% end

%%
B = [];
for i = 1:numel(binLeaf)
    binNo = binLeaf(i);
    binMinMax = world.BinBoundaries(binNo,:);
    vertices = binMinMax([1 2 3; 1 2 6; 1 5 6; 1 5 3; 4 2 3; 4 2 6; 4 5 6; 4 5 3])';
    objID = world.queryNN(vertices);
    if all(objID==objID(1)), continue; end
    centeroid = (binMinMax(1:3)+binMinMax(4:6))'/2;
    edges = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
    cut = objID(edges(:,1))~=objID(edges(:,2));
    voronoi_pts = (vertices(:,edges(cut,1)) + vertices(:,edges(cut,2)))/2;
    b = [repmat(centeroid, 1, size(voronoi_pts,2)); voronoi_pts];
    B=[B b];
end

%% file write
fid = fopen('voronoi.dat','wb');
fwrite(fid, size(world.Objects{1},2), 'integer*4');
fwrite(fid, [world.Objects{1}; world.Objects{2}],'float');

fwrite(fid, size(B,2), 'integer*4');
fwrite(fid, B, 'float');
fclose(fid);
fprintf('%d objects.\n %d boxes.\n', size(world.Objects{1},2), size(B,2));

function draw_barrier(p1, p2)
% mn = min(p1,p2);
% mx = max(p1,p2);
same_dimension = p1==p2;
assert(~isempty(find(same_dimension, 1)));
diff_dimension = ~same_dimension;

mn = min(p1(diff_dimension),p2(diff_dimension));
mx = max(p1(diff_dimension),p2(diff_dimension));
box = [mn(1) mn(1) mx(1) mx(1); ...
    mn(2) mx(2) mx(2) mn(2)];
barrier_pts = zeros(3,4);
barrier_pts(diff_dimension, :) = box;
barrier_pts(same_dimension, :) = p1(same_dimension);
plot3([barrier_pts(1,:); barrier_pts(1,[2 3 4 1])],[barrier_pts(2,:); barrier_pts(2,[2 3 4 1])],[barrier_pts(3,:); barrier_pts(3,[2 3 4 1])],'k')


%% visualize topology
% binLeaf = find(world.getIdxIsLeaf());
% E=zeros(2,0); V=zeros(2,0);
% for i = 1:numel(binLeaf)
%     binNo = binLeaf(i);
%     binMinMax = world.BinBoundaries(binNo,:);
%     pts = binMinMax([1 2; 1 4; 3 4; 3 2])';
%     objID = world.queryNN(pts);
%     if all(objID==objID(1)), continue; end
%
%     if numel(unique(objID))==1, continue; end;
%
%     faces = [1 2; 2 3; 3 4; 4 1];
%     cut = objID(faces(:,1))~=objID(faces(:,2));
%     voronoi_pts = (pts(:,faces(cut,1)) + pts(:,faces(cut,2)))/2;
%     E = [E voronoi_pts];
%     if numel(unique(objID))>=3,
%         V = [V mean(pts,2)];
%     end
% end
% scatter(E(1,:), E(2,:), 5, 'k', 'fill');
% scatter(V(1,:), V(2,:), 8, '^k', 'fill');


function show_objects(o,r)

cc = get(groot, 'DefaultAxesColorOrder');
n = size(o,2);
% figure; scatter3(o(1,:), o(2,:), o(3,:), 20, cc(1:n,:), 'fill');
figure;
[x,y,z] = sphere;
for i=1:n
    hold on; surf(r(i)*x+o(1,i),r(i)*y+o(2,i),r(i)*z+o(3,i),  'FaceColor',cc(i,:))
end
axis([0 1 0 1 0 1]);

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