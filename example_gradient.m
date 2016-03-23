function example_gradient
close all

%%
n = 2;
o = rand(2,n); r = 0.05+0.1*rand(1,n);

figure; scatter(o(1,:),o(2,:),'fill'); axis([0 1 0 1]); axis equal
t=0:pi/100:2*pi; 
hold on; plot(repmat(o(1,:),numel(t),1)+repmat(r,numel(t),1).*repmat(cos(t)',1,n), repmat(o(2,:),numel(t),1)+repmat(r,numel(t),1).*repmat(sin(t)',1,n))

[x,y] = meshgrid(0:0.025:1, 0:0.025:1);
z=zeros(size(x));
u=zeros(size(x));
v=zeros(size(x));
for i=1:numel(x)
    p = [x(i);y(i)];
    delta = repmat(p,1,n)-o;
    dist = max(sqrt(sum(delta.^2)) - r,0);

    % z(i) = sum((dist-meandist).^2); % original expression
%     z(i) = sum(dist.^2)- 1/n*sum(dist)^2;
    
    grad = sum(delta,2) -1/3*sum(dist)*sum(repmat(dist.^-1,2,1).*delta,2);
    u(i) = grad(1);
    v(i) = grad(2);
end
hold on;
contour(x,y,z);
quiver(x,y,u,v);

figure
[px,py] = gradient(z,0.025,0.025);    % numerical gradient
quiver(x,y,px,py);

function val = objective(p,o,r)

n = size(o,2);
delta = repmat(p,1,n)-o;
%     hold on; quiver(p(1)*ones(1,n), p(2)*ones(1,n), delta(1,:),delta(2,:),'r');
dist = sqrt(sum(delta.^2)) - r;
[~,ind] = sort(dist);           % indices of nearest 3
delta_3 = delta(:,ind(1:3));    % vectors to nearest 3
dist_3 = dist(ind(1:3));        % distances to nearest 3
val = sum(dist_3.^2)-1/3*sum(dist_3)^2;