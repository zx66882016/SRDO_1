function info = uav_eval_path(path, obstacles, cfg)

if nargin < 3, cfg = uav_config(); end
path = max(min(path,cfg.space),0);

seg = diff(path);
total_len = sum(sqrt(sum(seg.^2,2)));

centers = obstacles(:,1:3); radii = obstacles(:,4);
dx = path(:,1)-centers(:,1)';
dy = path(:,2)-centers(:,2)';
dz = path(:,3)-centers(:,3)';
dist = sqrt(dx.^2+dy.^2+dz.^2);
hard = dist < radii';
n_hard = sum(hard(:));
penalty = n_hard*cfg.pen_hard;
soft = (dist>=radii') & (dist<radii'+cfg.safe_dist);
if any(soft(:))
    penalty = penalty + cfg.pen_soft*sum((radii'+cfg.safe_dist-dist).*soft,'all');
end
min_clearance = min(min(dist-radii',[],2));

v1 = path(2:end-1,:)-path(1:end-2,:);
v2 = path(3:end,:)-path(2:end-1,:);
n1 = sqrt(sum(v1.^2,2)); n2 = sqrt(sum(v2.^2,2));
cos_a = sum(v1.*v2,2)./(n1.*n2+eps);
cos_a = max(-1,min(1,cos_a));
curv = sum(acos(cos_a).^2);

axis_u = [1 1 1]/sqrt(3);
prog = path*axis_u(:);
backtrack = sum(max(0, prog(1:end-1)-prog(2:end)));

cost = cfg.w_len*total_len + cfg.w_obs*penalty + cfg.w_smooth*curv + cfg.w_mono*backtrack;
info = struct('path',path,'len',total_len,'curv',curv,'backtrack',backtrack, ...
    'n_hard',n_hard,'feasible',n_hard==0,'min_clearance',min_clearance,'cost',cost);
end
