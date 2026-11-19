function [cost, info] = uav_cost3d(ctrl_pts, n_ctrl, start, goal, obstacles, cfg)


if nargin < 6, cfg = uav_config(); end

P = reshape(ctrl_pts(:), n_ctrl, 3);
P_full = [start; P; goal];              
n = size(P_full,1);
t  = 0:n-1;
tt = linspace(0,n-1,cfg.n_sample);
path = interp1(t, P_full, tt, 'spline');  
path = max(min(path, cfg.space), 0);      

seg = diff(path);
seg_len = sqrt(sum(seg.^2,2));
total_len = sum(seg_len);

centers = obstacles(:,1:3);
radii   = obstacles(:,4);
dx = path(:,1)-centers(:,1)';
dy = path(:,2)-centers(:,2)';
dz = path(:,3)-centers(:,3)';
dist = sqrt(dx.^2+dy.^2+dz.^2);            

hard = dist < radii';
n_hard = sum(hard(:));
penalty = n_hard * cfg.pen_hard;
soft = (dist >= radii') & (dist < radii' + cfg.safe_dist);
if any(soft(:))
    penalty = penalty + cfg.pen_soft * sum((radii'+cfg.safe_dist-dist).*soft,'all');
end

min_surf = min(dist - radii', [], 2);
min_clearance = min(min_surf);

v1 = path(2:end-1,:)-path(1:end-2,:);
v2 = path(3:end,:)-path(2:end-1,:);
n1 = sqrt(sum(v1.^2,2)); n2 = sqrt(sum(v2.^2,2));
cos_a = sum(v1.*v2,2)./(n1.*n2+eps);
cos_a = max(-1,min(1,cos_a));
curv = sum(acos(cos_a).^2);

axis_u = (goal-start)/norm(goal-start);
prog = path*axis_u(:);
backtrack = sum(max(0, prog(1:end-1)-prog(2:end)));

cost = cfg.w_len*total_len + cfg.w_obs*penalty + cfg.w_smooth*curv + cfg.w_mono*backtrack;

if nargout > 1
    info.path = path;
    info.cost = cost;
    info.len = total_len;
    info.curv = curv;
    info.backtrack = backtrack;
    info.n_hard = n_hard;
    info.feasible = (n_hard==0);
    info.min_clearance = min_clearance;
end
end
