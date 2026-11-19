function out = uav_astar3d(sc, cfg)


if nargin < 2, cfg = uav_config(); end
res = cfg.astar_res;
edges = 0:res:cfg.space;
ng = numel(edges);                 
[GX,GY,GZ] = ndgrid(edges,edges,edges);
gridPts = [GX(:),GY(:),GZ(:)];    % ng^3 × 3

% --------
occupied = false(ng,ng,ng);
for k = 1:size(sc.obs,1)
    d = sqrt(sum((gridPts-sc.obs(k,1:3)).^2,2));
    hit = reshape(d < sc.obs(k,4)+res/2, ng,ng,ng);
    occupied = occupied | hit;
end

% 
[~,sx] = min(abs(edges-sc.start(1))); [~,sy] = min(abs(edges-sc.start(2))); [~,sz] = min(abs(edges-sc.start(3)));
[~,gx_] = min(abs(edges-sc.goal(1))); [~,gy] = min(abs(edges-sc.goal(2))); [~,gz] = min(abs(edges-sc.goal(3)));
startNode = sub2ind([ng ng ng],sx,sy,sz);
goalNode  = sub2ind([ng ng ng],gx_,gy,gz);
occupied(sx,sy,sz)=false; occupied(gx_,gy,gz)=false;

% ---- ----
[dx,dy,dz] = ndgrid(-1:1,-1:1,-1:1);
dn = [dx(:),dy(:),dz(:)]; dn(14,:)=[];   
stepCost = sqrt(sum(dn.^2,2))*res;

gScore = inf(ng^3,1);
cameFrom = zeros(ng^3,1);
closed = false(ng^3,1);
gScore(startNode)=0;
heap = [];                         % [f, node]
heap = heapPush(heap,[h_(startNode),startNode]);
found = false;
t0=tic;
while ~isempty(heap)
    [item,heap]=heapPop(heap);
    cur=item(2);
    if closed(cur), continue; end
    closed(cur)=true;
    if cur==goalNode, found=true; break; end
    [ci,cj,ck]=ind2sub([ng ng ng],cur);
    for q=1:size(dn,1)
        ni=ci+dn(q,1); nj=cj+dn(q,2); nk=ck+dn(q,3);
        if ni<1||nj<1||nk<1||ni>ng||nj>ng||nk>ng, continue; end
        if occupied(ni,nj,nk), continue; end
        nb=sub2ind([ng ng ng],ni,nj,nk);
        if closed(nb), continue; end
        ng_=gScore(cur)+stepCost(q);
        if ng_<gScore(nb)
            gScore(nb)=ng_; cameFrom(nb)=cur;
            heap=heapPush(heap,[ng_+h_(nb),nb]);
        end
    end
end
searchTime=toc(t0);

out=struct('found',found,'res',res,'searchTime',searchTime);
if ~found, return; end

% --------
node=goalNode; pathIdx=node;
while node~=startNode
    node=cameFrom(node);
    if node==0, found=false; out.found=false; return; end
    pathIdx=[node;pathIdx]; 
end
[pi,pj,pk]=ind2sub([ng ng ng],pathIdx);
gridPath=[edges(pi)',edges(pj)',edges(pk)'];
gridPath(1,:)=sc.start; gridPath(end,:)=sc.goal;

% 
rawDense = resampleByArc(gridPath,cfg.n_sample);
out.raw = uav_eval_path(rawDense,sc.obs,cfg);
out.raw.gridPath = gridPath;

% ---- ----
t0=tic;
sp = splineSmooth(gridPath,cfg.n_sample);
out.spline = uav_eval_path(sp,sc.obs,cfg);
out.splineTime=toc(t0);

% ---------------- ----------------
function h = h_(node)
    [ii,jj,kk]=ind2sub([ng ng ng],node);
    h=norm([edges(ii)-edges(gx_),edges(jj)-edges(gy),edges(kk)-edges(gz)]);
end
end

% ====================================
function p = resampleByArc(path,n)
% 
seg=[0;cumsum(sqrt(sum(diff(path).^2,2)))];
if seg(end)==0, p=repmat(path(end,:),n,1); return; end
q=linspace(0,seg(end),n);
p=[interp1(seg,path(:,1),q,'linear')', ...
  interp1(seg,path(:,2),q,'linear')', ...
  interp1(seg,path(:,3),q,'linear')'];
end

function p = splineSmooth(path,n)
% 
seg=[0;cumsum(sqrt(sum(diff(path).^2,2)))];
q=linspace(0,seg(end),n);
p=[interp1(seg,path(:,1),q,'spline')', ...
  interp1(seg,path(:,2),q,'spline')', ...
  interp1(seg,path(:,3),q,'spline')'];
p(1,:)=path(1,:); p(end,:)=path(end,:);
end

function heap=heapPush(heap,item)
heap=[heap;item];
i=size(heap,1);
while i>1
    p=floor(i/2);
    if heap(p,1)<=heap(i,1), break; end
    heap([p i],:)=heap([i p],:); i=p;
end
end

function [item,heap]=heapPop(heap)
item=heap(1,:);
heap(1,:)=heap(end,:); heap(end,:)=[];
i=1; m=size(heap,1);
while true
    l=2*i; r=2*i+1; s=i;
    if l<=m && heap(l,1)<heap(s,1), s=l; end
    if r<=m && heap(r,1)<heap(s,1), s=r; end
    if s==i, break; end
    heap([s i],:)=heap([i s],:); i=s;
end
end
