function [fmin, best_pos, Convergence_curve, localInfo] = SRDO_uni(N, MaxIter, down, up, dim, fobj, opts, bo)

if nargin < 7 || isempty(opts), opts = struct(); end
if nargin < 8, bo = []; end
opts = set_defaults(opts);

alpha0 = 0.95; p0 = 0.85; elite_ratio = 0.25;
mutation_rate = 0.15; diversity_threshold = 0.02;

if isscalar(down), down = down*ones(1,dim); end
if isscalar(up),   up   = up*ones(1,dim);   end

% Mixed initialization
X = stable_initialization(N, dim, up, down);
Fitness = zeros(1,N);
for i = 1:N, Fitness(i) = fobj(X(i,:)); end

[~, idx] = sort(Fitness);
elite_num = max(3, round(N*elite_ratio));
elite_pop = X(idx(1:elite_num),:);
global_best = elite_pop(1,:);
global_best_fit = Fitness(idx(1));
Hunger = ones(N,1);
Convergence_curve = inf(1,MaxIter);
last_improvement = 0;
last_restart_fe = 0;
last_pattern_fe = 0;

% Route B: Cumulative budget for local fine-tuning (only enabled under budgetary scope)
localUsed = 0; localCalls = 0; localWins = 0;
patternFails = 0; cmaDone = false;
localCap = 0;
if opts.useLocalCMA && isa(bo,'BudgetObjective')
    localCap = opts.localFrac*bo.maxFEs;
end
localLambda = 4 + floor(3*log(dim));

% Ablation control: Random trigger points, used to verify the necessity of
% fixed schedule scheduling 
randTrigPoint = inf;
if opts.useLocalCMA && isfield(opts,'localRandTrigger') && opts.localRandTrigger
    randTrigPoint = 0.6 + 0.3*rand;
end

for t = 1:MaxIter
    % ---- Progress: Prioritize based on actual FEs, and revert to the iteration metric if no counter is available ----
    if isa(bo,'BudgetObjective')
        progress = bo.count / bo.maxFEs;
        fe_now   = bo.count;
    else
        progress = t/MaxIter; fe_now = N*t;
    end

    % ----Parametric curve----
    if strcmp(opts.mode,'high')
        sig = 1/(1+exp(-10*(progress-0.5)));
        alpha   = 0.1 + 0.85*(1-sig);
        p       = 0.2 + 0.65*(1-sig);
        lens_k  = 0.5 + 0.4*sig;
        levy_factor = 1.0 - 0.5*progress;
        sched_sig = sig;
    else
        switch opts.scheme
            case 1 % Exponential decay (fast mode)
                alpha = alpha0*exp(-4*progress^2);
                p     = p0*(1-progress^0.8);
                levy_factor = 1.0 - 0.5*progress;
                lens_k = 0.5 + 0.4*progress;
            case 2 % fixed parameter
                alpha=0.4; p=0.6; levy_factor=0.7; lens_k=0.7;
            case 3 % S-shaped
                alpha = 0.1 + 0.85/(1+exp(10*(progress-0.5)));
                p     = 0.1 + 0.75/(1+exp(8*(progress-0.5)));
                levy_factor = 0.2 + 0.8/(1+exp(6*(progress-0.5)));
                lens_k = 0.3 + 0.6/(1+exp(-8*(progress-0.5)));
            case 4 % Pure linear
                alpha = 0.9 - 0.8*progress;
                p     = 0.8 - 0.6*progress;
                levy_factor = 1.0 - 0.6*progress;
                lens_k = 0.5 + 0.5*progress;
        end
        sched_sig = [];
    end
    mutation_factor = 0.9*(1-progress^0.4);

    % Stagnation Enhancement
    stall_frac = 0.1;
    if t - last_improvement > stall_frac*MaxIter
        alpha = min(1.0, alpha*1.2);
        p     = min(0.95, p*1.1);
    end

    % ---- hunger factor----
    best_fitness = min(Fitness); worst_fitness = max(Fitness);
    for i = 1:N
        if Fitness(i) ~= best_fitness
            nf = (Fitness(i)-best_fitness)/(worst_fitness-best_fitness+eps);
            Hunger(i) = Hunger(i) + nf^1.5;
        end
    end
    mh = max(Hunger); if mh>0, Hunger = Hunger/mh; end

    % ---- Individual update ----
    for i = 1:N
        rnd = rand;
        if rnd < opts.split1
            % ===== Strategy A: Long-range Levy=====
            if opts.repl == 1
                step = gaussian_step(dim, progress);
            else
                step = stable_levy_step(dim, progress, Hunger(i), levy_factor);
            end
            if strcmp(opts.mode,'high')
                if rand < p, leader = elite_pop(randi(elite_num),:);
                    X_new = X(i,:) + alpha*step.*(leader-X(i,:));
                else
                    X_new = X(i,:) + alpha*step.*(up-down);
                end
            else
                if rand < (1-progress), X_new = X(i,:) + alpha*step.*(up-down);
                else, leader = elite_pop(randi(elite_num),:);
                    X_new = X(i,:) + alpha*step.*(leader-X(i,:));
                end
            end
        elseif rnd < opts.split2
            % ===== Strategy B: Differential Perturbation=====
            if opts.repl == 2
                X_new = X(i,:) + alpha*0.8*(up-down).*(2*rand(1,dim)-1);
            elseif ~opts.useDiff
                step = stable_levy_step(dim, progress, Hunger(i), levy_factor);
                X_new = X(i,:) + alpha*0.8*step.*(up-down);
            else
                r1=randi(N); r2=randi(N);
                while r1==i || r2==i || r1==r2
                    r1=randi(N); r2=randi(N);
                end
                if opts.useHunger
                    F = 0.5 + 0.4*Hunger(i)*(1-progress);
                else
                    F = 0.5;
                end
                X_new = X(i,:) + F*(global_best-X(i,:)) + F*(X(r1,:)-X(r2,:));
            end
        else
            % ===== Strategy C: Auditory Lens =====
            leader = elite_pop(randi(elite_num),:);
            if opts.repl == 3
                fixed_step = 0.05*(up-down);
                X_new = X(i,:) + fixed_step.*randn(1,dim) + 0.3*(leader-X(i,:));
            elseif ~opts.useLens
                step = stable_levy_step(dim, progress, Hunger(i), levy_factor);
                X_new = X(i,:) + alpha*0.6*step.*(leader-X(i,:));
            else
                X_new = auditory_lens_imaging(X(i,:), leader, down, up, lens_k);
            end
        end

        X_new = stable_bounds(X_new, down, up);

        % Gaussian micro-mutation
        if rand < mutation_rate*mutation_factor
            ms = (up-down).*(0.01 + 0.05*(1-progress)*rand(1,dim));
            X_new = X_new + ms.*randn(1,dim);
            X_new = stable_bounds(X_new, down, up);
        end

        Fnew = fobj(X_new);
        if Fnew < Fitness(i)
            X(i,:) = X_new; Fitness(i) = Fnew;
            if Fnew < global_best_fit
                global_best = X_new; global_best_fit = Fnew; last_improvement = t;
            end
        end
    end

    [~, idx] = sort(Fitness);
    elite_pop = X(idx(1:elite_num),:);
    Convergence_curve(t) = global_best_fit;

    % ---- Diversity maintenance ----
    if opts.useDiversity
        diversity = mean(std(X));
        need = (isa(bo,'BudgetObjective') && (fe_now-last_restart_fe >= 10*N)) ...
               || (~isa(bo,'BudgetObjective') && mod(t,10)==0) ...
               || diversity < diversity_threshold;
        if need
            worst_num = max(1, round(N*0.1));
            worst_idx = idx(end-worst_num+1:end);
            X(worst_idx,:) = stable_initialization(worst_num, dim, up, down);
            for k = 1:worst_num, Fitness(worst_idx(k)) = fobj(X(worst_idx(k),:)); end
            [~, idx] = sort(Fitness); elite_pop = X(idx(1:elite_num),:);
            last_restart_fe = fe_now;
        end
    end

    % ---- Final finishing (Option A: stall-triggered local finishing)----
    % stall_gen: The number of consecutive generations without improvement from the last global optimal improvement to the current iteration
    stall_gen = t - last_improvement;
    trigThresh = opts.localTrigger;
    if isfinite(randTrigPoint), trigThresh = randTrigPoint; end
    inLocalStage = opts.useLocalCMA && localCap > 0 && progress > trigThresh;
    % Local fine-tuning trigger condition: Entering the local phase, having not been executed before, and achieving no improvement continuously up to the threshold
    localTriggered = inLocalStage && ~cmaDone && stall_gen >= opts.localStallTrigger;

   % Mode Search: Maintain the current triggering rhythm of SRDO (diagnosis has proven that prematurely shutting down mode search can worsen the budget constraints of % F12 by two orders of magnitude); local refinement only consumes one's own maximum budget, with both coexisting
    if opts.usePattern
        gate = progress > 0.7 && (isa(bo,'BudgetObjective') ...
               && (fe_now-last_pattern_fe >= 8*N) ...
               || (~isa(bo,'BudgetObjective') && mod(t,8)==0));
        if gate
            step_size = 0.02*(up-down)*(1-progress);
            fitBefore = global_best_fit;
            [rs, rf] = pattern_search(global_best, down, up, fobj, step_size, 8);
            if rf < fitBefore
                global_best = rs; global_best_fit = rf;
                elite_pop(1,:) = rs; Fitness(idx(1)) = rf;
                Convergence_curve(t) = global_best_fit;
                patternFails = 0;
            else
                patternFails = patternFails + 1;
            end
            last_pattern_fe = fe_now;
        end
    end

    % Perform local fine-tuning (execute once when the trigger condition is met, and do not trigger again after cmaDone=true)
    if localTriggered
        % auto: If the remaining budget is sufficient for CMA to be used for more than 300 generations, then use CMA; otherwise, use Nelder-Mead
        solver = opts.localSolver;
        if strcmpi(solver,'auto')
            if isa(bo,'BudgetObjective'), remainFE = bo.maxFEs - bo.count;
            else, remainFE = N*(MaxIter-t); end
            if remainFE >= 300*localLambda, solver='cma'; else, solver='nm'; end
        end
        used = 0;
        switch solver
            case 'nm'
                [rs, rf, u1] = local_nm_refine(global_best, global_best_fit, ...
                    down, up, dim, fobj, bo, localCap-localUsed, ...
                    opts.localNMstep*(up-down), opts.localStall);
                used=u1;
                if isfield(opts,'debugLocal') && opts.debugLocal
                    fprintf('[localNM] feNow=%d used=%d stallGen=%d fBefore=%.6g fAfter=%.6g\n',...
                        fe_now,used,stall_gen,global_best_fit,rf);
                end
            otherwise
                k = max(1,opts.localMultiK);
                if k<=1
                    [rs, rf, u1] = local_cma_refine(global_best, global_best_fit, ...
                        down, up, dim, fobj, bo, localCap-localUsed, ...
                        opts.localSigma*(up-down), opts.localStall);
                    used=u1;
                else
                    % Multi-start: Select k non-overlapping elites from the current population, refine them independently, and allocate the budget equally
                    [~,six]=sort(Fitness);
                    starts=zeros(k,dim); sf=zeros(k,1); m=0;
                    for q=1:numel(six)
                        xq=X(six(q),:);
                        if m==0 || min(max(abs(starts(1:m,:)-repmat(xq,m,1))./max(up-down,eps),[],2)) > 1e-8
                            m=m+1; starts(m,:)=xq; sf(m)=Fitness(six(q));
                        end
                        if m>=k, break; end
                    end
                    perCap=max(floor((localCap-localUsed)/m),localLambda);
                    rs=global_best; rf=global_best_fit;
                    for q=1:m
                        [rq,fq,uq]=local_cma_refine(starts(q,:),sf(q), ...
                            down, up, dim, fobj, bo, perCap, ...
                            opts.localSigma*(up-down), opts.localStall);
                        used=used+uq;
                        if isfield(opts,'debugLocal') && opts.debugLocal
                            fprintf('[localCMA#%d/%d] feNow=%d used=%d fStart=%.6g fAfter=%.6g\n',...
                                q,m,fe_now,uq,sf(q),fq);
                        end
                        if fq<rf, rs=rq; rf=fq; end
                        if isa(bo,'BudgetObjective') && bo.count>=bo.maxFEs, break; end
                    end
                end
        end
        localUsed = localUsed + used; localCalls = localCalls+1; cmaDone = true;
        winNow = rf < global_best_fit;
        if winNow
            localWins = localWins+1;
            [~, idx] = sort(Fitness);
            global_best = rs; global_best_fit = rf;
            elite_pop(1,:) = rs; Fitness(idx(1)) = rf;
            Convergence_curve(t) = global_best_fit;
        end
    end
end

% ---- Final refinement (only reached when there is no budget counter or the main loop has not exhausted its budget) ----
if opts.useLocalCMA && localCap > 0
    remLocal = localCap - localUsed;
    if remLocal >= localLambda
        switch opts.localSolver
            case 'nm'
                [bp, fm] = local_nm_refine(global_best, global_best_fit, ...
                    down, up, dim, fobj, bo, remLocal, ...
                    0.5*opts.localNMstep*(up-down), opts.localStall+4);
            otherwise
                [bp, fm] = local_cma_refine(global_best, global_best_fit, ...
                    down, up, dim, fobj, bo, remLocal, ...
                    0.5*opts.localSigma*(up-down), opts.localStall+4);
        end
        if fm < global_best_fit, global_best = bp; global_best_fit = fm; end
    end
elseif opts.usePattern
    [bp, fm] = pattern_search(global_best, down, up, fobj, 0.001*(up-down), 12);
    if fm < global_best_fit, global_best = bp; global_best_fit = fm; end
end
fmin = global_best_fit; best_pos = global_best;
Convergence_curve(isinf(Convergence_curve)) = fmin;
Convergence_curve(end) = fmin;
if nargout>3
    localInfo = struct('cap',localCap,'used',localUsed,'calls',localCalls,'wins',localWins);
end
end

% ============================default parameter ============================
function opts = set_defaults(opts)
if ~isfield(opts,'mode'),   opts.mode='fast'; end
if ~isfield(opts,'scheme'), opts.scheme=1;   end
if ~isfield(opts,'useDiff'),      opts.useDiff=true;      end
if ~isfield(opts,'useLens'),      opts.useLens=true;      end
if ~isfield(opts,'useHunger'),    opts.useHunger=true;    end
if ~isfield(opts,'useDiversity'), opts.useDiversity=true; end
if ~isfield(opts,'usePattern'),   opts.usePattern=true;   end
if ~isfield(opts,'repl'),         opts.repl=0;            end
% Route B Candidate 1/2: Final local fine-tuning (local CMA-ES or Nelder-Mead), disabled by default, does not affect the current version
if ~isfield(opts,'useLocalCMA'),  opts.useLocalCMA=false; end
if ~isfield(opts,'localSolver'),  opts.localSolver='cma'; end  % 'cma'、'nm' or'auto'
if ~isfield(opts,'localMultiK'),  opts.localMultiK=1;     end  % CMA multi-starting point elite number
if ~isfield(opts,'localFrac'),    opts.localFrac=0.18;    end  % The proportion of the total budget allocated for local fine-tuning to the overall budget
if ~isfield(opts,'localSigma'),   opts.localSigma=0.02;   end  % CMA: sigma0=localSigma*(ub-lb)
if ~isfield(opts,'localNMstep'),  opts.localNMstep=0.02;  end  % NM
if ~isfield(opts,'localTrigger'), opts.localTrigger=0.8;  end  % Local fine-tuning triggers progress
if ~isfield(opts,'localStall'),   opts.localStall=200;    end  % Long-range fine-tuning continuously stagnates, and algebra is about to exit
% Option A: Add stall-triggered parameter
if ~isfield(opts,'localStallTrigger'), opts.localStallTrigger=50; end  % How many generations of continuous no improvement will trigger it
if ~isfield(opts,'localRandTrigger'),  opts.localRandTrigger=false; end % Ablation: Random trigger point comparison
% Sensitivity to segmentation points: Lévy segmentation split1, differential segmentation split2
if ~isfield(opts,'split1'), opts.split1=0.35; end
if ~isfield(opts,'split2'), opts.split2=0.70; end
end

% ============================ helper function ============================
function X = stable_initialization(N, dim, up, down)
X = zeros(N,dim); chaos_num = round(N*0.6);
for i = 1:chaos_num
    chaos = rand(1,dim);
    for k = 1:8, chaos = tent_map(chaos); end
    X(i,:) = down + (up-down).*chaos;
end
if N-chaos_num > 0
    X(chaos_num+1:end,:) = down + (up-down).*rand(N-chaos_num,dim);
end
end

function seq = tent_map(seq)
a = 0.7;
for j = 1:length(seq)
    if seq(j) < a, seq(j) = seq(j)/a;
    else, seq(j) = (1-seq(j))/(1-a); end
end
end

function vec = stable_levy_step(d, progress, hunger, levy_factor)
beta = 1.9 - 0.7*progress + 0.1*hunger;
beta = max(1.2, min(beta,1.95));
sigma = (gamma(1+beta)*sin(pi*beta/2)/ ...
        (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = randn(1,d)*sigma; v = randn(1,d); v(v==0)=1e-10;
step = u./abs(v).^(1/beta);
scale = 0.01 + 0.2*(1-progress)^1.5;
vec = scale*step*levy_factor;
end

function vec = gaussian_step(d, progress)
scale = 0.05 + 0.15*(1-progress)^1.5;
vec = scale*randn(1,d);
end

function X_new = auditory_lens_imaging(X, elite, lb, ub, k)
dlb = min(X,elite); dub = max(X,elite);
if any(dub-dlb < 1e-12)
    dlb = max(dlb-0.5,lb); dub = min(dub+0.5,ub);
end
dlb = max(dlb,lb); dub = min(dub,ub);
center = (dlb+dub)/2;
lens_point = center + (center-X)/k;
X_new = 0.6*lens_point + 0.4*elite;
end

function X = stable_bounds(X,Lb,Ub)
below = X<Lb; X(below) = 2*Lb(below)-X(below);
above = X>Ub; X(above) = 2*Ub(above)-X(above);
X = max(min(X,Ub),Lb);
X(X<Lb)=Lb(X<Lb); X(X>Ub)=Ub(X>Ub);
end

function [ps,pf] = pattern_search(sol,lb,ub,fobj,step_size,max_iter)
if nargin<6, max_iter=12; end
n=numel(sol); cur=sol; cur_fit=fobj(sol);
for it=1:max_iter
    improved=false;
    for i=1:n
        c=cur; c(i)=min(cur(i)+step_size(i),ub(i)); cf=fobj(c);
        if cf<cur_fit, cur=c; cur_fit=cf; improved=true; continue; end
        c=cur; c(i)=max(cur(i)-step_size(i),lb(i)); cf=fobj(c);
        if cf<cur_fit, cur=c; cur_fit=cf; improved=true; end
    end
    if ~improved
        step_size=0.5*step_size;
        if max(step_size)<1e-10, break; end
    end
end
ps=cur; pf=cur_fit;
end

function [bestx,bestf,used] = local_cma_refine(x0,x0fit,lb,ub,n,fobj,bo,feCap,sigma0,stallGen)
% LOCAL_CMA_REFINE -CMA-ES 
if nargin<10, stallGen=8; end
lb=lb(:); ub=ub(:);
lambda = 4 + floor(3*log(n));
mu     = floor(lambda/2);
w      = log(mu+0.5)-log(1:mu); w = w/sum(w); w=w(:);
mueff  = 1/sum(w.^2);
cc     = (4+mueff/n)/(n+4+2*mueff/n);
cs     = (mueff+2)/(n+mueff+5);
c1     = 2/((n+1.3)^2+mueff);
cmu    = min(1-c1, 2*(mueff-2+1/mueff)/((n+2)^2+mueff));
damps  = 1 + 2*max(0,sqrt((mueff-1)/(n+1))-1) + cs;
chiN   = sqrt(n)*(1-1/(4*n)+1/(21*n^2));

xmean = x0(:);
sigma = mean(sigma0(:));
pc=zeros(n,1); ps=zeros(n,1);
B=eye(n); D=ones(n,1); C=eye(n); invsqrtC=eye(n);
bestx=xmean; bestf=x0fit;
used=0; gen=0; lastImprove=0; stop=false;
rangeMean = mean(ub-lb);

while ~stop
    gen = gen+1;
    arx = zeros(n,lambda); arfit = inf(lambda,1);
    for k=1:lambda
        for attempt=1:10
            z=randn(n,1);
            x=xmean + sigma*(B*(D.*z));
            if all(x>=lb) && all(x<=ub), break; end
        end
        x=min(max(x,lb),ub);
        try
            fv = fobj(x.');
        catch ME
            if strcmp(ME.identifier,'BudgetObjective:Exhausted')
                stop=true; break;
            else
                rethrow(ME);
            end
        end
        used=used+1; arx(:,k)=x; arfit(k)=fv;
        if isfinite(fv) && fv<bestf
            bestf=fv; bestx=x; lastImprove=gen;
        end
        if used>=feCap, stop=true; break; end
    end
    if stop || any(isinf(arfit))
        keep=~isinf(arfit);
        if sum(keep)<mu, break; end
        arx=arx(:,keep); arfit=arfit(keep);
        if numel(arfit)<mu, break; end
    end
    [~,fi]=sort(arfit);
    sel=fi(1:mu);
    xold=xmean;
    xmean = arx(:,sel)*w;
    ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)*invsqrtC*((xmean-xold)/sigma);
    hsig = norm(ps)/sqrt(1-(1-cs)^(2*gen))/chiN < 1.4 + 2/(n+1);
    pc   = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(xmean-xold)/sigma;
    artmp = (arx(:,sel)-repmat(xold,1,mu))/sigma;
    C = (1-c1-cmu)*C ...
        + c1*(pc*pc.' + (1-hsig)*cc*(2-cc)*C) ...
        + cmu*artmp*diag(w)*artmp.';
    C = triu(C)+triu(C,1).';
    C = (C+C.')/2;
    C = C + eye(n)*1e-12*trace(C)/n;
    sigma = sigma*exp((cs/damps)*(norm(ps)/chiN - 1));
    [B,DD]=eig(C);
    dd=real(diag(DD)); dd=max(dd,1e-20);
    B=real(B);
    D=sqrt(dd); invsqrtC=B*diag(1./D)*B.';
    if gen-lastImprove>=stallGen, stop=true; end
    if sigma*max(D) < 1e-12*rangeMean, stop=true; end
    if max(D)>1e7*min(D), stop=true; end
    if isa(bo,'BudgetObjective') && bo.count>=bo.maxFEs, stop=true; end
end
bestx=bestx.';
end

function [bestx,bestf,used] = local_nm_refine(x0,x0fit,lb,ub,n,fobj,bo,feCap,step0,stallIter)
% LOCAL_NM_REFINE   Nelder-Mead 
if nargin<10, stallIter=200; end
lb=lb(:); ub=ub(:); x0=x0(:); step0=step0(:);
a=1; g=2; r=0.5; s=0.5;

S=repmat(x0,1,n+1);
for j=1:n
    v=x0; v(j)=min(x0(j)+step0(j),ub(j));
    if v(j)==x0(j), v(j)=max(x0(j)-step0(j),lb(j)); end
    S(:,j+1)=v;
end
F=inf(1,n+1); F(1)=x0fit; used=0; stop=false;
for j=2:n+1
    [F(j),ok]=nm_eval(fobj,S(:,j).');
    if ~ok, stop=true; break; end
    used=used+1;
    if used>=feCap, stop=true; break; end
end
bestx=x0; bestf=x0fit;
gen=0; lastImprove=0; rangeMean=mean(ub-lb);

while ~stop
    gen=gen+1;
    [F,ix]=sort(F); S=S(:,ix);
    if F(1)<bestf-1e-12*(abs(bestf)+1), bestf=F(1); bestx=S(:,1); lastImprove=gen; end
    if gen-lastImprove>=stallIter, stop=true; end
    diam=max(max(abs(S(:,2:n+1)-repmat(S(:,1),1,n)),[],2));
    if diam<1e-12*rangeMean, stop=true; end
    if used>=feCap, stop=true; end
    if isa(bo,'BudgetObjective') && bo.count>=bo.maxFEs, stop=true; end
    if stop, break; end

    cen=mean(S(:,1:n),2);
    xr=min(max(cen+a*(cen-S(:,n+1)),lb),ub);
    [fr,ok]=nm_eval(fobj,xr.'); if ~ok, break; end
    used=used+1;
    if fr<F(1)
        xe=min(max(cen+g*(xr-cen),lb),ub);
        [fe,ok2]=nm_eval(fobj,xe.');
        if ~ok2, if fr<bestf, bestf=fr; bestx=xr; end; break; end
        used=used+1;
        if fe<fr, S(:,n+1)=xe; F(n+1)=fe; else, S(:,n+1)=xr; F(n+1)=fr; end
    elseif fr<F(n)
        S(:,n+1)=xr; F(n+1)=fr;
    else
        if fr<F(n+1)
            xc=min(max(cen+r*(xr-cen),lb),ub);
        else
            xc=min(max(cen-r*(cen-S(:,n+1)),lb),ub);
        end
        [fc,ok]=nm_eval(fobj,xc.'); if ~ok, break; end
        used=used+1;
        if fc<min(fr,F(n+1))
            S(:,n+1)=xc; F(n+1)=fc;
        else
            for j=2:n+1
                S(:,j)=min(max(S(:,1)+s*(S(:,j)-S(:,1)),lb),ub);
                [F(j),ok]=nm_eval(fobj,S(:,j).'); if ~ok, stop=true; break; end
                used=used+1;
                if used>=feCap, stop=true; break; end
            end
            if stop, break; end
        end
    end
    if used>=feCap, stop=true; end
end
[~,ix]=sort(F);
if F(ix(1))<bestf, bestf=F(ix(1)); bestx=S(:,ix(1)); end
bestx=bestx.';
end

function [fv,ok]=nm_eval(fobj,xrow)
ok=true;
try
    fv=fobj(xrow);
catch ME
    if strcmp(ME.identifier,'BudgetObjective:Exhausted')
        fv=inf; ok=false;
    else
        rethrow(ME);
    end
end
end