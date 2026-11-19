function [fmin, best_pos, Convergence_curve] = LSHADE_uni(~, ~, down, up, dim, fobj, opts, bo)


if nargin < 7, opts = struct(); end
H = 6; pbest_frac = 0.1; Nmin = 4;
if isfield(opts,'H'), H=opts.H; end
if isfield(opts,'p'), pbest_frac=opts.p; end
if isscalar(down), down=down*ones(1,dim); end
if isscalar(up),   up=up*ones(1,dim); end

Ninit = max(10, round(20*log(dim)));
X = down + (up-down).*rand(Ninit,dim);
fit = zeros(Ninit,1);
for i=1:Ninit, fit(i)=fobj(X(i,:)); end
N = Ninit;
MCR = 0.5*ones(H,1); MF = 0.5*ones(H,1);
k = 1; A = [];
curve = [];

while bo.count < bo.maxFEs
    SCR=[]; SF=[]; df=[];
    [~, I] = sort(fit);
    pmin = max(2/N, pbest_frac);
    u = X; fu = fit;
    for i=1:N
        r = randi(H);
        CR = MCR(r) + 0.1*randn;
        CR = min(1, max(0, CR));
        if isnan(MCR(r)), CR = 1; end          
        F = MF(r) + 0.1*tan(pi*(rand-0.5));
        while F <= 0, F = MF(r) + 0.1*tan(pi*(rand-0.5)); end
        F = min(1, F);
        p = unifrnd(2/N, pmin);
        pnum = max(2, round(p*N));
        pbest = I(randi(pnum));
        r1 = randi(N); while r1==i, r1=randi(N); end
        SA = [X; A];
        r2 = randi(size(SA,1));
        while r2==i || r2==r1, r2=randi(size(SA,1)); end
        v = X(i,:) + F*(X(pbest,:)-X(i,:)) + F*(X(r1,:)-SA(r2,:));
    
        above = v>up;  v(above) = (up(above)+X(i,above))/2;
        below = v<down; v(below) = (down(below)+X(i,below))/2;
        v = min(max(v,down),up);
 
        jr = randi(dim); mask = rand(1,dim) < CR; mask(jr)=true;
        ui = X(i,:); ui(mask)=v(mask);
        fi = fobj(ui);
        if fi < fit(i)
            if ~isempty(A), A=[A; X(i,:)]; else, A=X(i,:); end
            SCR=[SCR; CR]; SF=[SF; F]; df=[df; fit(i)-fi];
            u(i,:)=ui; fu(i)=fi;
        end
    end
    X=u; fit=fu;

    while size(A,1) > N, A(randi(size(A,1)),:)=[]; end

    if ~isempty(SF)
        w = df/sum(df);
        if all(SCR==0) || any(isnan(MCR))
            MCR(k)=NaN;
        else
            MCR(k)=sum(w.*SCR);
        end
        MF(k) = sum(w.*SF.^2)/sum(w.*SF);
        k=mod(k,H)+1;
    end
    curve(end+1,1)=min(fit); 

    Nnew = max(Nmin, round(Nmin + (Ninit-Nmin)*(bo.maxFEs-bo.count)/bo.maxFEs));
    if Nnew < N
        [~, Iw] = sort(fit,'descend');
        keep = sort(Iw(1:Nnew));
        X=X(keep,:); fit=fit(keep); N=Nnew;
    end
end
[fmin, ib] = min(fit); best_pos = X(ib,:);
Convergence_curve = curve(:).';
end
