function [fmin, best_pos, Convergence_curve] = IPOP_CMAES_uni(~, ~, down, up, dim, fobj, opts, bo)


if nargin < 7, opts=struct(); end
if isscalar(down), down=down*ones(1,dim); end
if isscalar(up),   up=up*ones(1,dim);   end
sigma0 = 0.3*(up(1)-down(1));
lam_base = 4 + floor(3*log(dim));
if isfield(opts,'lambda0'), lam_base=opts.lambda0; end
max_restarts = 20;
if isfield(opts,'max_restarts'), max_restarts=opts.max_restarts; end

curve = [];
for restart = 0:max_restarts
    lambda = lam_base * 2^restart;
    mu = floor(lambda/2);
    w  = log(mu+0.5)-log(1:mu); w = w/sum(w); w=w(:);
    mueff = 1/sum(w.^2);
    cc   = (4+mueff/dim)/(dim+4+2*mueff/dim);
    cs   = (mueff+2)/(dim+mueff+5);
    c1   = 2/((dim+1.3)^2+mueff);
    cmu  = min(1-c1, 2*(mueff-2+1/mueff)/((dim+2)^2+mueff));
    damps= 1+2*max(0,sqrt((mueff-1)/(dim+1))-1)+cs;
    chiN = sqrt(dim)*(1-1/(4*dim)+1/(21*dim^2));

    xmean = (down+up)/2 + (up-down).*(rand(1,dim)-0.5);
    sigma = sigma0;
    pc=zeros(1,dim); ps=zeros(1,dim); B=eye(dim); D=ones(1,dim); C=eye(dim);
    gen=0; eigeneval=0;
    besthist=[]; arbest=inf;
    terminated=false;

    while ~terminated && bo.count < bo.maxFEs
        arx = zeros(lambda,dim); arfitness=inf(lambda,1);
        for k=1:lambda
            for attempt=1:10
                z=randn(dim,1);
                y=(B*diag(D)*z).';
                x=xmean+sigma*y;
                if all(x>=down) && all(x<=up), break; end
            end
            x=min(max(x,down),up);
            arx(k,:)=x;
            arfitness(k)=fobj(x);
        end
        [fitness, isort]=sort(arfitness);
        xold=xmean;
        xmean=w.'*arx(isort(1:mu),:);
        gen=gen+1;
        besthist(end+1)=fitness(1); 
        curve(end+1,1)=fitness(1);

        zmean = (1/sigma)*(xmean-xold)*B*diag(1./D)*B.'; 
        ps=(1-cs)*ps + sqrt(cs*(2-cs)*mueff)*(B*diag(D)*zmean.').';
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*gen))/chiN < 1.4+2/(dim+1);
        pc=(1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)*(xmean-xold)/sigma;
        artmp=(arx(isort(1:mu),:)-repmat(xold,mu,1))/sigma;
        C=(1-c1-cmu)*C + c1*(pc.'*pc + (1-hsig)*cc*(2-cc)*C) ...
            + cmu*artmp.'*diag(w)*artmp;
        sigma=sigma*exp((cs/damps)*(norm(ps)/chiN-1));

        if gen-eigeneval > 1/(c1+cmu)/dim/10
            eigeneval=gen;
            C=triu(C)+triu(C,1).';
            [B,DD]=eig(C); D=sqrt(max(real(diag(DD)),1e-20)).';
            C=real(C);
        end

        % ----
        win=ceil(10+30*dim/lambda);
        if numel(besthist) > win
            recent=besthist(end-win+1:end);
            if max(recent)-min(recent) <= 1e-12, terminated=true; end    
        end
        if sigma*max([abs(pc), max(D)]) < 1e-11*sigma0, terminated=true; end 
        if max(D)/min(D) > 1e7, terminated=true; end                       
        if all(xmean==xold), terminated=true; end                          
    end
    if bo.count >= bo.maxFEs, break; end
end
fmin=bo.best; best_pos=bo.bestX;
Convergence_curve=curve(:).';
end
