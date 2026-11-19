function check_engineering

probs = eng_problems();
npass=0; nfail=0;
opts=optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxFunctionEvaluations',3e5,'StepTolerance',1e-12,'OptimalityTolerance',1e-10);

for ip=1:numel(probs)
    p=probs(ip);
    f=p.rawf(p.xref); g=p.rawg(p.xref);
    feas = max(g./p.cScale(:)) <= 1e-3;
    relErr=abs(f-p.fref)/abs(p.fref);
    fprintf('P%d %-18s reference point f=%.6f (Literature %.4f, relative deviation %.2e), max(g)=%.3g %s\n', ...
        ip,p.tag,f,p.fref,relErr,max(g./p.cScale(:)),passStr(feas));
    if feas, npass=npass+1; else, nfail=nfail+1; end
end

% --------
p=probs(1); best=inf;
for s=1:120
    x0=p.lb+rand(1,p.dim).*(p.ub-p.lb);
    try
        [x,~,ef]=fmincon(p.rawf,x0,[],[],[],[],p.lb,p.ub,@(x)deal(p.rawg(x),[]),opts);
        if ef>0 && max(p.rawg(x))<=1e-6, best=min(best,p.rawf(x)); end
    catch, end
end
ok=abs(best-1.724852)<2e-3; report('焊接梁 fmincon',best,1.724852,ok); npass=npass+ok; nfail=nfail+~ok;

% --------
p=probs(2); best=inf;
starts=[[.8125 .4375 42.1 176.6];[1 1 50 100];rand(40,4).*[2 2 60 200]+[.06 .06 5 5]];
for s=1:size(starts,1)
    try
        [x,~,ef]=fmincon(p.rawf,starts(s,:),[],[],[],[],p.lb,p.ub,@(x)deal(p.rawg(x),[]),opts);
        if ef>0 && max(p.rawg(x)./p.cScale(:))<=1e-6, best=min(best,p.rawf(x)); end
    catch, end
end
ok=abs(best-6059.7143)<1; report('pressure vessel fmincon',best,6059.7143,ok); npass=npass+ok; nfail=nfail+~ok;

% -------
p=probs(3); best=inf;
for s=1:120
    x0=p.lb+rand(1,p.dim).*(p.ub-p.lb);
    try
        [x,~,ef]=fmincon(p.rawf,x0,[],[],[],[],p.lb,p.ub,@(x)deal(p.rawg(x),[]),opts);
        if ef>0 && max(p.rawg(x))<=1e-6, best=min(best,p.rawf(x)); end
    catch, end
end
ok=abs(best-2994.4245)/2994.4245<5e-3; report('reducer fmincon',best,2994.4245,ok); npass=npass+ok; nfail=nfail+~ok;

% ---- ----
p=probs(4); best=inf;
for s=1:60
    x0=p.lb+rand(1,p.dim).*(p.ub-p.lb);
    try
        [x,~,ef]=fmincon(p.rawf,x0,[],[],[],[],p.lb,p.ub,@(x)deal(p.rawg(x),[]),opts);
        if ef>0 && max(p.rawg(x))<=1e-6, best=min(best,p.rawf(x)); end
    catch, end
end
ok=abs(best-263.8958)<2e-2; report('Three-bar truss fmincon',best,263.8958,ok); npass=npass+ok; nfail=nfail+~ok;

% -------
p=probs(5); best=inf; bx=[];
As_set=p.discreteSet{1};
for iA=1:numel(As_set)
  for bb=28:40
    for hh=5:0.002:10
      x=[As_set(iA) bb hh]; g=p.rawg(x);
      if max(g)<=1e-4
        f=p.rawf(x); if f<best, best=f; bx=x; end
      end
    end
  end
end
ok=abs(best-359.2080)<5e-2;
report('RCB Discrete enumeration',best,359.2080,ok);
fprintf('     RCB discrete optimal point x=[%.2f %d %.4f]\n',bx(1),bx(2),bx(3));
npass=npass+ok; nfail=nfail+~ok;
pcont=p; pcont.lb=[0 28 5]; pcont.ub=[10 40 10];
bestc=inf;
for s=1:60
    x0=pcont.lb+rand(1,3).*(pcont.ub-pcont.lb);
    try
        [x,~,ef]=fmincon(pcont.rawf,x0,[],[],[],[],pcont.lb,pcont.ub,@(x)deal(pcont.rawg(x),[]),opts);
        if ef>0 && max(pcont.rawg(x))<=1e-6, bestc=min(bestc,pcont.rawf(x)); end
    catch, end
end
fprintf('    RCB continuous relaxation optimal f=%.4f (document 359.2080 is a mixed discrete optimal, normal)\n',bestc);

% --------
p=probs(6); best=inf;
for s=1:120
    x0=p.lb+rand(1,p.dim).*(p.ub-p.lb);
    try
        [x,~,ef]=fmincon(p.rawf,x0,[],[],[],[],p.lb,p.ub,@(x)deal(p.rawg(x),[]),opts);
        if ef>0 && max(p.rawg(x))<=1e-6, best=min(best,p.rawf(x)); end
    catch, end
end
ok=abs(best-26.4995)<5e-3; report('Tubular column fmincon',best,26.4995,ok); npass=npass+ok; nfail=nfail+~ok;

fprintf('\n================Validation result：PASS %d，FAIL %d ================\n',npass,nfail);
if nfail>0, error('There are failed formula fingerprint verifications. Please fix eng_problems.m'); end
end

function s=passStr(b)
if b, s='[PASS]'; else, s='[FAIL]'; end
end
function report(tag,got,ref,ok)
fprintf('%-22s get %.6f，reference %.6f，%s\n',tag,got,ref,passStr(ok));
end
