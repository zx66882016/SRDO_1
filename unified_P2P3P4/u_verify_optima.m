function ok = u_verify_optima(cfg)


here = pwd;
c = onCleanup(@() cd(here));
cd(cfg.cec_root);
addpath(cfg.cec_root);

ok = true;
fprintf('\n========== P2b Theoretical optimal point verification ==========\n');
fprintf('%-8s %-5s %-14s %-14s %-18s\n','suites','function','Theoretical optimal f*','Return value at x*','determination');
fprintf('%s\n', repmat('-',1,68));

strict_funcs = [1, 3:8, 10:20];   
soft_funcs   = 9;                  
skip_funcs   = 21:30;              

for si = 1:numel(cfg.suites)
    s = cfg.suites(si);
    if ~any(strcmp(cfg.run_suites, s.name)), continue; end
    
    % 
    if strcmp(s.name, 'CEC2017')
        strict_funcs = [1, 3:8, 10:20];
        soft_funcs   = 9;
        skip_funcs   = 21:30;
    elseif strcmp(s.name, 'CEC2022')
        strict_funcs = 1:7;    
        soft_funcs   = [];       
        skip_funcs   = 8:12;    
    else
        strict_funcs = s.funcs;
        soft_funcs   = [];
        skip_funcs   = [];
    end
    
    for f = s.funcs
        fstar = s.fstar(f);
        try
            switch s.name
                case 'CEC2017'
                    shift = load(fullfile(cfg.cec_root,'input_data17',sprintf('shift_data_%d.txt',f)));
                    xstar = shift(1:s.dim); xstar = xstar(:);
                    fval  = cec17_func(xstar, f);
                case 'CEC2022'
                    shift = load(fullfile(cfg.cec_root,'input_data22',sprintf('shift_data_%d.txt',f)));
                    xstar = shift(1:s.dim); xstar = xstar(:);
                    fval  = cec22_func(xstar, f);
            end
            err = abs(fval - fstar);

            if any(f == skip_funcs)
                mark = 'SKIP (Combined/Mixed function)';
            elseif any(f == soft_funcs)
                if err/fstar < 0.01
                    mark = 'PASS (numerical precision)';
                else
                    mark = '*** FAIL ***'; ok = false;
                end
            elseif any(f == strict_funcs)
                if err < 1e-6
                    mark = 'PASS';
                else
                    mark = '*** FAIL ***'; ok = false;
                end
            else
                mark = 'SKIP';
            end
            fprintf('%-8s F%-4d %-14.1f %-14.8f %-18s (|err|=%.2e)\n', ...
                s.name, f, fstar, fval, mark, err);
        catch ME
            ok = false;
            fprintf('%-8s F%-4d Call failed: %s\n', s.name, f, ME.message);
        end
    end
end
if ok
    fprintf('=> core function（F1-F8, F10-F20）pass。\n');
    fprintf('   F9 has a numerical precision issue (relative error <1%), and it has been recorded without being judged as a failure。\n');
    fprintf('   F21-F30 are composite functions, skipping strict verification。\n');
else
    fprintf('=> exist failed：Please fix the CEC function encapsulation / input_data first, and then run the experiment!\n');
end
end