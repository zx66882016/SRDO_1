classdef BudgetObjective < handle
    properties
        rawFcn
        maxFEs        = inf
        count         = 0
        best          = inf
        bestX         = []
        logFEs              
        logBest             
        logIdx        = 1
        cpFEs
        cpBest
        cpIdx         = 1
        exhausted     = false
        nOutBudget    = 0     
    end
    methods
        function obj = BudgetObjective(rawFcn, maxFEs, checkpoints, nLog)
            if nargin < 4, nLog = 200; end
            obj.rawFcn = rawFcn;
            obj.maxFEs = maxFEs;
            obj.cpFEs  = unique(checkpoints(:).');
            obj.cpBest = inf(size(obj.cpFEs));
           
            grid = unique(round(logspace(0, log10(maxFEs), nLog)));
            grid(grid < 1) = [];
            obj.logFEs  = unique([grid(:).', obj.cpFEs, maxFEs]);
            obj.logBest = inf(size(obj.logFEs));
        end

        function y = eval(obj, x)
            if isvector(x)
                y = obj.evalOne(x(:).');
            else
                m = size(x,1);
                y = zeros(m,1);
                for k = 1:m
                    y(k) = obj.evalOne(x(k,:));
                end
            end
        end

        function y = evalOne(obj, x)
            if obj.exhausted
                throw(MException('BudgetObjective:Exhausted','FEs budget exhausted'));
            end
            y = obj.rawFcn(x(:).');
            obj.count = obj.count + 1;
            if isfinite(y) && y < obj.best
                obj.best  = y;
                obj.bestX = x;
            end
            obj.record();
            if obj.count >= obj.maxFEs
                obj.exhausted = true;
                
                obj.fillRest();
                throw(MException('BudgetObjective:Exhausted','FEs budget exhausted'));
            end
        end

        function record(obj)
            while obj.cpIdx <= numel(obj.cpFEs) && obj.count >= obj.cpFEs(obj.cpIdx)
                obj.cpBest(obj.cpIdx) = obj.best;
                obj.cpIdx = obj.cpIdx + 1;
            end
            while obj.logIdx <= numel(obj.logFEs) && obj.count >= obj.logFEs(obj.logIdx)
                obj.logBest(obj.logIdx) = obj.best;
                obj.logIdx = obj.logIdx + 1;
            end
        end

        function fillRest(obj)
            obj.cpBest(isinf(obj.cpBest)) = obj.best;
            obj.logBest(isinf(obj.logBest)) = obj.best;
        end

        function finalize(obj)
            obj.fillRest();
        end
    end
end
