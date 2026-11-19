function [Best_score,Best_pos, Convegence] = GLO(N, mt, lb, ub, D, fobj)
    
     for i = 1 : N
            Sol(i).X = unifrnd(lb,ub);   
            Sol(i).Cost = fobj(Sol(i).X);   
            X(i,:) = Sol(i).X;
            F(i) = Sol(i).Cost;
    end
    [~,index] = min([Sol.Cost]);
    Best = Sol(index);



    %% Setting parameters
    FP_om = 0.4;   % Female percentage in one-male group
    FP_mm = 0.3;   % Female percentage in multi-male group
    MP_mm = 0.3;   % Male   percentage in multi-male group
    omega = 0.05;  % Migration percentage
    phi = 5;       % Minimum member count

    %% Initialization
    Group = Initialization_GLO(Sol, N);

    %% Optimization
    for ct = 1 : mt

        ASF = 0.1 + (1 - ct/mt)^2 * .9;   % Adaptive Scaling factor
        
        %% Group 1: One-male (including women and children)
        S_om = numel(Group(1).Member); % Size of all male group
        Nf = floor(S_om*FP_om);        % Eq. (11) number of female langurs
        
        for i = 1 : S_om
            if i == 1
                % The first langur is alpha
                new.X = update_alpha(Group(1).Member(i).X, Best.X, ASF, ct, mt, lb, ub, D);
            elseif i <= Nf
                % The next Nf langures are females
                new.X = update_female(Group(1).Member, 1, i, D, ASF); 
            else
                % The ramaining langures are children
                new.X = update_child(Group(1).Member, 1, Nf, D, ASF, lb, ub);
            end
            
            % Evaluation and Elitism
            [Group(1), Best] = Elitism(Group(1), Best, new, fobj, i, lb, ub);
        end
        
        % Maintaning the social structure (reordering the langurs)
        [~,sind] = sort([Group(1).Member.Cost]);
        Group(1).Member = Group(1).Member(sind);

        %% Group 2: Multi-male (women and children)
        S_mm = numel(Group(2).Member);    % Size of multi-male group
        Nm = max(floor(S_mm*MP_mm),1);    % Number of males Eq. (16)
        Nf = max(floor(S_mm*FP_mm),1);    % Number of females Eq. (17)
        
        for i = 1 : S_mm
            if i == 1
                % The first langur is alpha
                new.X = update_alpha(Group(2).Member(i).X, Best.X, ASF, ct, mt, lb, ub, D);
            elseif i <= Nm
                % The next Nm langurs are males
                new.X = update_male(Group(2).Member(i).X, Best, ct, mt, lb, ub);
            elseif i <= Nm+Nf
                % The next Mf langurs are females
                new.X = update_female(Group(2).Member, Nm, i, D, ASF);
            else
                % The remaining langures are children
                new.X = update_child(Group(2).Member, Nm, Nf, D, ASF, lb, ub);
            end

            % Evaluation and Elitism
            [Group(2), Best] = Elitism(Group(2), Best, new, fobj, i, lb, ub);
        end
        
        % Maintaning the social structure (reordering the langurs)
        [~,sind] = sort([Group(2).Member.Cost]);
        Group(2).Member = Group(2).Member(sind);
        
        %% Group 3: All-male
        S_am = numel(Group(3).Member);   % Size of all-male group
        
        for i = 1 : S_am
            if i==1
                % The first langur is alpha
                new.X = update_alpha(Group(3).Member(i).X, Best.X, ASF, ct, mt, lb, ub, D);
            else
                % The remaining langurs are male
                new.X = update_male(Group(3).Member(i).X, Best, ct, mt, lb, ub);
            end
            
            % Evaluation and Elitism
            [Group(3), Best] = Elitism(Group(3), Best, new, fobj, i, lb, ub);
        end
        
        % Maintaning the social structure (reordering the langurs)
        [~,sind] = sort([Group(3).Member.Cost]);
        Group(3).Member = Group(3).Member(sind);

        %% Dynamism (Migration)   
        Group = Dynamism(Group, omega, phi);
        
        %% Autonomy (roaming around)
        [Group, Best] = Autonomy(Group, Best, ct, mt, lb, ub, D, fobj);
        
        %% Reporting
        Convegence(ct) = Best.Cost;  
    end
Best_score = Best.Cost;

Best_pos=Best.X;



end

