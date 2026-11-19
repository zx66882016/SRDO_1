function [Group, Best] = Elitism(Group, Best, new, fobj, i, lb, ub)

    new.X = min(max(new.X, lb), ub);
    new.Cost = fobj(new.X);

    if new.Cost < Group.Member(i).Cost
        Group.Member(i) = new;
    end
    if new.Cost < Group.Alpha.Cost
        Group.Alpha = new;
    end
    if new.Cost < Best.Cost
        Best = new;
    end
            
end

