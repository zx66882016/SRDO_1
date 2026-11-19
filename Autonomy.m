function [Group, Best] = Autonomy(Group, Best, ct, mt, lb, ub, D, fobj)

    for i = 1 : 3
       for j = 1 : numel(Group(i).Member)
            levy_scale = rand .* (ub-lb) * exp(-20 * (ct/mt)^2);
            levy = 1 ./ (abs(randn(1, D)) .^ (1/1.5));
            step = levy_scale .* levy * .001;
            new.X = Group(i).Member(j).X .* step;
           [Group(i), Best] = Elitism(Group(i), Best, new, fobj, j, lb, ub);
       end
    end

end
