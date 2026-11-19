function newX = update_male(X, Best, ct, mt, lb, ub)
    
    C2 = 1 + randn;
    inertia = max(.4, .9 * (1 - ct / mt));
    newX = X + inertia * (X - Best.X) + C2 * (randn*(ub-lb));

end
