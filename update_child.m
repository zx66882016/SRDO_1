function newX = update_child(Member, Nm, Nf, D, ASF, lb, ub)

    mother = randi([Nm+1 max(Nm+Nf,2)]);
    father = randi([1 Nm]);
    fi = randn(1,D) .* ASF;
    mi = 1-fi;
    newX = (fi .* Member(father).X) + (mi .* Member(mother).X) + randn(1, D) .* 0.05*ASF .* (ub-lb);  % Eq. (12)

end

