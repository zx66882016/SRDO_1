function newX = update_female(Member, Nm, i, D, ASF)

    mate = RouletteWheelSelection([Member(1:Nm).Cost]);
    newX = Member(i).X + (randn(1, D) * .8 * ASF) .* (Member(mate).X - Member(i).X);

end

