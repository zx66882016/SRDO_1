function Group = Initialization_GLO(Sol, N)

 
    Groupsizes = cumsum([floor(N/3), floor(N/3), N-2*floor(N/3)]);
    Group(1).Member = Sol(1:Groupsizes(1));
    Group(2).Member = Sol(Groupsizes(1)+1:Groupsizes(2));
    Group(3).Member = Sol(Groupsizes(2)+1:Groupsizes(3));
    
    for i = 1 : 3
        [~,ind] = sort([Group(i).Member.Cost]);
        Group(i).Member = Group(i).Member(ind);
        Group(i).Alpha = Group(i).Member(1);
    end
    
end

