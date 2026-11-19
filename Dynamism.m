function Group = Dynamism(Group, omega, phi)

    for i = 1 : 3
        members = [];
        for j = 1 : numel(Group(i).Member)
            if rand <= omega && numel(Group(i).Member) > phi
                DG = randi([setdiff(1:3,i)]);  

                
                PlaceIndex = numel(Group(DG).Member) + 1;
                Group(DG).Member(PlaceIndex).X = Group(i).Member(j).X;
                Group(DG).Member(PlaceIndex).Cost = Group(i).Member(j).Cost;

                members = [members j]; 
                if numel(Group(i).Member)-numel(members) < 5
                    break;
                end
            end
        end

        if numel(members) > 0
            Group(i).Member(members) = [];
        end
    end

end

