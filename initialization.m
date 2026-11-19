function X=initialization(N,Dim,UB,LB)

B_no= size(UB,2); 

if B_no==1
    X=rand(N,Dim).*(UB-LB)+LB;
end

if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
    end
