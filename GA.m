function [Best_score,Best_pos,curve]=GA(pop,Max_iter,lb,ub,dim,fobj)
popsize=pop;            
lenchrom=dim;              
fun = fobj;  
pc=0.6;                  
pm=0.3;                  
if(max(size(ub)) == 1)
   ub = ub.*ones(dim,1);
   lb = lb.*ones(dim,1);  
end
maxgen=Max_iter; 

bound=[lb(:),ub(:)];  

for i=1:popsize
 
    GApop(i,:)=Code(lenchrom,bound);       
    [fitness(i)]=fun(GApop(i,:));           
end

%Find the best chromosome
[bestfitness bestindex]=min(fitness);
zbest=GApop(bestindex,:);   
gbest=GApop;                
fitnessgbest=fitness;       
fitnesszbest=bestfitness;   
%% Iterative optimization
for i=1:maxgen
        GApop=Select2(GApop,fitness,popsize);
        GApop=Cross(pc,lenchrom,GApop,popsize,bound);
        GApop=Mutation(pm,lenchrom,GApop,popsize,[i maxgen],bound);
        pop=GApop;
      for j=1:popsize
        %fitness value
        [fitness(j)]=fun(pop(j,:));
        %Individual optimal update
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %Group optimal update
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
        
    end
    
    curve(i)=fitnesszbest;    
        display(['GA:  At iteration ', num2str(i), ' the fitness is ', num2str(fitnesszbest)]);
end
Best_score = fitnesszbest;
Best_pos = zbest;
end
%% Selection function
function ret=Select2(individuals,fitness,sizepop)
fitness= 1./(fitness);
sumfitness=sum(fitness);
sumf=fitness./sumfitness;
index=[];
for i=1:sizepop 
    pick=rand;
    while pick==0
        pick=rand;
    end
    for j=1:sizepop
        pick=pick-sumf(j);
        if pick<0
            index=[index j];
            break;  
        end
    end
end
individualsTemp=individuals(index,:);
fitnessTemp=fitness(index);
if(size(individualsTemp,1) == 0)
    ret=individuals;
else
    ret=individualsTemp;
end
end
%% Cross function
function ret=Mutation(pmutation,lenchrom,chrom,sizepop,pop,bound)
for i=1:sizepop 
    pick=rand;
    while pick==0
        pick=rand;
    end
    index=ceil(pick*sizepop);
    pick=rand;
    if pick>pmutation
        continue;
    end
    flag=0;
    while flag==0
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick*sum(lenchrom));  
        if pos<=0 
            pos = 1;
        end
        if pos>size(bound,1)
            pos = size(bound,1);
        end
        v=chrom(i,pos);
        v1=v-bound(pos,1);
        v2=bound(pos,2)-v;
        pick=rand; 
        if pick>0.5
            delta=v2*(1-pick^((1-pop(1)/pop(2))^2));
            chrom(i,pos)=v+delta;
        else
            delta=v1*(1-pick^((1-pop(1)/pop(2))^2));
            chrom(i,pos)=v-delta;
        end   
        flag=test(lenchrom,bound,chrom(i,:));    
    end
end
ret=chrom;
end

%% variogram
function ret=Cross(pcross,lenchrom,chrom,sizepop,bound)
for i=1:sizepop 
 
    pick=rand(1,2);
    while prod(pick)==0
        pick=rand(1,2);
    end
    index=ceil(pick.*sizepop);
    pick=rand;
    while pick==0
        pick=rand;
    end
    if pick>pcross
        continue;
    end
    flag=0;
    while flag==0
        pick=rand;
        while pick==0
            pick=rand;
        end
        pos=ceil(pick.*sum(lenchrom)); 
        pick=rand;
        v1=chrom(index(1),pos);
        v2=chrom(index(2),pos);
        chrom(index(1),pos)=pick*v2+(1-pick)*v1;
        chrom(index(2),pos)=pick*v1+(1-pick)*v2; 
        flag1=test(lenchrom,bound,chrom(index(1),:));  
        flag2=test(lenchrom,bound,chrom(index(2),:));  
        if   flag1*flag2==0
            flag=0;
        else flag=1;
        end   
    end
end
ret=chrom;
end
function flag=test(lenchrom,bound,code)


flag=1;
[n,m]=size(code);

for i=1:n
    if code(i)<bound(i,1) || code(i)>bound(i,2)
        flag=0;
    end
end
end

function ret=Code(lenchrom,bound)


flag=0;
while flag==0
    pick=rand(1,lenchrom);
    ret=bound(:,1)'+(bound(:,2)-bound(:,1))'.*pick; 
    flag=test(lenchrom,bound,ret);             
end
end