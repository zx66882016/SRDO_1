
function [fMin , bestX, Convergence_curve ] = DBO(pop, M,c,d,dim,fobj  )

P_percent = 0.2;    
pNum = round( pop *  P_percent );   


lb= c.*ones( 1,dim );    
ub= d.*ones( 1,dim );    
%Initialization
for i = 1 : pop

    x( i, : ) = lb + (ub - lb) .* rand( 1, dim );
    fit( i ) = fobj( x( i, : ) ) ;
end

pFit = fit;
pX = x;
XX=pX;
[ fMin, bestI ] = min( fit );      % fMin denotes the global optimum fitness value
bestX = x( bestI, : );             % bestX denotes the global optimum position corresponding to fMin

for t = 1 : M

    [fmax,B]=max(fit);
    worse= x(B,:);
    r2=rand(1);
    for i = 1 : pNum
        if(r2<0.8)
            r1=rand(1);
            a=rand(1,1);
            if (a>0.1)
                a=1;
            else
                a=-1;
            end
            x( i , : ) =  pX(  i , :)+0.2*abs(pX(i , : )-worse)+a*0.1*(XX( i , :)); 
        else

            aaa= randperm(180,1);
            if ( aaa==0 ||aaa==90 ||aaa==180 )
                x(  i , : ) = pX(  i , :);
            end
            theta= aaa*pi/180;

            x(  i , : ) = pX(  i , :)+tan(theta).*abs(pX(i , : )-XX( i , :));    

        end

        x(  i , : ) = Bounds( x(i , : ), lb, ub );
        fit(  i  ) = fobj( x(i , : ) );
    end
    [ fMMin, bestII ] = min( fit );      % fMin denotes the current optimum fitness value
    bestXX = x( bestII, : );             % bestXX denotes the current optimum position

    R=1-t/M;                         
    Xnew1 = bestXX.*(1-R);
    Xnew2 =bestXX.*(1+R);                   
    Xnew1= Bounds( Xnew1, lb, ub );
    Xnew2 = Bounds( Xnew2, lb, ub );
    Xnew11 = bestX.*(1-R);
    Xnew22 =bestX.*(1+R);                    
    Xnew11= Bounds( Xnew11, lb, ub );
    Xnew22 = Bounds( Xnew22, lb, ub );
    
    for i = ( pNum + 1 ) :round(pop*0.45)                 
        x( i, : )=bestXX+((rand(1,dim)).*(pX( i , : )-Xnew1)+(rand(1,dim)).*(pX( i , : )-Xnew2));
        x(i, : ) = Bounds( x(i, : ), Xnew1, Xnew2 );
        fit(i ) = fobj(  x(i,:) ) ;
    end

    for i = round(pop*0.45) +1: round(pop*0.75)                 


        x( i, : )=pX( i , : )+((randn(1)).*(pX( i , : )-Xnew11)+((rand(1,dim)).*(pX( i , : )-Xnew22)));
        x(i, : ) = Bounds( x(i, : ),lb, ub);
        fit(i ) = fobj(  x(i,:) ) ;

    end

    for j = round(pop*0.75) +1 : pop                
        x( j,: )= bestX/5+randn(1,dim).*((abs(( pX(j,:  )-bestXX)))+(abs(( pX(j,:  )-bestX))))./2;
        x(j, : ) = Bounds( x(j, : ), lb, ub );
        fit(j ) = fobj(  x(j,:) ) ;
    end
    
    XX=pX;
    for i = 1 : pop
        if ( fit( i ) < pFit( i ) )
            pFit( i ) = fit( i );
            pX( i, : ) = x( i, : );
        end

        if( pFit( i ) < fMin )
        
            fMin= pFit( i );
            bestX = pX( i, : );
       
        end
    end

    Convergence_curve(t)=fMin;




end


function s = Bounds( s, Lb, Ub)

temp = s;
I = temp < Lb;
temp(I) = Lb(I);


J = temp > Ub;
temp(J) = Ub(J);
s = temp;
function S = Boundss( SS, LLb, UUb)
temp = SS;
I = temp < LLb;
temp(I) = LLb(I);

J = temp > UUb;
temp(J) = UUb(J);

S = temp;
