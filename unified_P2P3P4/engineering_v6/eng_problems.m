function P = eng_problems()
P = struct([]);

% --------------------------------
q = struct();
q.tag='welded_beam'; q.name_en='Welded beam design'; q.name_cn='Welded beam design';
q.dim=4;
q.lb=[0.125 0.1 0.1 0.1]; q.ub=[5 10 10 5];
q.mask=[0 0 0 0]; q.discreteSet=[];
q.cScale=ones(1,6);
q.varNames={'h weld leg size', 'l weld seam length', 't beam height', 'b beam width'};
q.fref=1.724852;
q.xref=[0.205730 3.470489 9.036624 0.205730];
q.rawf=@weld_f; q.rawg=@weld_g;
q.note='Ragsdell-Phillips(1976)/Coello Standard Model，polar moment of inertia J=2*sqrt(2)*h*l*(l^2/12+(h+t)^2/4)';
P=addp(P,q);

% ---------------- 2. Pressure vessel design----------------
q = struct();
q.tag='pressure_vessel'; q.name_en='Pressure vessel design'; q.name_cn='Pressure vessel design';
q.dim=4;
q.lb=[0.0625 0.0625 0 0]; q.ub=[10 10 100 240];
q.mask=[0.0625 0.0625 0 0]; q.discreteSet=[];
q.cScale=[1 1 1.296e6];
q.varNames={'T_s Shell thickness','T_h head thickness','R inner diameter','L cylinder length'};
q.fref=6059.7143;
q.xref=[0.8125 0.4375 42.098446 176.636596];
q.rawf=@ves_f; q.rawg=@ves_g;
q.note='T_s and T_h are integer multiples of 0.0625 inch';
P=addp(P,q);

% ---------------- 3. Speed reducer design----------------
q = struct();
q.tag='speed_reducer'; q.name_en='Speed reducer design'; q.name_cn='Speed reducer design';
q.dim=7;
q.lb=[2.6 0.7 17 7.3 7.3 2.9 5]; q.ub=[3.6 0.8 28 8.3 8.3 3.9 5.5];
q.mask=[0 0 1 0 0 0 0]; q.discreteSet=[];
q.cScale=ones(1,11);
q.varNames={'x_1 face width','x_2 modulus','x_3 number of teeth on the pinion','x_4 Bearing spacing of shaft 1','x_5 Axle 2 bearing spacing', ...
            'Diameter of x_6 axis 1','Diameter of x_7 axis 2'};
q.fref=2994.4245;
q.xref=[3.5 0.7 17 7.3 7.715319 3.350215 5.286654];
q.rawf=@red_f; q.rawg=@red_g;
q.note='x_3 The number of teeth shall be an integer；According to the text, the optimal continuous numerical value of the constant is approximately 2994.471, while the literature reports it as 2994.4245';
P=addp(P,q);

% ---------------- 4. Three-bar truss design----------------
q = struct();
q.tag='three_bar_truss'; q.name_en='Three-bar truss design'; q.name_cn='Three-bar truss design';
q.dim=2;
q.lb=[0 0]; q.ub=[100 100];
q.mask=[0 0]; q.discreteSet=[];
q.cScale=ones(1,3);
q.varNames={'x_1=x_3 Cross-sectional area of the inclined bar','x_2 Vertical bar cross-sectional area'};
q.fref=263.8958;
q.xref=[0.788675 0.408248];
q.rawf=@truss_f; q.rawg=@truss_g;
q.note='symmetrical structure，P=2，allowable stress sigma=2';
P=addp(P,q);

% ---------------- 5. Reinforced concrete beam design----------------
q = struct();
q.tag='rc_beam'; q.name_en='Reinforced concrete beam design'; q.name_cn='Reinforced concrete beam design';
q.dim=3;
q.lb=[6.0 28 5]; q.ub=[8.4 40 10];
q.mask=[0 1 0];
As_set=[6.0 6.16 6.32 6.6 7.0 7.11 7.2 7.8 7.9 8.0 8.4];
q.discreteSet={As_set, [], []};
q.cScale=ones(1,3);
q.varNames={'A_s Area of tensile steel bar (in^2)','b Liang Kuan (in)','h Liang Gao (in)'};
q.fref=359.2080;
q.xref=[6.32 34 8.5];
q.rawf=@rcb_f; q.rawg=@rcb_g;
q.note='Amir-Hasegawa Mixed discrete problems：A_s Take standard rebar specifications，b Round to the nearest whole inch；fc=5ksi, fy=50ksi, d=0.85h, Required bending resistance 1350 kip-in';
P=addp(P,q);

% ---------------- 6. Tubular column design----------------
q = struct();
q.tag='tubular_column'; q.name_en='Tubular column design'; q.name_cn='Tubular column design';
q.dim=2;
q.lb=[2 0.2]; q.ub=[14 0.8];
q.mask=[0 0]; q.discreteSet=[];
q.cScale=ones(1,6);
q.varNames={'d Average pipe diameter','t Tube wall thickness'};
q.fref=26.5314;
q.xref=[5.451153 0.291967];
q.rawf=@col_f; q.rawg=@col_g;
q.note='P=2500 kg, sigma_y=500 kg/cm^2, E=0.85e6 kg/cm^2, L=250 cm；Optimality in continuous numerical values 26.4995';
P=addp(P,q);

end

function P=addp(P,q), q.id=numel(P)+1; if isempty(P), P=q; else, P(end+1)=q; end, end

% ============ 1. Welded beam ============
function f=weld_f(x)
h=x(1);l=x(2);t=x(3);b=x(4);
f=1.10471*h^2*l+0.04811*t*b*(14+l);
end
function g=weld_g(x)
h=x(1);l=x(2);t=x(3);b=x(4);
P=6000;L=14;E=3e7;G=1.2e7;
M=P*(L+l/2);
R=0.5*sqrt(l^2+(h+t)^2);
J=2*sqrt(2)*h*l*(l^2/12+0.25*(h+t)^2);
tau1=P/(sqrt(2)*h*l);
tau2=M*R/J;
tau=sqrt(tau1^2+tau1*tau2*l/R+tau2^2);
sig=6*P*L/(b*t^2);
Pc=4.013*E/(6*L^2)*t*b^3*(1-0.25*t/L*sqrt(E/G));
delta=4*P*L^3/(E*b*t^3);
f=weld_f(x);
g=[tau/13600-1; sig/30000-1; h-b; ...
   0.10471*h^2+0.04811*t*b*(14+l)-5; delta-0.25; P-Pc];
end

% ============ 2. pressure vessel ============
function f=ves_f(x)
Ts=round(x(1)/0.0625)*0.0625;
Th=round(x(2)/0.0625)*0.0625;
R=x(3);L=x(4);
f=0.6224*Ts*R*L+1.7781*Th*R^2+3.1661*Ts^2*L+19.84*Ts^2*R;
end
function g=ves_g(x)
Ts=round(x(1)/0.0625)*0.0625;
Th=round(x(2)/0.0625)*0.0625;
R=x(3);L=x(4);
g=[-Ts+0.0193*R; -Th+0.00954*R; -pi*R^2*L-(4/3)*pi*R^3+1296000];
end

% ============ 3. reducer ============
function f=red_f(x)
x1=x(1);x2=x(2);x3=round(x(3));x4=x(4);x5=x(5);x6=x(6);x7=x(7);
f=0.7854*x1*x2^2*(3.3333*x3^2+14.9334*x3-43.0934) ...
  -1.508*x1*(x6^2+x7^2)+7.4777*(x6^3+x7^3)+0.7854*(x4*x6^2+x5*x7^2);
end
function g=red_g(x)
x1=x(1);x2=x(2);x3=round(x(3));x4=x(4);x5=x(5);x6=x(6);x7=x(7);
g=[ 27/(x1*x2^2*x3)-1;
    397.5/(x1*x2^2*x3^2)-1;
    1.93*x4^3/(x2*x3*x6^4)-1;
    1.93*x5^3/(x2*x3*x7^4)-1;
    sqrt((745*x4/(x2*x3))^2+16.9e6)/(110*x6^3)-1;
    sqrt((745*x5/(x2*x3))^2+157.5e6)/(85*x7^3)-1;
    x2*x3/40-1;
    5*x2/x1-1;
    x1/(12*x2)-1;
    (1.5*x6+1.9)/x4-1;
    (1.1*x7+1.9)/x5-1];
end

% ============ 4. Three-bar truss ============
function f=truss_f(x)
x1=x(1);x2=x(2);
f=100*(2*sqrt(2)*x1+x2);
end
function g=truss_g(x)
x1=x(1);x2=x(2);P=2;sig=2;
den=sqrt(2)*x1^2+2*x1*x2;
g=[ (sqrt(2)*x1+x2)/den*P-sig;
    x2/den*P-sig;
    1/(x1+x2)*P-sig];
end

% ============ 5. Reinforced concrete beam ============
function f=rcb_f(x)
As=x(1);b=x(2);h=x(3);
f=29.4*As+0.6*b*h;
end
function g=rcb_g(x)
As=x(1);b=x(2);h=x(3);
fc=5; fy=50;                
d=0.85*h;                    
Mu_req=1350;                 
a=As*fy/(0.85*fc*b);
phiMn=0.9*As*fy*(d-a/2);
g=[ Mu_req/phiMn-1;
    (180+7.375*As^2/h)/(As*b)-1;
    b/h-4];
end

% ============ 6. Tubular column ============
function f=col_f(x)
d=x(1);t=x(2);
f=9.8*d*t+2*d;
end
function g=col_g(x)
d=x(1);t=x(2);
P=2500;sy=500;E=0.85e6;L=250;
g=[ P/(pi*d*t*sy)-1;
    8*P*L^2/(pi^3*E*d*t*(d^2+t^2))-1;
    2/d-1;
    d/14-1;
    0.2/t-1;
    t/0.8-1];
end
