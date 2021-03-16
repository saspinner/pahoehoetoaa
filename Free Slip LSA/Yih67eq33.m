%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%solving system of equations for%%%%%%%%%%%
%%%%%%%Yih instability with free slip%%%%%%%%
%%%%%%%%%%%%top condition%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
syms m n c du U_p K_p mu_u ...
    Au Bu Cu Du Al Bl Cl Dl c0p Ra r alpha ...
    au al F DeltaBu DeltaCu DeltaDu DeltaBl DeltaCl DeltaDl...
    bu bl hu_1 dhu_1 d2hu_1 hl_n dhl_n
% syms c du mu_u ...
%     Au Bu Cu Du Al Bl Cl Dl c0p Ra alpha ...
%     au al DeltaBu DeltaCu DeltaDu DeltaBl DeltaCl DeltaDl...
%     bu bl hu_1 dhu_1 d2hu_1 hl_n dhl_n
fs = 1;   %free slip? 1 is yes 0 is no 



%%%%%%%%%%%%%%%%%%%%%%%%%%equation 7 and 8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fs == 1
    fs_7 = [2 1 0 0 0 0];
    b_7 = [0; 0; 0; 0; 0; (-du^2*K_p)/(2*U_p*mu_u)];
else
    fs_7 = [1 1 1 0 0 0];
    b_7 = [0; 0; 0; 1; 0; (-du^2*K_p)/(2*U_p*mu_u)];
end
%steady state solution
%       Au au bu Al   al bl
A_7  = [0   0  1  0   0 -1;
        0   1  0  0   -m 0;
        0   0  0  n^2 -n 1;
        fs_7; %<=free slip?
        1   0  0  -m   0 0;
        1   0  0   0   0 0];
    
M_7 = A_7\b_7;
Au = M_7(1); au = M_7(2); bu = M_7(3); 
Al = M_7(4); al = M_7(5); bl = M_7(6);
%%%%%%%%%%%%%%%%%%%%%%%%%%equation 33%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%short form where we assume E=1, 
%we assume E=1 because E does not appear in the eigenvalue problem
if fs ==1
    fs_33 = [0  1 3  0   0    0 0];
else
    fs_33 = [1 2 3 0 0 0 0];
end
%%%%Bu Cu Du Bl Cl   Dl E
A_33 = [0  1 0  0  -m    0 0;
     0  0 1  0   0   -m 0;
     1  1 1  0   0    0 1;
     fs_33; %<free slip condition
     0  0 0 -n   n^2 -n^3 1;
     0  0 0  1  -2*n 3*n^2 0];
%for what c' would give you null results
%    1  0 0 -1   0     0 -(al-au)/c]; 
b_33 = -A_33(:,end);

M_33 = A_33(:,1:end-1)\b_33;

Bu = M_33(1); Cu = M_33(2); Du = M_33(3);
Bl = M_33(4); Cl = M_33(5); Dl = M_33(6);
 
c0p = (al-au)/(Bu-Bl);
%%%%%%%%%%%%%%%%%%%%%%%%%%equation 39 and 40%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%short form where we assume DeltaE=0%%%%%%%%%%%%%%%%%%%%%%%
%we assume DeltaE=0 because DeltaE does not appear in the eigenvalue
%problem
%%%%%%upper
hu_1   = Au*Du/210 + au*Du/60 + (-Au*Bu+au*Cu-3*c0p*Du)/60 -(Cu*c0p + Au)/12;
dhu_1   = Au*Du/30 + au*Du/10 + (-Au*Bu+au*Cu-3*c0p*Du)/12 -(Cu*c0p + Au)/3;
d2hu_1 = Au*Du/5   + au*Du/2  + (-Au*Bu+au*Cu-3*c0p*Du)/3  -(Cu*c0p + Au);
%%%%%%lower
hl_n   = Al*Dl*(-n)^7/210 + al*Dl*(-n)^6/60 + (-Al*Bl+al*Cl-3*c0p*Dl)*(-n)^5/60 -(Cl*c0p + Al)*(-n)^4/12;
dhl_n   = Al*Dl*(-n)^6/30 + al*Dl*(-n)^5/10 + (-Al*Bl+al*Cl-3*c0p*Dl)*(-n)^4/12 -(Cl*c0p + Al)*(-n)^3/3;

if fs == 1
    fs_3840 = [0       2       6        0       0        0   1i*alpha*Ra*d2hu_1];
else
    fs_3840 = [1       2       3        0       0        0   1i*alpha*Ra*dhu_1];
end
%%%%%%%% DeltaBu DeltaCu DeltaDu DeltaBl DeltaCl DeltaDl contants
A_3840 = [ 1       1       1        0       0        0   1i*alpha*Ra*hu_1;
           fs_3840; 
           0       0       0       -n      n^2     -n^3  1i*alpha*Ra*r/m*hl_n;
           0       0       0        1      -2*n    3*n^2 1i*alpha*Ra*r/m*dhl_n;
           0       -1      0        0       m        0          0;
           0       0       -6        0      0        6*m -1i*alpha*Ra*((F^2/(c0p))-(r-1)*(c0p*Bu+au))];
%for what c' would give you null results
%    1  0 0 -1   0     0 -(al-au)delc/c^2]; 
b_3840 =  -A_3840(:,end);

M_3840 = A_3840(:,1:end-1)\b_3840;

DeltaBu = M_3840(1); DeltaCu = M_3840(2); DeltaDu = M_3840(3);
DeltaBl = M_3840(4); DeltaCl = M_3840(5); DeltaDl = M_3840(6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculating ci

ci = -(DeltaBu - DeltaBl)*c0p^2/(1i*(al-au));
J = ci/(Ra*alpha);
G = J*(au-al)*m/c0p^2;
simplify(G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 1.1:.1:70;
n = 1.25;
r = 1;
K_p = 10;
du = 1;
mu_u = 10^3;
U_p = K_p*du^2/mu_u;
F = 0;


JJ = subs(J*10^4);
plot(m, JJ,'-b')
hold on
n = 2.5;
JJ = subs(J*10^4);
plot(m, JJ,'-b')
hold on
n = 5;
JJ = subs(J*10^4);
plot(m, JJ,'-b')
hold on
n = 10;
JJ = subs(J*10^4);
plot(m, JJ,'-b')