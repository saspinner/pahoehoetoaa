%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%solving system of equations for%%%%%%%%%%%
%%%%%%%Yih instability with free slip%%%%%%%%
%%%%%%%%%%%%top condition%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
syms m n a_u a_l c du U_p K_p mu_u ...
    Au Bu Cu Du Al Bl Cl Dl c0p Ra r alpha ...
    au al F


%%%%%%%%%%%%%%%%%%%%%%%%%%equation 7 and 8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%steady state solution
A_7  = [0 0 1 0 0 -1;
      0 1 0 0 -m 0;
      0 0 0 n^2 -n 1;
      2 1 0 0 0 0;
      1 0 0 -m 0 0;
      1 0 0 0 0 0];
b_7 = [0; 0; 0; 0; 0; (-du^2*K_p)/(2*U_p*mu_u)]; %steady state solution
%%%%%%%%%%%%%%%%%%%%%%%%%%equation 33%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%short form where we assume E=1, 
%we assume E=1 because E does not appear in the eigenvalue problem
%%%%Bu Cu Du Bl Cl   Dl E
A_33 = [0  1 0  0  -m    0 0;
     0  0 1  0   0   -m 0;
     1  1 1  0   0    0 1;
     0  1 3  0   0    0 0; %<free slip condition
     0  0 0 -n   n^2 -n^3 1;
     0  0 0  1  -2*n 3*n^2 0];
%for what c' would give you null results
%    1  0 0 -1   0     0 -(a_l-a_u)/c]; 
b_33 = -A_33(:,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%equation 39 and 40%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%short form where we assume DeltaE=0%%%%%%%%%%%%%%%%%%%%%%%
%we assume DeltaE=0 because DeltaE does not appear in the eigenvalue
%problem
%%%%%%upper
hu_1   = Au*Du/210 + au*Du/60 + (-Au*Bu+au*Cu-3*c0p*Du)/60 -(Cu*c0p + Au)/12;
d2hu_1 = Au*Du/5   + au*Du/2  + (-Au*Bu+au*Cu-3*c0p*Du)/3  -(Cu*c0p + Au);
%%%%%%lower
hl_n   = Al*Dl/210 + al*Dl/60 + (-Al*Bl+al*Cl-3*c0p*Dl)/60 -(Cl*c0p + Al)/12;
dhl_n   = Al*Dl*(-n)^6/30 + al*Dl*(-n)^5/20 + (-Al*Bl+al*Cl-3*c0p*Dl)*(-n)^4/12 -(Cl*c0p + Al)*(-n)^3/3;

%%%%%%%% DeltaBu DeltaCu DeltaDu DeltaBl DeltaCl DeltaDl contants
A_3840 = [ 1       1       1        0       0        0   1i*alpha*Ra*hu_1;
           0       2       6        0       0        0   1i*alpha*Ra*d2hu_1;
           0       0       0       -n      n^2     -n^3  1i*alpha*Ra*r/m*hl_n;
           0       0       0        1      -2*n    3*n^2 1i*alpha*Ra*r/m*dhl_n;
           0       -1      0        0       m        0          0;
           0       0       6        -1      0        6*m 1i*alpha*Ra*(1/c0p*F^2-(r-1)*(c0p*Bu+au))];
%for what c' would give you null results
%    1  0 0 -1   0     0 -(a_l-a_u)/c]; 
b_3840 =  -A_3840(:,end);






A = A_3840;
b = b_3840;

inv(A(:,1:end-1))*b