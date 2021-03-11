%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%solving system of equations for%%%%%%%%%%%
%%%%%%%Yih instability with free slip%%%%%%%%
%%%%%%%%%%%%top condition%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms m n a_u a_l c du U_p K_p mu_u

%equation 33
%long form
%%%%%%%  Eu Bu Cu Du Fu El            Bl Cl    Dl   Fl
%  A = [ 0  0   0  0  1 0             0  0    0     0; %1
%        0  0   0  0  0 0             0  0    0     1; %2
%        0  0   2  0  0 0             0  -2*m 0     0; %3
%        0  0   0  1  0 0             0  0    -m    0; %4
%        1  1   1  1  1 0             0  0    0     0; %5
%        0  0   2  6 12 0             0  0    0     0; %6
%        0  0   0  0  0 1             -n n^2  -n^3  n^4; %7
%        0  0   0  0  0 0             1 -2*n  3*n^2 -4*n^3; %8
%        1  0   0  0  0 -1             0  0     0      0; %9
%-(a_l-a_u)/c  1   0  0  0 0  -1 0     0      0] %10 

%short form
%%%%Bu Cu Du Bl Cl   Dl
A = [0  1 0  0  -m    0;
      0  0 1  0   0   -m;
     1  1 1  0   0    0;
      0  1 3  0   0    0;
     0  0 0 -n   n^2 -n^3;
     0  0 0  1  -2*n 3*n^2;
     1  0 0 -1   0     0];

%steady state solution
% A  = [0 0 1 0 0 -1;
%       0 1 0 0 -m 0;
%       0 0 0 n^2 -n 1;
%       2 1 0 0 0 0;
%       1 0 0 -m 0 0;
%       1 0 0 0 0 0];


% b = [0; 0; 0; 0; 0; (-du^2*K_p)/(2*U_p*mu_u)]; %steady state solution
%long form
%     1  2  3  4  5   6  7   8  9   10
% b = [ 0; 0; 0; 0; 0;  0; 0;  0; 0;  0];

%short form
b = [0; 0; -1; 0; -1; 0; (a_l-a_u)/c];
inv(A)*b