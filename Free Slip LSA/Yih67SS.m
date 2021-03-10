%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Steady state %%%%%%%%%%%%%%%%%
%%%%%%% Yih instability with free slip%%%%%%%
%%%%%%%%%%%% top condition %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu_l  = 10^3; 
mu_u  = 2*10^3;
m     = mu_l/mu_u;
rho_l = 3000;
rho_u = 2900;
r     = rho_l/rho_u;
du    = 1;
dl    = 2;
g     = 9.8;
n     = dl/du;

 U_p   = mu_u/rho_u/du;
K_p    = mu*U_p/du^2; 
% U_p   = K_p*du^2/mu_u;

y = -n:.01:1;

A_l =-du^2*K_p/(2*U_p*mu_l);
A_u = m*A_l;
%free slip
a_l = -2*A_l; %du^2*K_p/U_p*mu_l;
a_u = m*a_l;
b = -A_l*(n^2+2*n);

%no slip
% a_l = (1+A_l*(n^2-m))/(m+n);
% a_u = m*a_l;
% b   = (1-A_u*(1+n))*n/(m+n); 

U_u = A_u.*y.^2 +a_u.*y + b;
U_l = A_l.*y.^2 + a_l.*y + b;

U = zeros(size(y));
for i = 1:numel(y)
    if y(i)<0
        U(i) = U_l(i);
    else
        U(i) = U_u(i);
    end
end

figure
plot(U,y,'b')