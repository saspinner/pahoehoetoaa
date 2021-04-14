%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%solving system of equations for%%%%%%%%%%%
%%%%%%%Yih instability with free slip%%%%%%%%
%%%%%%%%%%%%top condition%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
syms m n c du U_p K_p mu_u Kl_p Ku_p dl mu_l...
    Au Bu Cu Du Al Bl Cl Dl c0p Ra r alpha ...
    au al F DeltaBu DeltaCu DeltaDu DeltaBl DeltaCl DeltaDl...
    bu bl hu_1 dhu_1 d2hu_1 hl_n dhl_n

fs = 0;   %free slip? 1 is yes 0 is no 



%%%%%%%%%%%%%%%%%%%%%%%%%%equation 7 and 8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Steady State Solution
if fs == 1
    fs_7 = [2 1 0 0 0 0];
    b_7 = [0; 0; 0; 0; (-du^2*Ku_p)/(2*U_p*mu_u); (-du^2*Kl_p)/(2*U_p*mu_l)];
    A_7  = [0   0  1  0   0 -1; 
        0   1  0  0   -m 0;
        0   0  0  n^2 -n 1;
        fs_7; %<=free slip?
        1   0  0  0   0 0;
        0   0  0   1   0 0];
else
    fs_7 = [1 1 1 0 0 0];
    b_7 = [0; 0; 0; 1; 0; (-du^2*K_p)/(2*U_p*mu_u)];
    %       Au au bu Al   al bl
    A_7  = [0   0  1  0   0 -1; 
        0   1  0  0   -m 0;
        0   0  0  n^2 -n 1;
        fs_7; %<=free slip?
        1   0  0  -m   0 0;
        1   0  0   0   0 0];
end

    
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
%reproduce Yih?

% m     = 1.1:1:70;
% n     = 1.25;
% r     = 1;
% K_p   = 0;
% du    = 1;
% mu_u  = 10^3;
% rho_l = 3000;
% rho_u = 3000;%2900;
% r     = rho_l/rho_u;
% U_p   = 1;%K_p*du^2/mu_u;
% g     = 9.8;
% F     = sqrt(rho_l-rho_u)/rho_u*g*du/U_p^2;

%or
if fs == 1
    n     = 5;
    rho_l = 3000;
    rho_u = 2970;
    r     = rho_l/rho_u;
    du    = 1;
    dl    = du*n;
    K_p   = 0;
    mu_l  = 10^3;
    U_p   = .2:.2:20;
    g     = 9.8;
    
    %m1
    m     = .1;
    mu_u  = mu_l./m;
    nn = 0;
    K_p   = 100:10:10000;
    Ku_p   = K_p; 
    U_p   = Ku_p.*du^2./mu_u;
    F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        gamma = n.^3.*r./m.^2;
        Kl_p   = Ku_p./gamma;
        JJ1(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ1 = double(JJ1);
    %m2
    m     = .2;
    mu_u  = mu_l./m;
    nn = 0;
    K_p   = 100:10:10000;
    Ku_p   = K_p; 
    U_p   = Ku_p.*du^2./mu_u;
    F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    n  = 0;
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        gamma = n.^3.*r./m.^2;
        Kl_p   = Ku_p./gamma;
        JJ2(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ2 = double(JJ2);
        %m3
    m     = .33;
    mu_u  = mu_l./m;
    nn = 0;
    K_p   = 100:10:10000;
    Ku_p   = K_p; 
    U_p   = Ku_p.*du^2./mu_u;
    F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    n  = 0;
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        gamma = n.^3.*r./m.^2;
        Kl_p   = Ku_p./gamma;
        JJ3(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ3 = double(JJ3);
    
    
else

    n     = 5;
    rho_l = 3000;
    rho_u = 2970;
    r     = rho_l/rho_u;
    du    = 1;
    dl    = du*n;
    K_p   = 0;
    mu_l  = 10^3;
    U_p   = .2:.2:20;
    g     = 9.8;
    F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    %m1
    m     = .1;
    mu_u  = mu_l./m;
    nn = 0;
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        JJ1(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ1 = double(JJ1);
    %m2
    m     = .2;
    mu_u  = mu_l./m;
    nn = 0;
    n  = 0;
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        JJ2(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ2 = double(JJ2);
        %m3
    m     = .33;
    mu_u  = mu_l./m;
    nn = 0;
    n  = 0;
    for n = .2:.2:6
        nn = nn + 1;
        dl    = du*n;
        JJ3(nn,:) = subs(J*10^4);
    end
    nn = 0;
    n = .2:.2:6;
    dl    = du.*n;
    JJJ3 = double(JJ3);

end




figure

[M,NN] = meshgrid(m,n);
[FF,N] = meshgrid(F,n);



% set printing options
if fs == 1
str       = {'Regime_FreeSurface'};
else
    str   = {'Regime_NoSlip'};
end
figname   = char(strcat(str(1)));
format    = '-dpng';
resl      = '-r200';
rend      = '-opengl';
printfig  = true;


% prepare formating options
HA = {'HorizontalAlignment','left','center','right'};
VA = {'VerticalAlignment','bottom','middle','top'};
UN = {'Units','Normalized','Inches'};
TX = {'Interpreter','Latex'};
TL = {'TickLabelInterpreter','Latex'};
LW = {'LineWidth',1,1.25,1.5,2};
FS = {'FontSize',10,15,18,21,24};
MS = {'MarkerSize',6,8,12};
LS = {'LineStyle','-','--','-.',':'};
% LC = {'Color',color};

% prepare axes/borders dimensions
axh = 3*2;
axw = axh;
ahs = 0.15;
avs = 0.15;
axb = 0.7;
axt = 0.2;
axl = 0.8;
axr = 0.5;
cbh = axh; cbw = 0.2;
fh = axb + 1*axh + 1*avs +           axt;
fw = axl + 1*axw + 1*ahs + 1.5*cbw + axr;

% initialize figure and axes
f = figure;
set(f,'Units','Inches','Position',[0.7 12 fw fh]);
set(f,'PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off', 'MenuBar','none');
set(f,'Resize','off','Toolbar','none');
ax(1) = axes('Units','Inches','position',[axl         axb         axw axh]);

IJ1 = JJJ1./abs(JJJ1);
IJ2 = JJJ2./abs(JJJ2);
IJ2(IJ2<0) = 0;
IJ3 = JJJ3./abs(JJJ3);
IJ3(IJ3<0) = 0;
IJ = IJ1 + IJ2 + IJ3;
axes(ax(1)); hold on;
contourf(N,FF,IJ,'edgecolor','none');
colormap([ 33,102,172;253,219,199;239,138,98;178,24,43]./255)
hold on


set(gca,'YScale', 'log')
xlabel('Height ration, n',TX{:},FS{[1,4]},UN{[1,3]},'Position',[axw/2+ahs/2,-axb/2]);
ylabel('Froude number, F',TX{:},FS{[1,4]},UN{[1,3]},'Position',[-axl/2,axh/2]);

print(f,format,resl,rend,figname,'-loose');