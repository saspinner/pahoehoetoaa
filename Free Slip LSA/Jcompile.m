%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compiler for plotting J
clear all
close all
%fs = 1 is free slip
%fs = 2 is no slip
%fs = 0 is constant slip
fs = 2;
[J,GG,U_u,U_l] = Jvalue(fs);
if fs == 1
    n     = 5;
    rho_l = 3000;
    rho_u = 2970;
    r     = rho_l/rho_u;
    du    = 1;
    dl    = du*n;
    mu_u  = 10^3;
    g     = 9.8;
    mm     = [.95 .9 .85];
    K_p   = logspace(2,6,10);
    Ku_p   =  K_p ;
    U_p   = Ku_p.*du^2./mu_u;
    F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    
    for ll = 1:numel(mm)
        m     = mm(ll);
        mu_l  = mu_u.*m;
        nn = 0;
        
        for n = .2:.2:6
            nn = nn + 1;
            dl    = du*n;
            gamma = n.^3.*r./m.^2;
            Kl_p   =  K_p./gamma ;
            if ll==1
                JJ1(nn,:) = subs(J*10^4);
            elseif ll == 2
                JJ2(nn,:) = subs(J*10^4);
            else
                JJ3(nn,:) = subs(J*10^4);
            end
        end
        nn = 0;
        n = .2:.2:6;
        dl    = du.*n;
        if ll ==1
            JJJ1 = double(JJ1);
        elseif ll == 2
            JJJ2 = double(JJ2);
        else
            JJJ3 = double(JJ3);
        end
    end
elseif fs == 2
    n     = 5;
    rho_l = 3000;
    rho_u = 2970;
    r     = rho_l/rho_u;
    dl    = 1;
    du    = dl/n;
    mu_l  = 10^3;
    g     = 9.8;
    mm     = [.1 .2 .33];
    K_p   = logspace(2,6,10);
    Kl_p   =  K_p ;
    U_p   = Kl_p.*dl^2./mu_l;
    F     = sqrt((rho_l-rho_u)./rho_u*g*dl./U_p.^2);
    
    for ll = 1:numel(mm)
        m     = mm(ll);
        mu_u  = mu_l./m;
        nn = 0;
        
        for n = .2:.2:6
            nn = nn + 1;
            du    = dl/n;
            Ku_p   =  K_p;
            if ll==1
                JJ1(nn,:) = subs(J*10^4);
            elseif ll == 2
                JJ2(nn,:) = subs(J*10^4);
            else
                JJ3(nn,:) = subs(J*10^4);
            end
        end
        nn = 0;
        n = .2:.2:6;
        du    = dl./n;
        if ll ==1
            JJJ1 = double(JJ1);
        elseif ll == 2
            JJJ2 = double(JJ2);
        else
            JJJ3 = double(JJ3);
        end
    end
elseif fs == 0
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
[FF,N] = meshgrid(F,n);
IJ1 = JJJ1./abs(JJJ1);
IJ2 = JJJ2./abs(JJJ2);
IJ2(IJ2<0) = 0;
IJ3 = JJJ3./abs(JJJ3);
IJ3(IJ3<0) = 0;
IJ = IJ1 + IJ2 + IJ3;
JplotRD(N,FF,IJ,fs)

%plot steady state
    n     = 5;
    rho_l = 3000;
    rho_u = 2970;
    r     = rho_l/rho_u;
    du    = 1;
    dl    = du*n;
    K_p   = 0;
    mu_l  = 10^3;
    g     = 9.8;
    mm     = [.1 .2 .33];
    if fs == 1 
       K_p   = 10^5;
       Ku_p   =  K_p ;
       U_p   = Ku_p.*du^2./mu_u;
       F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    elseif fs == 2
       K_p   = 10^5;
       Ku_p  =  K_p ;
       Kl_p  = K_p;
       U_p   = Kl_p.*dl^2./mu_l;
       F     = sqrt((rho_l-rho_u)./rho_u*g*dl./U_p.^2);
    elseif fs == 0
        K_p  = 0;
        U_p  = 10;
        F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
    end

    y     = -n:.01:1;
    for ll = 1:numel(mm)
        m     = mm(ll);
        mu_u  = mu_l./m;
        gamma = n.^3.*r./m.^2;
        if fs == 1 
            Kl_p   =  K_p./gamma ;
        elseif fs == 2 || fs ==0
            Kl_p = Ku_p;
        end
        
        U_u_(ll,:) = subs(U_u);
        U_l_(ll,:) = subs(U_l);
    end
  JplotSS(U_u_,U_l_,y,fs)


