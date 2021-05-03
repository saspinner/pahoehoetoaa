%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%compiler for plotting J


fs = 1;
[J,GG] = Jvalue(fs);
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
    mm     = [.1 .2 .33];
    
    for ll = 1:numel(mm)
        m     = mm(ll);
        mu_u  = mu_l./m;
        nn = 0;
        K_p   = 1*10e2:10e2:10e4;
        Ku_p   =  K_p ;
        U_p   = Ku_p.*du^2./mu_u;
        F     = sqrt((rho_l-rho_u)./rho_u*g*du./U_p.^2);
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
[FF,N] = meshgrid(F,n);
IJ1 = JJJ1./abs(JJJ1);
IJ2 = JJJ2./abs(JJJ2);
IJ2(IJ2<0) = 0;
IJ3 = JJJ3./abs(JJJ3);
IJ3(IJ3<0) = 0;
IJ = IJ1 + IJ2 + IJ3;
JplotRD(N,FF,IJ,fs)

