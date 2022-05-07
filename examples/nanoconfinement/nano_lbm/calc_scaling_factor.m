%clear all


lbmfactor  = [];
kfactor    = [];
%for pressure=[ 1,2,5,10,20 ]
for pressure=[ 20 ]
    
    dx=0.5e-9;
    %dx=5.7e-6;
    sim_size = 256;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p0_P=pressure*1.0e6;        %pressure, [Pa]
    T0_P=400.0;                 %temperature, [k]
    
    H_L=sim_size-1;             %characterisitic length in lattice, [lu]   %changefordiffsize:)
    H_P=dx*H_L;                 %characteristic length in physical unit, [m]
    %=============================================
    pressgra_P=0.1E6;           %pressure gradient, [Pa/m]
    
    %%------------ physical unit and lattice unit ----------------------------
    Rig=8.314;                  %ideal gas constant, [J/(mol·K)]
    d_CH4=0.375e-9;             %diameter of methane molecule
    M_CH4=16.04e-3;             %molar weigth, [kg/mol]
    NA=6.02e23;                 %Avogadro’s number,[1/mol]
    mass_CH4=M_CH4/NA;          %mass of methane molecule
    w=0.01142;                  %accentic factor
    
    %------------ critical parameters for methane --------------------
    Tcr_P=190.564;              %critical temperature in physical unit,NIST, [K]
    rhocr_P=162.66;             %critical density in physical unit,NIST,     [kg/m^3]
    pcr_P=4.5992E6;             %critical pressure in physical unit,NIST,    [Pa]
    
    %===============================================================
    
    % function to obtain density, based on P-R EOS
    [Z,fhi,rho0_P] = PengRobinson(T0_P,p0_P,Tcr_P,pcr_P,w,M_CH4,0); 
    Num=rho0_P/mass_CH4;              % molecular number density of methane
    
    %dynamic viscousity from Lee et al.1966, [Pa*s]
    K_Lee=(9.379+0.01607*M_CH4)*(T0_P^1.5)/(209.2+19.26*M_CH4+T0_P);
    X_Lee=3.448+986.4/T0_P+0.01009*M_CH4;
    Y_Lee=2.447-0.2224*X_Lee;
    mu_P=1.0d-7*K_Lee*exp(X_Lee*((rho0_P/1000.0)^Y_Lee));
    
    %kinetic viscosity
    nu_P=mu_P/rho0_P;                   % Kinematic viscosity, [m^2/s]
    
    rho0_L=1.0;                         % free gas density in lattice, [lu]
    rho0_PL=rho0_P/rho0_L;              % density scale
    
    H_PL=H_P/H_L;                         % length scale, [m/lu]
    MFP_P=1.0/(sqrt(2.0)*Num*pi*d_CH4^2); % mean free path in physical unit, [m]
    
    MFP_L=MFP_P/H_PL;                   % mean free path in lattice unit, [lu]
    Kn=MFP_P/H_P;                       % Knudsen number
    
    nu_L=sqrt(2.0/pi/3.0)*MFP_L;        % Kinetic viscosity in lattice unit,[lu]
    nu_PL=nu_P/nu_L;                    % Kinetic viscosity scale, [(m^2/m)/lu]
    
    usurf_PL=nu_PL/H_PL;                % velocity scale, [(m/s)/lu]
    usurf_P=0.0;
    usurf_L=usurf_P/usurf_PL;           % velocity in lattice unit, [lu]
    
    accelz_P=pressgra_P/rho0_P;         % accelation in z direction, [m/s^2]
    accelz_PL=usurf_PL^2/H_PL;          % accelation scale, [(m/s^2)/lu]
    accelz_L=accelz_P/accelz_PL;        % accelation in lattice unit, [lu]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    lbmfactor(end+1) = 1*usurf_PL/(1e-8/accelz_L);
    kfactor(end+1)   = 1*mu_P/pressgra_P;
    
end