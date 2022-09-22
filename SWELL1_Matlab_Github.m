% Included: cytosol, F-actin, G-actin, and ions
% Steady-state
% Unknowns:
% X = [pc, vc, vn, theta_n, theta_c, cNa, cK, cCl, pH, cA, cBuf, phi, v0]
% Velocities and equations are written in the frame of the moving cell
% Parameters are in units: nm, s, Pa, mM, and mV

% cn =  [cNa,  cK,  cCl,  cH,  cHCO3,  cA,  cBuf, cHBuf];
%        (1)   (2)  (3)   (4)  (5)     (6)  (7)   (8)
% Zn =  [ 1,    1,  -1,    1,   -1,    -1,  -1,   0];
% cn0 = [cNa0, cK0, cCl0, cH0, cHCO30, cA0, cG0];
%        (1)   (2)  (3)   (4)  (5)     (6)  (7)
% Zn0 = [ 1,    1,  -1,    1,   -1,     -1,  0];

% Matlab Program for: 
% Polarized NHE1 and SWELL1 Regulate Migration Direction, Efficiency and Metastasis
% Qestions? Contact: Yizeng Li at liyizeng52@hotmail.com

clear
clc

%% Parameters are in units: nm, s, Pa, mM, and mV

L = 50.d3;              % (nm) cell length
b = 3d3;                % (nm) cell width
w = 10d3;               % (nm) cell depth
h = 500;                % (nm) membrane/cortex thickness
S = b*w;
V = L*S;
R = 8.31451;        % (J/mol K) Ideal gas constant
T = 310;            % (K) absolute temperature
F = 96485.3*(1d-3);        % (C/mol) Faraday's constant, 1d-3 comes from the rescaling of Vm

% Chemicals
N_sp   = 8;   % number of ion species, same length of cn, Zn, and Dn
n_Na   = 1;
n_K    = 2;
n_Cl   = 3;
n_H    = 4;
n_HCO3 = 5;
n_A    = 6;
n_Buf  = 7;
n_HBuf = 8;

pKc = 6.1;      % pK for bicarbonate-carbonic acid pair
pKB = 7.5;      % pK for intracellular buffer
kH = 29.76/1d3;     % (atm/mM) Henory constant at 25 degree C
PCO2 = 0.05;    % (atm) Partial pressure of CO2 concentration

NA = 9.5d-14*(1d27);  % (aMol) the 1d27 scaling comes from the scaling of the cell
NBufNHBuf = 9.5d-14*(1d27);  % (aMol) the 1d27 scaling comes from the scaling of the cell
Zn0 = [1,1,-1,1,-1,-1,0];
Zn =  [1,1,-1,1,-1,-1,-1,0];
Dn =  [1,1, 1,1, 1, 1, 1,1]*1.d8;     % (nm^2/s) diffusion constant for each species

phi0f = 0;
phi0b = 0;

cNa0f = 145;
cNa0b = 145;
cK0f  = 9;
cK0b  = 9;
cCl0f = 105;
cCl0b = 105;
cHCO30f = 35;
cHCO30b = 35;
cG0f  = 25;      % (mM) external glucose concentration
cG0b  = 25;      % (mM) external glucose concentration
cH0f = 10^(-pKc)*1d3*(PCO2/kH/cHCO30f); % in mM, not in M
cH0b = 10^(-pKc)*1d3*(PCO2/kH/cHCO30b); % in mM, not in M
pH0f = -log10(cH0f/1d3);
pH0b = -log10(cH0b/1d3);

cA0f = -1/Zn0(n_A)*(Zn0(n_Na)*cNa0f + Zn0(n_K)*cK0f + Zn0(n_Cl)*cCl0f ...
    + Zn0(n_H)*cH0f + Zn0(n_HCO3)*cHCO30f);
cA0b = -1/Zn0(n_A)*(Zn0(n_Na)*cNa0b + Zn0(n_K)*cK0b + Zn0(n_Cl)*cCl0b ...
    + Zn0(n_H)*cH0b + Zn0(n_HCO3)*cHCO30b);

cn0f = [cNa0f,cK0f,cCl0f,cH0f,cHCO30f,cA0f,cG0f];
cn0b = [cNa0b,cK0b,cCl0b,cH0b,cHCO30b,cA0b,cG0b];
c0f = sum(cn0f);
c0b = sum(cn0b);

NHEratio0 = 2.3;              % NHE polarization ratio, front to back (> 1)
SWELLratio0 = 0.6;             % SWELL polariztion ratio, front to back (< 1)

PNaf = 1d-1;
PNab = 1d-1;
PKf  = 2d1;    % (mol/J)(mM)(nm/s) = mol^2/J/micron^2/s
PKb  = 2d1;
PClf = 7d1;
PClb = PClf/SWELLratio0;
aNKEf   = 10.d6/R/T;
aNKEb   = 10.d6/R/T;
aNHE0f  = 5d2;
aNHE0b  = aNHE0f/NHEratio0;
aAE20f  = 1d2;
aAE20b  = 2.8d2;

aNaKNa  = 1d-1;
aNaKK   = 1d-2;
beta1 = 1/h;     % in GT
beta2 = h;    % in GT
beta3 = 30*(1d-3);      % in GVNKE
beta4 = -0.150*(1d3);   % in GVNKE
beta5 = 15;             % in GNHE
beta6 = 7.2;            % in GNHE
beta7 = 10;             % in GAE2
beta8 = 7.1;            % in GAE2
beta9 = 100*(1d-3);     % in GVNa
beta10 = -0.04*(1d3);   % in GVNa


% Mechanical
p0f = 0*1d5;            % (Pa) external pressure at the front
p0b = 0*1d5;            % (Pa) external pressure at the back
fextf = 0d2;            % (Pa) external force per unit area at the front of the cell
fextb = 0d2;            % (Pa) external force per unit area at the back of the cell

Thetac = 0.1;           % (mM) reference value of G-actin
Thetan = 0.2;           % (mM) reference value of F-actin
Theta  = 0.3;           % (mM) reference total G- and F-actin, \int_0^L (thetan + thetac)dx/L
thetacc  = 0.2d-3;      % (mM) Critical value for actin polymerization
Dtc =  1.d7;     % (nm^2/s) diffusion constant for G-actin

Jactinf0 = 1.0;           % (nm mM/s) Jactinf = Jactinf0*thetac^f/(thetacc + thetac^f)
Jactinf0_fixed = 2.2;

ksigman = 1d3;          % (Pa /mM) Coefficient of passive actin network pressure
gamma0 = 5d-1;          % (1/s) constant rate of actin depolymerization

alphaf = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
alphab = 1.d-1;     % (nm/Pa/s) coefficient of water permeation
dg     = 3.0d-2;         % (Pa s/nm) drag coefficient from the displaced water
kad = 5.0d-2;           % (Pa s/nm) adhesive force, Fad^b = kad*v0

eta    = 1d-9;           % (Pa s/nm^2/mM)
etast0 = 4d-5;           % (Pa s/nm^2/mM) drag coefficient between the actin and the substrate


%%
Iter = 21;
N  = 81;
dx = L/(N-1);
x  = linspace(0,L,N);

% number of variables of the model
n_var = 13; % pc, v_c, v_n, theta_n, theta_c, cNa, cK, cCl, pH, cA, cBuf, Vm, v0
var_legh = [N,1,N,N,N,N,N,N,N,N,N,N,1]; % length of eavh variable
var_strt = zeros(1,n_var); % start position of each variable
var_strt(1) = 1;
for in = 2:n_var
    var_strt(in) = var_strt(in-1) + var_legh(in-1);
end
N_var = var_strt(n_var) + var_legh(n_var) - 1; % total size of the matrix

s_pc  = var_strt(1);
s_vc  = var_strt(2);
s_vn  = var_strt(3);
s_tn  = var_strt(4);
s_tc  = var_strt(5);
s_Na  = var_strt(6);
s_K   = var_strt(7);
s_Cl  = var_strt(8);
s_pH  = var_strt(9);
s_A   = var_strt(10);
s_Buf = var_strt(11);
s_phi = var_strt(12);
s_v0  = var_strt(13);

NV = 3;
V0 = zeros(NV,1);
for ih = NV:-1:1

    if ih == 1
        NHEratio = NHEratio0;
        aNHE0b  = aNHE0f/NHEratio;
        SWELLratio = SWELLratio0;
        PClb = PClf/SWELLratio;
        Jactinf0 = Jactinf0_fixed;
    elseif ih == 2
        NHEratio = NHEratio0;
        aNHE0b  = aNHE0f/NHEratio;
        SWELLratio = 1;
        PClb = PClf/SWELLratio;
        Jactinf0 = Jactinf0_fixed;
    elseif ih == 3
        NHEratio = NHEratio0;
        aNHE0b  = aNHE0f/NHEratio;
        SWELLratio = 1.6;
        PClb = PClf/SWELLratio;
        Jactinf0 = Jactinf0_fixed;
    end

    dgf = dg/2;
    dgb = dg/2;
    etast = etast0*ones(N,1);

    % initial guess
    if ih == NV || Solution == 0
        cNa = 6*ones(N,1);
        cK  = 150.83*ones(N,1);
        cCl = 20*ones(N,1);
        pH  = 7.36*ones(N,1);
        phi = -65*ones(N,1);
        cA  = NA/V*ones(N,1);
        cBuf = NBufNHBuf/V./(1+10.^(pKB-pH));

        cH    = 10^3*10.^(-pH);
        cHCO3 = PCO2/kH*10.^(pH-pKc);
        cHBuf = cBuf.*10.^(pKB-pH);

        cin = cNa + cK + cCl + cH + cHCO3 + cA + cBuf + cHBuf;

        pc = (mean(cin)-(c0b+c0f)/2)*R*T*ones(N,1);
        thetan = linspace(Thetan*1,Thetan*1,N)';
        thetac = linspace(Thetac*1,Thetac*1,N)';
        Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
        DJactinfDtcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;
        Jwaterf = -alphaf*(pc(N)-p0f-R*T*(cin(N)-c0f));
        v0 = (etast(1)*Jactinf+dg/L*Jwaterf)/(Thetan*etast(1)+kad/L+dg/L);
        vc = v0 - Jwaterf;
        vn = linspace(v0,v0-Jactinf,N)';
    end

    X = [pc; vc; vn; thetan; thetac; cNa; cK; cCl; pH; cA; cBuf; phi; v0];

    iter = 0;
    ITER = true;
    while ITER
        iter = iter + 1;

        DF = zeros(N_var,N_var);
        Fn = zeros(N_var,1);

        %% Derived and Derivatives Terms
        % Mechancial
        sigma_n = ksigman*thetan;     % N by 1 vector
        sigma   = sigma_n;
        dsigmadtn =  ksigman;

        gamma = gamma0*logspace(-0,0,N)';
        dgammadtn = zeros(N,1);

        tauf = b/2*(sigma(N) + pc(N) - p0f - dgf*(vc+v0) + fextf);
        taub = b/2*(sigma(1) + pc(1) - p0b + dgb*(vc+v0) + fextb);

        Jactinf = Jactinf0*thetac(N)/(thetacc+thetac(N));
        DJactinfDtcN = Jactinf0*thetacc/(thetacc+thetac(N))^2;

        % Chemical
        cH    = 10^3*10.^(-pH);
        cHCO3 = PCO2/kH*10.^(pH-pKc);
        cHBuf = cBuf.*10.^(pKB-pH);

        Gnf = [PNaf, PKf, PClf, 0,0,0,0,0];
        Gnb = [PNab, PKb, PClb, 0,0,0,0,0];
        cn  = [cNa, cK, cCl, cH, cHCO3, cA, cBuf,cHBuf];
        cin = cNa + cK + cCl + cH + cHCO3 + cA + cBuf + cHBuf;

        dcHdpH = -log(10)*cH;
        dcHCO3dpH = log(10)*cHCO3;
        dcHBufdpH = -log(10)*cHBuf;
        dcdcBuf   = 1 + 10.^(pKB-pH);
        dcdpH  = dcHdpH + dcHCO3dpH + dcHBufdpH;

        % Passive channels. Tm is a function of pc, vc, thn, and v0
        Tmf = 1/(1+exp(-beta1*(tauf-beta2)));
        Tmb = 1/(1+exp(-beta1*(taub-beta2)));
        jpf = Gnf(1:3)*Tmf.*(R*T*log(cn0f(1:3)./cn(N,1:3))-Zn(1:3)*F*(phi(N)-phi0f));
        jpb = Gnb(1:3)*Tmb.*(R*T*log(cn0b(1:3)./cn(1,1:3))-Zn(1:3)*F*(phi(1)-phi0b));
        djpfdcnf  = -Gnf(1:3)*Tmf*R*T./cn(N,1:3);
        djpbdcnb  = -Gnb(1:3)*Tmb*R*T./cn(1,1:3);
        djpfdphif = -Gnf(1:3)*Tmf*F.*Zn(1:3);
        djpbdphib = -Gnb(1:3)*Tmb*F.*Zn(1:3);
        dTmfdtauf = beta1*exp(-beta1*(tauf-beta2))*Tmf^2;
        dTmfdtnf  = b/2*dsigmadtn*dTmfdtauf;
        dTmfdpcf  = b/2*dTmfdtauf;
        dTmfdvc   = -b/2*dgf*dTmfdtauf;
        dTmfdv0   = -b/2*dgf*dTmfdtauf;
        dTmbdtaub = beta1*exp(-beta1*(taub-beta2))*Tmb^2;
        dTmbdtnb  = b/2*dsigmadtn*dTmbdtaub;
        dTmbdpcb  = b/2*dTmbdtaub;
        dTmbdvc   = b/2*dgb*dTmbdtaub;
        dTmbdv0   = b/2*dgb*dTmbdtaub;

        % Na-K pump
        Gvf = 2/(1+exp(-beta3*(phi(N)-phi0f -beta4))) - 1;
        Gvb = 2/(1+exp(-beta3*(phi(1)-phi0b -beta4))) - 1;
        dGvfdphif = 2*beta3*exp(-beta3*(phi(N)-phi0f -beta4))/(1+exp(-beta3*(phi(N)-phi0f -beta4)))^2;
        dGvbdphib = 2*beta3*exp(-beta3*(phi(1)-phi0b -beta4))/(1+exp(-beta3*(phi(1)-phi0b -beta4)))^2;
        jNKEf  = -aNKEf*R*T*Gvf/(1 + aNaKNa*cNa0f/cNa(N))^3 /(1 + aNaKK*cK(N)/cK0f)^2;
        jNKEb  = -aNKEb*R*T*Gvb/(1 + aNaKNa*cNa0b/cNa(1))^3 /(1 + aNaKK*cK(1)/cK0b)^2;
        djNKEfdphif = jNKEf/Gvf*dGvfdphif;
        djNKEbdphib = jNKEb/Gvb*dGvbdphib;
        djNKEfdcNaf = 3*aNKEf*R*T*Gvf/(1 + aNaKNa*cNa0f/cNa(N))^4 ...
            /(1 + aNaKK*cK(N)/cK0f)^2 *(-aNaKNa*cNa0f/cNa(N)^2);
        djNKEbdcNab = 3*aNKEb*R*T*Gvb/(1 + aNaKNa*cNa0b/cNa(1))^4 ...
            /(1 + aNaKK*cK(1)/cK0b)^2 *(-aNaKNa*cNa0b/cNa(1)^2);
        djNKEfdcKf  = 2*aNKEf*R*T*Gvf/(1 + aNaKNa*cNa0f/cNa(N))^3 ...
            /(1 + aNaKK*cK(N)/cK0f)^3 *(aNaKK/cK0f);
        djNKEbdcKb  = 2*aNKEb*R*T*Gvb/(1 + aNaKNa*cNa0b/cNa(1))^3 ...
            /(1 + aNaKK*cK(1)/cK0b)^3 *(aNaKK/cK0b);

        % Na-H Exchanger
        GNHEf = 1/(1+exp(beta5*(pH(N)-beta6)));
        GNHEb = 1/(1+exp(beta5*(pH(1)-beta6)));
        dGNHEfdpHf  = -beta5*exp(beta5*(pH(N)-beta6))*GNHEf^2;
        dGNHEbdpHb  = -beta5*exp(beta5*(pH(1)-beta6))*GNHEb^2;
        jNHEf = aNHE0f*GNHEf*R*T*(log(cNa0f/cNa(N))-(pH(N)-pH0f)*log(10));
        jNHEb = aNHE0b*GNHEb*R*T*(log(cNa0b/cNa(1))-(pH(1)-pH0b)*log(10));
        djNHEfdpHf  = -aNHE0f*GNHEf*R*T*log(10) + jNHEf/GNHEf*dGNHEfdpHf;
        djNHEbdpHb  = -aNHE0b*GNHEb*R*T*log(10) + jNHEb/GNHEb*dGNHEbdpHb;
        djNHEfdcNaf = -aNHE0f*GNHEf*R*T/cNa(N);
        djNHEbdcNab = -aNHE0b*GNHEb*R*T/cNa(1);


        % Cl-HCO3 Exchanger
        GAE2f = 1/(1+exp(-beta7*(pH(N)-beta8)));
        GAE2b = 1/(1+exp(-beta7*(pH(1)-beta8)));
        dGAE2fdpHf  = beta7*exp(-beta7*(pH(N)-beta8))*GAE2f^2;
        dGAE2bdpHb  = beta7*exp(-beta7*(pH(1)-beta8))*GAE2b^2;
        jAE2f = aAE20f*GAE2f*R*T*(log(cCl0f/cCl(N))+(pH(N)-pH0f)*log(10));
        jAE2b = aAE20b*GAE2b*R*T*(log(cCl0b/cCl(1))+(pH(1)-pH0b)*log(10));
        djAE2fdpHf  = aAE20f*GAE2f*R*T*log(10) + jAE2f/GAE2f*dGAE2fdpHf;
        djAE2bdpHb  = aAE20b*GAE2b*R*T*log(10) + jAE2b/GAE2b*dGAE2bdpHb;
        djAE2fdcClf = -aAE20f*GAE2f*R*T/cCl(N);
        djAE2bdcClb = -aAE20b*GAE2b*R*T/cCl(1);

        Jn = zeros(N_sp,N-1);
        for in = 1:N_sp
            Jn(in,:) = -Dn(in)*(cn(2:N,in)-cn(1:N-1,in))/dx ...
                + vc*(cn(2:N,in)+cn(1:N-1,in))/2 ...
                - Dn(in)*Zn(in)*F/R/T*(cn(2:N,in)+cn(1:N-1,in))/2.*(phi(2:N)-phi(1:N-1))/dx;
        end

        %% Equations for pc and the Derivatives
        Fn(s_pc) = -pc(2)+pc(1) + dx*eta*thetan(1)*(vn(1)-vc);
        DF(s_pc,[s_pc,s_pc+1]) = [1, -1];
        DF(s_pc,s_vc) = -dx*eta*thetan(1);
        DF(s_pc,s_vn) = dx*eta*thetan(1);
        DF(s_pc,s_tn) = dx*eta*(vn(1) - vc);

        Fn(s_pc+1:s_pc+N-2) = -pc(3:N)+pc(1:N-2) ...
            + 2*dx*eta*thetan(2:N-1).*(vn(2:N-1)-vc);
        for i = 2:N-1
            DF(s_pc+i-1,[s_pc+i-2, s_pc+i]) = [1, -1];
            DF(s_pc+i-1,s_vc)     = -2*dx*eta*thetan(i);
            DF(s_pc+i-1,s_vn+i-1) =  2*dx*eta*thetan(i);
            DF(s_pc+i-1,s_tn+i-1) = -2*dx*eta*(vc - vn(i));
        end

        Fn(s_pc+N-1) = -alphab/alphaf*(pc(1)-p0b) - (pc(N)-p0f)...
            + (dgf-alphab/alphaf*dgb)*(vc + v0) ...
            + alphab/alphaf*R*T*(cin(1)-c0b) + R*T*(cin(N)-c0f);

        DF(s_pc+N-1,[s_pc,s_pc+N-1]) = -[alphab/alphaf, 1];
        DF(s_pc+N-1,s_vc) = dgf - alphab/alphaf*dgb;
        DF(s_pc+N-1,[s_Na, s_Na+N-1])  = R*T*[alphab/alphaf, 1];
        DF(s_pc+N-1,[s_K , s_K+N-1])   = R*T*[alphab/alphaf, 1];
        DF(s_pc+N-1,[s_Cl, s_Cl+N-1])  = R*T*[alphab/alphaf, 1];
        DF(s_pc+N-1,[s_pH, s_pH+N-1])  = R*T*[alphab/alphaf*dcdpH(1), dcdpH(N)];
        DF(s_pc+N-1,[s_Buf,s_Buf+N-1]) = R*T*[alphab/alphaf*dcdcBuf(1), dcdcBuf(N)];
        DF(s_pc+N-1,[s_A , s_A+N-1])   = R*T*[alphab/alphaf, 1];
        DF(s_pc+N-1,s_v0) = dgf - alphab/alphaf*dgb;

        %% Equations for vc and the Derivatives
        Fn(s_vc) = -(pc(N)-p0f) + (dgf + 1/alphaf)*vc + dgf*v0 ...
            + R*T*(cin(N)-c0f);

        DF(s_vc,s_pc+N-1) = -1;
        DF(s_vc,s_vc)     = dgf + 1/alphaf;
        DF(s_vc,s_Na+N-1) = R*T;
        DF(s_vc,s_K+N-1)  = R*T;
        DF(s_vc,s_Cl+N-1) = R*T;
        DF(s_vc,s_pH+N-1) = R*T*dcdpH(N);
        DF(s_vc,s_Buf+N-1)= R*T*dcdcBuf(N);
        DF(s_vc,s_A+N-1)  = R*T;
        DF(s_vc,s_v0) = dgf;

        %% Equations for vn and the Derivatives
        Fn(s_vn) = - (sigma(2)-sigma(1)) + dx*eta*thetan(1)*(vc-vn(1)) ...
            - dx*etast(1)*thetan(1)*(vn(1) + v0);
        DF(s_vn,s_vc) = dx*eta*thetan(1);
        DF(s_vn,s_vn) = -dx*eta*thetan(1) - dx*etast(1)*thetan(1);
        DF(s_vn,s_tn) = dsigmadtn + dx*eta*(vc-vn(1)) - dx*etast(1)*(vn(1)+v0);
        DF(s_vn,s_tn+1) = -dsigmadtn;
        DF(s_vn,s_v0) = -dx*etast(1)*thetan(1);

        Fn(s_vn+1:s_vn+N-2) = -(sigma(3:N)-sigma(1:N-2)) ...
            + 2*dx*eta*thetan(2:N-1).*(vc-vn(2:N-1)) ...
            - 2*dx*etast(2:N-1).*thetan(2:N-1).*(vn(2:N-1)+v0);
        for i = 2:N-1
            DF(s_vn+i-1,s_vc) = 2*dx*eta*thetan(i);
            DF(s_vn+i-1,s_vn+i-1) = -2*dx*thetan(i)*(eta+etast(i));
            DF(s_vn+i-1,[s_tn+i-2,s_tn+i]) = [1,-1]*dsigmadtn;
            DF(s_vn+i-1,s_tn+i-1) = 2*dx*eta*(vc-vn(i)) - 2*dx*etast(i)*(vn(i)+v0);
            DF(s_vn+i-1,s_v0) = -2*dx*etast(i)*thetan(i);
        end

        Fn(s_vn+N-1) = thetan(N)*vn(N) + Jactinf;
        DF(s_vn+N-1,s_vn+N-1) = thetan(N);
        DF(s_vn+N-1,s_tn+N-1) = vn(N);
        DF(s_vn+N-1,s_tc+N-1) = DJactinfDtcN;

        %% Equations for thetan and the Derivatives
        Fn(s_tn) = thetan(1)*vn(1);
        DF(s_tn,s_vn) = thetan(1);
        DF(s_tn,s_tn) = vn(1);

        Fn(s_tn+1:s_tn+N-2) = (thetan(3:N).*vn(3:N)-thetan(1:N-2).*vn(1:N-2))...
            + 2*dx*gamma(2:N-1).*thetan(2:N-1);
        for i = 2:N-1
            DF(s_tn+i-1,[s_vn+i-2,s_vn+i]) = [-thetan(i-1), thetan(i+1)];
            DF(s_tn+i-1,s_tn+i-2) = -vn(i-1);
            DF(s_tn+i-1,s_tn+i-1) = 2*dx*gamma(i) + 2*dx*dgammadtn(i)*thetan(i);
            DF(s_tn+i-1,s_tn+i)   =  vn(i+1);
        end

        Fn(s_tn+N-1) = (thetan(N)*vn(N)-thetan(N-1)*vn(N-1)) ...
            + dx*gamma(N)*thetan(N);
        DF(s_tn+N-1,[s_vn+N-2,s_vn+N-1]) = [-thetan(N-1), thetan(N)];
        DF(s_tn+N-1,s_tn+N-2) = -vn(N-1);
        DF(s_tn+N-1,s_tn+N-1) =  vn(N) + dx*gamma(N) + dx*thetan(N)*dgammadtn(N);

        %% Equations for thetac and the Derivatives
        Fn(s_tc) = thetac(1)*vc - Dtc/dx*(thetac(2)-thetac(1));
        DF(s_tc,s_vc)   = thetac(1);
        DF(s_tc,s_tc)   = vc + Dtc/dx;
        DF(s_tc,s_tc+1) = - Dtc/dx;

        Fn(s_tc+1:s_tc+N-2) = vc * (thetac(3:N)-thetac(1:N-2)) ...
            - 2*Dtc/dx*(thetac(1:N-2) - 2*thetac(2:N-1) + thetac(3:N))...
            - 2*dx*gamma(2:N-1).*thetan(2:N-1);
        for i = 2:N-1
            DF(s_tc+i-1,s_vc) = (thetac(i+1) - thetac(i-1));
            DF(s_tc+i-1,s_tn+i-1) = -2*dx*gamma(i) - 2*dx*thetan(i)*dgammadtn(i);
            DF(s_tc+i-1,s_tc+i-2) = -vc - 2*Dtc/dx;
            DF(s_tc+i-1,s_tc+i-1) = 4*Dtc/dx;
            DF(s_tc+i-1,s_tc+i)   =  vc - 2*Dtc/dx;
        end

        Fn(s_tc+N-1) = 1/2*(sum(thetan(1:N-1)+thetac(1:N-1)) ...
            + sum(thetan(2:N)+thetac(2:N)))...
            - L*(Thetan + Thetac)/dx;
        DF(s_tc+N-1,[s_tn,s_tn+N-1]) = 1/2;
        DF(s_tc+N-1,s_tn+1:s_tn+N-2) = 1;
        DF(s_tc+N-1,[s_tc,s_tc+N-1]) = 1/2;
        DF(s_tc+N-1,s_tc+1:s_tc+N-2) = 1;

        %% Equations for cNa and the Derivatives
        n_ion = n_Na;
        JNa = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJNa = zeros(N-1,N_var);
        for in = 1:N-1
            DJNa(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJNa(in,s_Na+in-1)  =  Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJNa(in,s_Na+in)    = -Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJNa(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJNa(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_Na) = JNa(1) - jpb(n_ion) - jNKEb - jNHEb;
        DF(s_Na,s_pc)  = - dTmbdpcb*jpb(n_ion)/Tmb;
        DF(s_Na,s_vc)  = - dTmbdvc*jpb(n_ion)/Tmb;
        DF(s_Na,s_tn)  = - dTmbdtnb*jpb(n_ion)/Tmb;
        DF(s_Na,s_v0)  = - dTmbdv0*jpb(n_ion)/Tmb;
        DF(s_Na,s_Na)  = - djpbdcnb(n_ion) - djNKEbdcNab - djNHEbdcNab;
        DF(s_Na,s_K)   = - djNKEbdcKb;
        DF(s_Na,s_pH)  = - djNHEbdpHb;
        DF(s_Na,s_phi) = - djpbdphib(n_ion) - djNKEbdphib;
        DF(s_Na,:) = DF(s_Na,:) + DJNa(1,:);

        Fn(s_Na+1:s_Na+N-2)   = JNa(1:N-2)    - JNa(2:N-1);
        DF(s_Na+1:s_Na+N-2,:) = DJNa(1:N-2,:) - DJNa(2:N-1,:);

        Fn(s_Na+N-1) = JNa(N-1) + jpf(n_ion) + jNKEf + jNHEf;
        DF(s_Na+N-1,s_pc+N-1)  = dTmfdpcf*jpf(n_ion)/Tmf;
        DF(s_Na+N-1,s_vc)      = dTmfdvc*jpf(n_ion)/Tmf;
        DF(s_Na+N-1,s_tn+N-1)  = dTmfdtnf*jpf(n_ion)/Tmf;
        DF(s_Na+N-1,s_v0)      = dTmfdv0*jpf(n_ion)/Tmf;
        DF(s_Na+N-1,s_Na+N-1)  = djpfdcnf(n_ion) + djNKEfdcNaf + djNHEfdcNaf;
        DF(s_Na+N-1,s_K+N-1)   = djNKEfdcKf;
        DF(s_Na+N-1,s_pH+N-1)  = djNHEfdpHf;
        DF(s_Na+N-1,s_phi+N-1) = djpfdphif(n_ion) + djNKEfdphif;
        DF(s_Na+N-1,:) = DF(s_Na+N-1,:) + DJNa(N-1,:);

        %% Equations for cK and the Derivatives
        n_ion = n_K;
        JK = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJK = zeros(N-1,N_var);
        for in = 1:N-1
            DJK(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJK(in,s_K+in-1)   =  Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJK(in,s_K+in)     = -Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJK(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJK(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_K) = JK(1) - jpb(n_ion) + 2/3*jNKEb;
        DF(s_K,s_pc)  = - dTmbdpcb*jpb(n_ion)/Tmb;
        DF(s_K,s_vc)  = - dTmbdvc*jpb(n_ion)/Tmb;
        DF(s_K,s_tn)  = - dTmbdtnb*jpb(n_ion)/Tmb;
        DF(s_K,s_v0)  = - dTmbdv0*jpb(n_ion)/Tmb;
        DF(s_K,s_Na)  = 2/3*djNKEbdcNab;
        DF(s_K,s_K)   = - djpbdcnb(n_ion)  + 2/3*djNKEbdcKb;
        DF(s_K,s_phi) = - djpbdphib(n_ion) + 2/3*djNKEbdphib;
        DF(s_K,:) = DF(s_K,:) + DJK(1,:);

        Fn(s_K+1:s_K+N-2)   = JK(1:N-2)    - JK(2:N-1);
        DF(s_K+1:s_K+N-2,:) = DJK(1:N-2,:) - DJK(2:N-1,:);

        Fn(s_K+N-1) = JK(N-1) + jpf(n_ion) - 2/3*jNKEf;
        DF(s_K+N-1,s_pc+N-1)  = dTmfdpcf*jpf(n_ion)/Tmf;
        DF(s_K+N-1,s_vc)      = dTmfdvc*jpf(n_ion)/Tmf;
        DF(s_K+N-1,s_tn+N-1)  = dTmfdtnf*jpf(n_ion)/Tmf;
        DF(s_K+N-1,s_v0)      = dTmfdv0*jpf(n_ion)/Tmf;
        DF(s_K+N-1,s_Na+N-1)  = - 2/3*djNKEfdcNaf;
        DF(s_K+N-1,s_K+N-1)   = djpfdcnf(n_ion)  - 2/3*djNKEfdcKf;
        DF(s_K+N-1,s_phi+N-1) = djpfdphif(n_ion) - 2/3*djNKEfdphif;
        DF(s_K+N-1,:) = DF(s_K+N-1,:) + DJK(N-1,:);

        %% Equations for cCl and the Derivatives
        n_ion = n_Cl;
        JCl = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJCl = zeros(N-1,N_var);
        for in = 1:N-1
            DJCl(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJCl(in,s_Cl+in-1)  =  Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJCl(in,s_Cl+in)    = -Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJCl(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJCl(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_Cl) = JCl(1) - jpb(n_ion) - jAE2b;
        DF(s_Cl,s_pc)  = - dTmbdpcb*jpb(n_ion)/Tmb;
        DF(s_Cl,s_vc)  = - dTmbdvc*jpb(n_ion)/Tmb;
        DF(s_Cl,s_tn)  = - dTmbdtnb*jpb(n_ion)/Tmb;
        DF(s_Cl,s_v0)  = - dTmbdv0*jpb(n_ion)/Tmb;
        DF(s_Cl,s_Cl)  = - djpbdcnb(n_ion) - djAE2bdcClb;
        DF(s_Cl,s_pH)  = - djAE2bdpHb;
        DF(s_Cl,s_phi) = - djpbdphib(n_ion);
        DF(s_Cl,:) = DF(s_Cl,:) + DJCl(1,:);

        Fn(s_Cl+1:s_Cl+N-2)   = JCl(1:N-2)    - JCl(2:N-1);
        DF(s_Cl+1:s_Cl+N-2,:) = DJCl(1:N-2,:) - DJCl(2:N-1,:);

        Fn(s_Cl+N-1) = JCl(N-1) + jpf(n_ion) + jAE2f;
        DF(s_Cl+N-1,s_pc+N-1)  = dTmfdpcf*jpf(n_ion)/Tmf;
        DF(s_Cl+N-1,s_vc)      = dTmfdvc*jpf(n_ion)/Tmf;
        DF(s_Cl+N-1,s_tn+N-1)  = dTmfdtnf*jpf(n_ion)/Tmf;
        DF(s_Cl+N-1,s_v0)      = dTmfdv0*jpf(n_ion)/Tmf;
        DF(s_Cl+N-1,s_Cl+N-1)  = djpfdcnf(n_ion) + djAE2fdcClf;
        DF(s_Cl+N-1,s_pH+N-1)  = djAE2fdpHf;
        DF(s_Cl+N-1,s_phi+N-1) = djpfdphif(n_ion);
        DF(s_Cl+N-1,:) = DF(s_Cl+N-1,:) + DJCl(N-1,:);


        %% Equations for pH and the Derivatives
        n_ion = n_HCO3;
        JHCO3 = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJHCO3 = zeros(N-1,N_var);
        for in = 1:N-1
            DJHCO3(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJHCO3(in,s_pH+in-1)  = dcHCO3dpH(in)*(Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHCO3(in,s_pH+in)    = dcHCO3dpH(in+1)*(-Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHCO3(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJHCO3(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        n_ion = n_Buf;
        JBuf  = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJBuf = zeros(N-1,N_var);
        for in = 1:N-1
            DJBuf(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJBuf(in,s_Buf+in-1)  =  Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJBuf(in,s_Buf+in)    = -Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJBuf(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJBuf(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        n_ion = n_H;
        JH  = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJH = zeros(N-1,N_var);
        for in = 1:N-1
            DJH(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJH(in,s_pH+in-1)  = dcHdpH(in)*(Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJH(in,s_pH+in)    = dcHdpH(in+1)*(-Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJH(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJH(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_pH) = JHCO3(1) + JBuf(1) - JH(1) + jAE2b - jNHEb;
        DF(s_pH,s_Na)  = -djNHEbdcNab;
        DF(s_pH,s_Cl)  =  djAE2bdcClb;
        DF(s_pH,s_pH)  =  djAE2bdpHb - djNHEbdpHb;
        DF(s_pH,:) = DF(s_pH,:) + DJHCO3(1,:) + DJBuf(1,:) - DJH(1,:);

        Fn(s_pH+1:s_pH+N-2)   = JHCO3(1:N-2) - JHCO3(2:N-1) ...
            + JBuf(1:N-2) - JBuf(2:N-1) - JH(1:N-2) + JH(2:N-1);
        DF(s_pH+1:s_pH+N-2,:) = DJHCO3(1:N-2,:) - DJHCO3(2:N-1,:)...
            + DJBuf(1:N-2,:) - DJBuf(2:N-1,:) - DJH(1:N-2,:) + DJH(2:N-1,:);

        Fn(s_pH+N-1) = JHCO3(N-1) + JBuf(N-1) - JH(N-1) - jAE2f + jNHEf;
        DF(s_pH+N-1,s_Na+N-1)  =  djNHEfdcNaf;
        DF(s_pH+N-1,s_Cl+N-1)  = -djAE2fdcClf;
        DF(s_pH+N-1,s_pH+N-1)  = -djAE2fdpHf + djNHEfdpHf;
        DF(s_pH+N-1,:) = DF(s_pH+N-1,:) + DJHCO3(N-1,:) + DJBuf(N-1,:) - DJH(N-1,:);

        %% Equations for cA and the Derivatives
        n_ion = n_A;
        JA = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJA = zeros(N-1,N_var);
        for in = 1:N-1
            DJA(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJA(in,s_A+in-1)   =  Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJA(in,s_A+in)     = -Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx;
            DJA(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJA(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_A)   = JA(1);
        DF(s_A,:) = DJA(1,:);

        Fn(s_A+1:s_A+N-2)   = JA(1:N-2)    - JA(2:N-1);
        DF(s_A+1:s_A+N-2,:) = DJA(1:N-2,:) - DJA(2:N-1,:);

        Fn(s_A+N-1) = sum(cA(1:N-1)+cA(2:N))/2*dx/L - NA/S/L;
        DF(s_A+N-1,[s_A,s_A+N-1]) = 1/2*dx/L;
        DF(s_A+N-1,s_A+1:s_A+N-2) = dx/L;

        %% Equations for cBuf and the Derivatives
        n_ion = n_HBuf;
        JHBuf = Jn(n_ion,:);   % a 1 by (N-1) vector
        DJHBuf = zeros(N-1,N_var);
        for in = 1:N-1
            DJHBuf(in,s_vc)       = (cn(in+1,n_ion)+cn(in,n_ion))/2;
            DJHBuf(in,s_pH+in-1)  = dcHBufdpH(in)*(Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHBuf(in,s_pH+in)    = dcHBufdpH(in+1)*(-Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHBuf(in,s_Buf+in-1) = 10^(pKB-pH(in))*(Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHBuf(in,s_Buf+in)   = 10^(pKB-pH(in+1))*(-Dn(n_ion)/dx + vc/2 ...
                - Dn(n_ion)*Zn(n_ion)*F/R/T/2*(phi(in+1)-phi(in))/dx);
            DJHBuf(in,s_phi+in-1) =  Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
            DJHBuf(in,s_phi+in )  = -Dn(n_ion)*Zn(n_ion)*F/R/T...
                *(cn(in+1,n_ion)+cn(in,n_ion))/2/dx;
        end

        Fn(s_Buf) = JBuf(1) + JHBuf(1);
        DF(s_Buf,:) = DJBuf(1,:) + DJHBuf(1,:);

        Fn(s_Buf+1:s_Buf+N-2)   = JBuf(1:N-2) - JBuf(2:N-1) ...
            + JHBuf(1:N-2) - JHBuf(2:N-1);
        DF(s_Buf+1:s_Buf+N-2,:) = DJBuf(1:N-2,:) - DJBuf(2:N-1,:)...
            + DJHBuf(1:N-2,:) - DJHBuf(2:N-1,:);

        Fn(s_Buf+N-1) = (sum(cBuf(1:N-1)+cBuf(2:N))...
            +sum(cHBuf(1:N-1)+cHBuf(2:N)))/2*dx/L - NBufNHBuf/L/S;
        DF(s_Buf+N-1,[s_pH,s_pH+N-1]) = [dcHBufdpH(1),dcHBufdpH(N)]/2*dx/L;
        DF(s_Buf+N-1,s_pH+1:s_pH+N-2) = dcHBufdpH(2:N-1)*dx/L;
        DF(s_Buf+N-1,[s_Buf,s_Buf+N-1]) = [1+10^(pKB-pH(1)),1+10^(pKB-pH(N))]/2*dx/L;
        DF(s_Buf+N-1,s_Buf+1:s_Buf+N-2) = (1+10.^(pKB-pH(2:N-1)))*dx/L;


        %% Equations for phi and the Derivatives
        for in = 1:N
            Fn(s_phi+in-1) = sum(Zn.*cn(in,:));
            DF(s_phi+in-1,s_Na+in-1) = Zn(n_Na);
            DF(s_phi+in-1,s_K+in-1)  = Zn(n_K);
            DF(s_phi+in-1,s_Cl+in-1) = Zn(n_Cl);
            DF(s_phi+in-1,s_pH+in-1) = Zn(n_H)*dcHdpH(in) + Zn(n_HCO3)*dcHCO3dpH(in);
            DF(s_phi+in-1,s_A+in-1)  = Zn(n_A);
            DF(s_phi+in-1,s_Buf+in-1) = Zn(n_Buf);
        end

        %% Equation for v0 and the derivatives
        Fn(s_v0) = (fextf-fextb) + (p0f-p0b) + (dgf + dgb)*(vc + v0) + kad*v0 ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1).*vn(1:N-1)) ...
            + sum(etast(2:N).*thetan(2:N).*vn(2:N)))...
            + dx/2*v0*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));
        DF(s_v0,s_vc) = dgf + dgb;
        DF(s_v0,[s_vn,s_vn+N-1]) = dx/2*[etast(1)*thetan(1),etast(N)*thetan(N)];
        DF(s_v0,s_vn+1:s_vn+N-2) = dx*etast(2:N-1).*thetan(2:N-1);
        DF(s_v0,[s_tn,s_tn+N-1]) = dx/2*[etast(1)*(vn(1)+v0),etast(N)*(vn(N)+v0)];
        DF(s_v0,s_tn+1:s_tn+N-2) = dx*etast(2:N-1).*(vn(2:N-1)+v0);
        DF(s_v0,s_v0) = (dgf + dgb) + kad ...
            + dx/2*(sum(etast(1:N-1).*thetan(1:N-1)) + sum(etast(2:N).*thetan(2:N)));


        temp_Fn = Fn;

        %% Solve for the matrix
        DF = sparse(DF);
        X = X - DF\Fn;

        if sum(isnan(X)) >= 1 || iter == Iter || norm(temp_Fn) < norm(Fn)
            Solution = 0;
            pc     = NaN(N,1);
            vc     = NaN;
            vn     = NaN(N,1);
            thetan = NaN(N,1);
            thetac = NaN(N,1);
            cNa    = NaN(N,1);
            cK     = NaN(N,1);
            cCl    = NaN(N,1);
            pH     = NaN(N,1);
            cA     = NaN(N,1);
            cBuf   = NaN(N,1);
            phi    = NaN(N,1);
            v0     = NaN;
            break
        else
            Solution = 1;
        end

        NaNcount = isnan(X);
        if sum(sum(NaNcount)) >= 1
            Solution = 0;
        end

        pc = X(s_pc:s_pc+N-1);
        vc = X(s_vc);
        vn = X(s_vn:s_vn+N-1);
        thetan = X(s_tn:s_tn+N-1);
        thetac = X(s_tc:s_tc+N-1);
        cNa  = X(s_Na:s_Na+N-1);
        cK   = X(s_K:s_K+N-1);
        cCl  = X(s_Cl:s_Cl+N-1);
        pH   = X(s_pH:s_pH+N-1);
        cA   = X(s_A:s_A+N-1);
        cBuf = X(s_Buf:s_Buf+N-1);
        phi  = X(s_phi:s_phi+N-1);
        v0   = X(s_v0);

        if iter > 1
            error = abs((X-temp_X)./(X+eps));
            error = sum(error)/(N_var);
            if error < 1d-6 || iter == Iter
                ITER = false;
            end
        end
        temp_X = X;
    end

    V0(ih) = v0;
end

%% Plotting

figure(1)
X = categorical({'0.6 SWELL1','1 SWELL1','1.6 SWELL1'});
X = reordercats(X,{'0.6 SWELL1','1 SWELL1','1.6 SWELL1'});
bar(X, V0*3.6, 0.6)
set(gca,'fontsize',15);
ylabel('v_0 ({\mu}m/h)','fontsize',15)
box off