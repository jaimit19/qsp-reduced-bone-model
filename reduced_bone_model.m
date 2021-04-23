%-------------------------------------------------------------------------
% This code is based on the final reduced model 
% in Hasegawa, C. and Duffull, S.B. (2018), Automated Scale Reduction of 
% Nonlinear QSP Models With an Illustrative Application to a Bone Biology 
% System. CPT Pharmacometrics Syst. Pharmacol., 7: 562-572. 
% https://doi.org/10.1002/psp4.12324
% This "lumped" reduced model is based on the larger models of Peterson 
% MC, Riggs MM (2010) Bone 46:49-63 and Peterson MC, Riggs MM (2012) CPT 
% Pharmacometrics Syst Pharmacol 1:e14, which in turn are based on prior
% mathematical models of others for calcium homeostasis and bone biology. 
%-------------------------------------------------------------------------
% First draft of code by Helen Moore April 18, 2021, based on the NONMEM
% mod file provided by Hasegawa and Duffull and available on GitHub:
% https://github.com/metrumresearchgroup/OpenBoneMin/tree/master/inst/community/lump
% Code further modified by Jaimit Parikh to generate the preliminary
% results
%-------------------------------------------------------------------------
% Ref [1] Hasegawa and Duffull (2018)
% Ref [2] Peterson and Riggs   (2010)
% Ref [3] Peterson and Riggs   (2012)

%%% Parameter values
% Params 1-33 from Table 1 in Ref [1]
% Param 34 calculated from Eq. (1) in Ref [3], with parameter name R_BMD
%  used instead of k_inBMD, and d_BMD instead of k_out,BMD, to match the 
%  to match the naming conventions for params 1-33. Also, BSAP (related 
%  to OB) and CTx (related to OC) must be assumed 100% of baseline in  
%  order to match the results of Ref. [1]. To estimate R_BMD, BMD_LS 
%  (assumed to be the same as BMD in Ref. [1]) is assumed to have initial
%  value equal to 100% of baseline. Thus R_BMD = d_BMD*100 (see Ref. [3]).
%  In the end, just used d_BMD*100 in ODE #9, rather than R_BMD. 
% Params 35-37 from below Eq. (1) in Ref [3]

parameters = boneModelPars();
IC = getInitialConditions();
scale = 365*24; t0 = 0; tfinal = 4*scale;
[T,Y] = ode45(@(t, y)odefun(t, y, parameters),[t0 tfinal],IC);


f = figure('DefaultAxesFontSize',16, 'Position', [20 20 1800 900]);
subplot(3, 2, 1);
plot(T / scale,  (Y(:, 9) - 100), 'color', 'k', 'LineWidth', 2);
xlabel('Time, Years'); ylabel('% change BMD');
subplot(3, 2, 2);
plot(T /scale, drugConc(T), 'color', 'k', 'LineWidth', 2); 
xlabel('Time, Years'); ylabel('[Drug]');
subplot(3, 2, 3);
plot(T /scale, 100 * (1 + (Y(:,5) - IC(5))./ IC(5)), 'color', 'k', 'LineWidth', 2); 
xlabel('Time, Years'); ylabel('TGFb');
subplot(3, 2, 4);
changeBSAP = 100 * (1 + ((Y(:,7) + Y(:,8)) - (IC(7) + IC(8))) ./ (IC(7) + IC(8)));
plot(T / scale, changeBSAP, 'color', 'k','LineWidth', 2); 
xlabel('Time, Years'); ylabel('BSAP');
subplot(3, 2, 5);
plot(T /scale, Y(:, 4), 'color', 'k', 'LineWidth', 2); 
xlabel('Time, Years'); ylabel('OC');
subplot(3, 2, 6);
plot(T / scale, Y(:, 3), 'color', 'k','LineWidth', 2); 
xlabel('Time, Years'); ylabel('CMX');


function C = drugConc(T)
%C = 0 + (T - T); 
C = 50 + 30*(sin(0.0008*T));
end

function p = boneModelPars()
p.R_L1        = 75.0;     % "production rate of L1";
p.k_OBL1      = 55.3;   % "rate constant for effect of OB on L1 productn";
p.k_L2L1      = 160;      % 3  rate constant from L2 to L1
p.k_CMXL      = 0.112;    % 4  rate constant from CMX to L1 or L2
p.d_L1        = 0.970;    % 5  degradation rate constant of L1
p.R_L2        = 0.000337; % 6  production rate of L2
p.k_OBRANKL   = 0.234;    % 7  rate constant for effect of OB on RANKL prodn
p.d_L2        = 0.00110;  % 8  degradation rate const of L2
p.d_RANKL     = 0.00290;  % 9  degradation rate const of RANKL
p.k_int       = 0.00795;  % 10 eliminatn rate const of denosumab-RANKL cmplx
p.K_SS        = 63.4;     % 11 steady-state deno-RANKL binding affinity ng/mL
p.k_L2CMX     = 0.0000190;% 12 rate constant from L2 to CMX
p.R_OC        = 0.00000298;%13 production rate of OC
p.d_OC        = 0.0898;   % 14 degradation rate constant of OC
p.a_1         = 2.18;     % 15 max response of TGF to degradation of OC
p.rho_1       = 0.200;    % 16 min response of TGF to degradation of OC
p.del_1       = 16.2;     % 17 TGF amt producing half-max resp to OC degradn
p.gam_1       = 1;        % 18 sigmoidicity term for TGF effect on OC degrdn
p.a_2         = 3.80;     % 19 max response of CMX to degradation of OC
p.rho_2       = 0.470;    % 20 min response of CMX to degradation of OC
p.del_2       = 0.000013; % 21 CMX amt producing half-max resp to OC degradn
p.gam_2       = 3.09;     % 22 sigmoidicity term for CMX effect on OC degrdn
p.k_OCTGF     = 5.66;     % 23 rate const for effect of OC on TGF production
p.d_TGF       = 0.0298;   % 24 degradation rate constant of TGF
p.R_ROB       = 0.000003; % 25 production rate of ROB
p.a_3         = 4.18;     % 26 max response of TGF to production of ROB
p.rho_3       = 0.202;    % 27 min response of TGF to production of ROB
p.del_3       = 34.0;     % 28 TGF amt producing half-max resp to ROB prodn
p.gam_3       = 1.81;     % 29 sigmoidicity term for TGF effect on ROB prodn
p.k_ROBOB     = 0.003;    % 30 rate constant from ROB to OB
p.f           = 0.06;     % 31 fraction converting from ROB to SOB
p.d_FOB       = 0.01;     % 32 degradation rate constant of FOB
p.d_SOB       = 0.000001; % 33 degradation rate constant of SOB
%var         = 0.382;   %BMD additive resid error var; only for nlme in [1]
p.R_BMD       = 0.0146;   % 34 BMD production rate - ended up not using
p.d_BMD       = 0.000146; % 35 BMD degradation rate constant (1/h)
p.gam_OB      = 0.0739;   % 36 sigmoidicity term for OB effect on BMD prodtn
p.gam_OC      = 0.0679;   % 37 sigmoidicity term for OC effect on BMD degradn
p.OC0 = 0.0012;
p.FOB0 = 0.0040;
p.SOB0 = 0.0010;
p.C = 0;
p.BMD0 = 100;
end

function IC = getInitialConditions()
%%% Initial conditions; first 8 from Eqs. (13) - (20) in Ref [1]
% Initial condition for BMD is below Eq. (1) in Ref [3]
L10 = 30082;       % L1
L20 = 1.4001;      % L2
CMX0 = 0.00022283;  % CMX
OC0 = 0.0012;      % OC
TGF0 = 0.2281;      % TGF
ROB0 = 0.0010;      % ROB
FOB0 = 0.0040;      % FOB
SOB0 = 0.0010;      % SOB
BMD0 = 100;         % BMD (%)
IC = [L10, L20, CMX0, OC0, TGF0, ROB0, FOB0, SOB0, BMD0]; 
end

%%% ODE system
function dydt = odefun(t,y,p)
dydt = zeros(9,1); %setting up empty vector for the differential equations
% Equations 1-8 below are from Eqs. (21) - (28) in Ref [1]
% Equation 9 is from Eq. (1) in Ref [3], using R_BMD in place of k_in,BMD 
%   and d_BMD in place of k_out,BMD, to match the conventions in Ref. [1].
%   Also, BSAP (related to OB) and CTx (related to OC) must both be assumed 
%   equal to 100% of baseline at time zero to match results in Ref. [1].

L1 = y(1);     %RANK-related lumped state
L2 = y(2);    %RANKL-related lumped state
CMX = y(3);    %RANK-RANKL complex
OC = y(4);     %osteoclasts
TGF = y(5);    %active TGF-beta
ROB = y(6);    %responding osteoblasts
FOB = y(7);    %fast osteoblasts
SOB = y(8);    %slow osteoblasts
BMD = y(9);    %bone mineral density
%%%% Pharmacokinetics (PK) for denosumab
%ABS = y(10);   %ABS,DEFDOSE    absorption/dosing compartment
%CEN = y(11);   %central compartment
%PER = y(12);   %peripheral compartment
%RAN = y(13) ;   %RANKL levels
CC = drugConc(t);
dL1dt = p.R_L1 + p.k_OBL1 * (FOB + SOB) + p.k_L2L1 * L2 + ...
    p.k_CMXL * CMX - p.d_L1 * L1;

dL2dt = p.R_L2 + p.k_OBRANKL * (FOB + SOB) + p.k_CMXL * CMX - ...
    (p.d_L2 + ((p.k_int - p.d_RANKL) * CC / (p.K_SS + CC)) / 3) * L2; % check for C

dCMXdt = p.k_L2CMX * L2 - p.k_CMXL * CMX;                                        

dOCdt = p.R_OC - p.d_OC * (p.rho_1 + (p.a_1 - p.rho_1) * TGF^p.gam_1 / ...
    (p.del_1 ^ p.gam_1 + TGF^p.gam_1)) * (p.a_2 - (p.a_2 - p.rho_2) * ...
    (CMX / 10)^p.gam_2 / (p.del_2^p.gam_2 + (CMX / 10)^p.gam_2)) * OC;                              

dTGFdt = p.k_OCTGF * OC - p.d_TGF * TGF;                                    

dROBdt = p.R_ROB*(p.rho_3 + (p.a_3 - p.rho_3) * TGF^p.gam_3 / ...
    (p.del_3^p.gam_3 + TGF^p.gam_3)) -p.k_ROBOB * ROB;                                     

dFOBdt = p.k_ROBOB * (1 - p.f) * ROB - p.d_FOB * FOB;                                        

dSOBdt = p.k_ROBOB * p.f * ROB - p.d_SOB * SOB;                                        

dBMDdt = p.d_BMD * p.BMD0 * ((FOB + SOB) / (p.FOB0 + p.SOB0))^p.gam_OB - ...
    p.d_BMD * (OC / p.OC0)^p.gam_OC * BMD; % check IC(9)



dydt(1)  = dL1dt;     %RANK-related lumped state
dydt(2)  = dL2dt;    %RANKL-related lumped state
dydt(3)  = dCMXdt;    %RANK-RANKL complex
dydt(4)  = dOCdt;     %osteoclasts
dydt(5)  = dTGFdt;    %active TGF-beta
dydt(6)  = dROBdt;    %responding osteoblasts
dydt(7)  = dFOBdt;    %fast osteoblasts
dydt(8)  = dSOBdt;    %slow osteoblasts
dydt(9)  = dBMDdt;    %bone mineral density




end
