parameters = boneModelParameters();
IC =boneModelIC();
t0 = 0; tfinal = 4*12;
[T,Y] = ode45(@(t, y)odefun(t, y, parameters),[t0 tfinal],IC);

function p =  boneModelParameters()
p.RIN4 = 0.003;      % THETA1 FIX          ; RIN4        1   (0 .60552)
p.K40 = 0.003;       % THETA2  FIX          ; K40         2
p.FRAC = 0.06;       % THETA3 FIX           ; FRAC        3
p.K10 = 0.01;         % THETA4 FIX           ; K10         4
p.K20 = 0.001;       % THETA5 FIX          ; K20         5
p.RIN3 = 0.00298;  % THETA6 FIX        ; RIN3        6
p.K30 = 0.0111;     % THETA7         ; K30         7   .0292 FIX
p.RIN5 = 0.08;           % THETA8  (0, 0.16);  RIN5        8  ###
p.K50 = 0.0011;     % THETA9  FIX         ; K50         9
p.K75 = 0.1120;     % THETA10 FIX         ; K75         10
p.K60 = 0.0298;     % THETA11 FIX         ; k60         11
p.K57 = 0.000019;  % THETA12 FIX       ; k57         12  as a function of RANK[0]?
p.RIN8 = 75;          % THETA13 FIX            ; RIN8        13
p.K18 = 55.3;         % THETA14  FIX          ; k18         14
p.K58 = 160;          % THETA15 FIX           ; k58         15
p.K80 = 0.97;         % THETA16 FIX           ; k80         16
p.K15 = 0.234;       % THETA17 FIX          ; K15 and K25 17
p.kout_BMD =0 ;     % THETA18  (0, 0.146);          ; kout_BMD    18
p.POWOB = .0739;   %THETA19 FIX         ; POWOB       19  (0 .0739)
p.POWOC = .0779;  % THETA20 FIX         ; POWOC       20  (0 .0779)
p.KSS = 50;              % THETA21 (0, 100) (0 100)           ; KSS         21 ###
p.CL = 3.06;            % THETA22 FIX          ; CL 22
p.V1 = 2490;            % THETA23 FIX          ; V1 23
p.Q = 37.9;              % THETA24 FIX          ; Q  24
p.V2 = 1360;            % THETA25 FIX          ; V2 25
p.KA = 0.212;          %  THETA26 FIX         ; KA 26
p.F1 = 0.638;          %  THETA27 FIX         ; F1 27
end

function IC = boneModelIC()
FOB0 = 0.0040; %IC1
SOB0 = 0.0010; %IC2
OC0 = 0.0012; %IC3
ROB0 = 0.0010; %IC4
L20 = 1.4001; % IC5
TGF0 = 0.2281; %IC6
CMX0 = 0.00022283; %IC7
L10 = 30082; % IC8
BMD0 = 100; %IC9
IC10 = 5; % ?
IC11 = 5; % ?
IC12 = 4; % ?
IC13 = 2; % R00 ?
IC = [FOB0, SOB0, OC0, ROB0, L20, TGF0,...
    CMX0, L10, BMD0, IC10, IC11, IC12, IC13];
end

function dydt = odefun(t,y,p)
dydt = zeros(13,1);

% FOB = y(1);   %fast osteoblasts
% SOB = y(2);    %slow osteoblasts
% OC = y(3);    %osteoclasts
% ROB = y(4);   %responding osteoblasts
% L2 = y(5);    %RANKL-related lumped state
% TGF = y(6);    %active TGF-beta
% CMX = y(7);    %RANK-RANKL complex
% L1 = y(8);    %RANK-related lumped state
% BMD = y(9); % Bone mineral density
% Y10 = y(10);
% Y11 = y(11);
% Y12 = y(12);
% Y13 = y(13);

IC = boneModelIC();
THETA = struct2array(p);

% ; 4 ROB
PIC0 = IC(6);                      %; 0.228142
RIN4 = THETA(1)*24*30/1000;          %; *(IC1+IC2)/PIC0
K40 = THETA(2)*24*30;
E0PICROB = 0.8838*PIC0;
EMAXPICROB = 3.9745;
GAMMA = 1.810;
EC50PICROBPAREN = (EMAXPICROB*IC(6) / (PIC0 - E0PICROB)) - IC(6);
EC50PICROB = exp(log(EC50PICROBPAREN));

%; 1 fast OB
K41 = K40*(1-THETA(3));
K10 = THETA(4)*24*30; %       ; k17aD
K20 = THETA(5)*24*30/1000;

% ; 2 slow OB
K42 = K40*THETA(3);


% ; 3 OC
RIN3 = THETA(6)*24*30/1000;
K30 = THETA(7)*24*30;
E0PICOC = 0.878215*PIC0;
EMAXPICOC = 2;
EC50PICOCPAREN = (EMAXPICOC*IC(6)/(PIC0 - E0PICOC)) - IC(6);
EC50PICOC = exp(log(EC50PICOCPAREN));

E0RANKL = 3.8034;
EMAXL = 0.4698;
PIL0 = 0.000022283;
LSURVOCGAM = 3.0902;
EC50SURVINPAR = (E0RANKL - EMAXL)*(PIL0^LSURVOCGAM/(E0RANKL - 1))...
    - PIL0^LSURVOCGAM;
EC50SURV = exp(log(EC50SURVINPAR)/LSURVOCGAM);


%; 5 RANKL-RELATED LUMPED STATE
RIN5 = THETA(8)*24*30/1000;
K50 = THETA(9)*24*30;
K75 = THETA(10)*24*30;
K15 = THETA(17)*24*30;
K25 = K15;
KSS = THETA(21);


%; 6 ACTIVE TGF_BETA
K60 = THETA(11)*24*30;
K36 = K60*IC(6)/IC(3);

%; 7 RANK-RANKL COMPLEX
K57 = THETA(12)*24*30;         %; as a function of RANK[0]?
K70 = K75;

% 8 RANK-RELATED LUMPED STATE
RIN8 = THETA(13)*24*30;
K18 = THETA(14)*24*30;
K28 = K18;
K58 = THETA(15)*24*30;
K80 = THETA(16)*24*30;
K78 = K75;


%; 9 BMD
KOUT_BMD = THETA(18)*24*30/1000;
POWOB = THETA(19);
POWOC = THETA(20);
BMD0 = 100;

%; 10-13 PK
KDEG = 2; KSSINI = 2; %?
CL   = THETA(22)/1000*24*30; %                  ; L/month
V1   = THETA(23)/1000;
Q    = THETA(24)/1000*24*30;
V2   = THETA(25)/1000;
TVKA = THETA(26)/24*24*30;
KA   = TVKA*(66/71.5)^(-0.577);
F10  = THETA(27);
KSYN = IC(13)*KDEG; %R00



KINT  = 0.00795*24*30;
KOUTL = 0.00290*24*30;

C = 0.5*((y(11)/V1 - y(13) - KSSINI) + sqrt((y(11)/V1 - y(13) - KSSINI)^2 +...
    4*KSSINI*y(11)/V1));
CLTOT = CL + KINT*V1*y(13)/(KSSINI+C);
MIC = (KINT-KDEG)*C/(KSSINI+C);

dydt(10) = -KA*y(10);
dydt(11) = KA*y(10) - CLTOT*C - Q*(C-y(12)/V2);
dydt(12) = Q*(C-y(12)/V2);
dydt(13) = KSYN - (KDEG+MIC)*y(13);

DRUG = (KINT-KOUTL)*C/(KSS+C);

%; H+2016
PICROB = E0PICROB +...
    EMAXPICROB*y(6)^GAMMA/(y(6)^GAMMA...
    + EC50PICROB^GAMMA);

%; H+2018D
PICOC = E0PICOC + EMAXPICOC*y(6)/(y(6) + EC50PICOC);

%; H-2218D
PIL = y(7)/10;
LSURVOC = E0RANKL - (E0RANKL - EMAXL)*(PIL^LSURVOCGAM/(PIL^LSURVOCGAM ...
    + EC50SURV^LSURVOCGAM));

dydt(1) =K41*y(4) - K10*y(1);     %; fast OB 1
dydt(2) =K42*y(4) - K20*y(2);          %                                ; slow OB 2
dydt(3) =RIN3 - K30*PICOC*LSURVOC*y(3); %                               ; OC      3
dydt(4) =RIN4*PICROB - K40*y(4);                %                      ; ROB     4
dydt(5) =RIN5 + K15*y(1) + K25*y(2) + K75*y(7) - (K50+DRUG/3)*y(5); %  ; RANKL-RELATED LUMPED STATE 5
dydt(6) =K36*y(3) - K60*y(6);         %                                ; TGFBETA 6
dydt(7) =K57*y(5) - K70*y(7);           %                              ; RANK-RANKL COMPLEX 7
dydt(8) =RIN8 + K18*y(1) + K28*y(2) + K58*y(5) + K78*y(7) - K80*y(8) ; % RANK-RELATED LUMPED STATE 8

BSAP = y(1) + y(2);
SCTX = y(3);
RIN_BMD = KOUT_BMD*BMD0*(BSAP/(IC(1)+IC(2)))^POWOB;
K_BMD = KOUT_BMD*(SCTX/IC(3))^POWOC;
dydt(9)=RIN_BMD - K_BMD*y(9) ;                                    %   ; BMD 9


end
