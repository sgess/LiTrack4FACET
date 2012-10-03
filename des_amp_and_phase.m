function MACH = des_amp_and_phase()
% MACH = des_amp_and_phase() retrieves kylstron amplitude and phase
% information with calculated feedback phases and LEM fudge factors

global PARAM;

lattice   = PARAM.MACH.LTC;   % Phasing lattice, "uniform" or "decker"

E0        = PARAM.ENRG.E0;    % GeV ... initial energy
E1        = PARAM.ENRG.E1;    % GeV ... energy at LBCC
E2        = PARAM.ENRG.E2;    % GeV ... energy at FACET

TOTL_phas = PARAM.MACH.RAMP;  % Phase ramp
LONE_phas = PARAM.LONE.PHAS;  % 2-10 phase, should be ~-21 for lattice = unifrom, ~-11.5 for lattice = decker
LTWO_phas = PARAM.LTWO.PHAS;  % 11-20 phase
LONE_gain = PARAM.LONE.GAIN;  % Energy gain in 2-10, set automatically if 0

MACH.SECT.Z    = [101.6000; % Start LI02
                  203.2000; % Start LI03
                  304.8000; % Start LI04
                  406.4000; % Start LI05
                  508.0000; % Start LI06
                  609.6000; % Start LI07
                  711.2000; % Start LI08
                  812.8000; % Start LI09
                  914.4000; % Start LI10
                  1016.000; % Start LI11
                  1117.600; % Start LI12
                  1219.200; % Start LI13
                  1320.800; % Start LI14
                  1422.400; % Start LI15
                  1524.000; % Start LI16
                  1625.600; % Start LI17
                  1727.200; % Start LI18
                  1828.800; % Start LI19
                  1930.400];% Start LI20
           
MACH.SECT.LEFF = [75.47370; % Length LI02
                  96.01200; % Length LI03
                  84.01050; % Length LI04
                  97.41120; % Length LI05
                  97.41120; % Length LI06
                  97.41120; % Length LI07
                  97.41120; % Length LI08
                  97.41120; % Length LI09
                  66.97020; % Length LI10
                  82.78580; % Length LI11
                  96.53650; % Length LI12
                  97.41120; % Length LI13
                  97.41120; % Length LI14
                  97.41120; % Length LI15
                  97.41120; % Length LI16
                  97.41120; % Length LI17
                  97.41120; % Length LI18
                  84.36010];% Length LI19
              
MACH.SECT.ENLD = [1.594234; % E gain, no load LI02
                  1.852992; % E gain, no load LI03
                  1.717310; % E gain, no load LI04
                  1.880000; % E gain, no load LI05
                  1.880000; % E gain, no load LI06
                  0.000000; % E gain, no load LI07
                  0.000000; % E gain, no load LI08
                  1.880000; % E gain, no load LI09
                  1.438510; % E gain, no load LI10
                  1.492416; % E gain, no load LI11
                  1.863119; % E gain, no load LI12
                  1.880000; % E gain, no load LI13
                  1.880000; % E gain, no load LI14
                  1.880000; % E gain, no load LI15
                  1.880000; % E gain, no load LI16
                  1.880000; % E gain, no load LI17
                  1.880000; % E gain, no load LI18
                  1.628119];% E gain, no load LI19

MACH.SECT.PHAS  = zeros(18,1);
MACH.SECT.AMPL  = zeros(18,1);

if strcmp(lattice,'uniform')
    
    % Allow 2-10 egain to set as parameter for comparing to MJH file
    % If 0, calculate using 2-10 phase
    if LONE_gain == 0
        Egain = E1 - E0;
        Eampl = Egain/cosd(LONE_phas+TOTL_phas);
    else
        Eampl = LONE_gain;
    end
    
    % Sum 2-10 cavity length
    L0210 = sum(MACH.SECT.LEFF(1:9));
    MACH.SECT.PHAS(1:9) = LONE_phas+TOTL_phas;
    % Allocate egain by length
    MACH.SECT.AMPL(1:9) = Eampl*MACH.SECT.LEFF(1:9)/L0210; 
    
    % Energy gain for 11-19 calculated with phase
    Egain = E2 - E1;
    Eampl = Egain/cosd(LTWO_phas+TOTL_phas);
    
    % Sum 2-10 cavity length
    L1120 = sum(MACH.SECT.LEFF(10:18));
    MACH.SECT.PHAS(10:18) = LTWO_phas+TOTL_phas;
    % Allocate egain by length
    MACH.SECT.AMPL(10:18) = Eampl*MACH.SECT.LEFF(10:18)/L1120; 
    
end

if strcmp(lattice,'decker')
    
    % 2-10 staged phasing
    MACH.SECT.PHAS(1) = 0 + TOTL_phas;           % LI02 PDES = 0
    MACH.SECT.PHAS(2) = LONE_phas + TOTL_phas;   % LI03 PDES ~ -11.5
    MACH.SECT.PHAS(3) = 2*LONE_phas + TOTL_phas; % LI04 PDES ~ -23
    MACH.SECT.PHAS(4) = 3*LONE_phas + TOTL_phas; % LI05 PDES ~ -34.5
    MACH.SECT.PHAS(5) = 3*LONE_phas + TOTL_phas; % LI06 PDES ~ -34.5
    MACH.SECT.PHAS(6) = -90;                     % LI07 OFF
    MACH.SECT.PHAS(7) = -90;                     % LI08 OFF
    MACH.SECT.PHAS(8) = 0 + TOTL_phas;           % LI09 PDES = 0
    MACH.SECT.PHAS(9) = 5*LONE_phas + TOTL_phas; % LI10 PDES ~ -57.5
    
    % add up energy in 2-10                                               
    NOLEM = sum(MACH.SECT.ENLD(1:9).*cosd(MACH.SECT.PHAS(1:9)));
    % fractional contribution of sectors to total energy
    FRAC  = MACH.SECT.ENLD(1:9).*cosd(MACH.SECT.PHAS(1:9))/NOLEM;
    % NOLEM is greater than E1 - E0 because "all" klystrons are on
    Egain = E1 - E0;
    DIFF  = NOLEM - Egain;
    % LEM will be subtracted off sector by sector, such that LEM*cos(phi)
    % is equal to the proportional energy gain of the sector
    LEM   = FRAC*DIFF./cosd(MACH.SECT.PHAS(1:9));
    LEM(6:7) = 0;
    MACH.SECT.AMPL(1:9) = MACH.SECT.ENLD(1:9) - LEM(1:9);
   
    % Energy gain for 11-19 calculated with phase
    Egain = E2 - E1;
    Eampl = Egain/cosd(LTWO_phas+TOTL_phas);
    
    % Sum 2-10 cavity length
    L1120 = sum(MACH.SECT.LEFF(10:18));
    MACH.SECT.AMPL(10:18) = Eampl*MACH.SECT.LEFF(10:18)/L1120;
    MACH.SECT.PHAS(10:18) = LTWO_phas;
    
    % Energy gain for 11-19 calculated with phase
    Egain = E2 - E1;
    Eampl = Egain/cosd(LTWO_phas+TOTL_phas);
    
    % Sum 2-10 cavity length
    L1120 = sum(MACH.SECT.LEFF(10:18));
    MACH.SECT.PHAS(10:18) = LTWO_phas+TOTL_phas;
    % Allocate egain by length
    MACH.SECT.AMPL(10:18) = Eampl*MACH.SECT.LEFF(10:18)/L1120;
    
end              