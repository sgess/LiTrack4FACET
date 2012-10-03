  global LINAC; % Linac amplitude and phase
  
  ampl = LINAC.SECT.AMPL;
  phas = LINAC.SECT.PHAS;
  leff = LINAC.SECT.LEFF;

  global PARAM; % Initial beam and machine parameters  
  % FACET energy set points. We are pretty sure about these, I think. . .
  E0        = PARAM.ENRG.E0;    % GeV ... initial energy
  E1        = PARAM.ENRG.E1;    % GeV ... energy at LBCC
  E2        = PARAM.ENRG.E2;    % GeV ... energy at FACET
  
  % NRTL compressor klystron and R56 #s
  NRTL_ampl = PARAM.NRTL.AMPL;  % AMPL DR13 11 VDES
  NRTL_phas = PARAM.NRTL.PHAS;  % on the zero-crossing
  NRTL_leff = PARAM.NRTL.LEFF;  % cavity length
  NRTL_R56  = PARAM.NRTL.R56;   % This is design val
  NRTL_T566 = PARAM.NRTL.T566;  % Design val?
  NRTL_ELO  = PARAM.NRTL.ELO;   % NRTL low energy cut
  NRTL_EHI  = PARAM.NRTL.EHI;   % NRTL high energy cut
  
  % Phase and length of 02-10
  LONE_leff = PARAM.LONE.LEFF;  % Length of LI02-LI10 (m)
  LONE_phas = PARAM.LONE.PHAS;  % Chirp phase
  LONE_gain = PARAM.LONE.GAIN;  % energy gain in LI02-LI10 (GeV)
  LONE_ampl = PARAM.LONE.FBAM;  % feedback amplitude (GV)
  
  % S10 chcn #s
  LI10_R56  = PARAM.LI10.R56;   % Measured val?
  LI10_T566 = PARAM.LI10.T566;  % Measured val?
  LI10_ISR  = PARAM.LI10.ISR;   %
  
  % Energy gain and length of 02-10
  LTWO_leff = PARAM.LTWO.LEFF;  % Length of LI02-LI10 (m)
  LTWO_phas = PARAM.LTWO.PHAS;  % Chirp phase
  LTWO_ampl = PARAM.LTWO.FBAM;  % feedback amplitude (GV)
  
  % S20 chcn R56 #s
  LI20_R56  = PARAM.LI20.R56;   % Measured val?
  LI20_T566 = PARAM.LI20.T566;  % Measured val?
  LI20_ISR  = PARAM.LI20.ISR;   %
  LI20_ELO  = PARAM.LI20.ELO;   % S20 low energy cut
  LI20_EHI  = PARAM.LI20.EHI;   % S20 high energy cut
  
  % Initial beam parameters
  inp       = 'G';		        % gaussian Z and dE/E (see sigz0 =..., sigd0 =...)
  
  % 6mm bunches coming out of the ring teensy energy spread.
  sigz0     = PARAM.INIT.SIGZ0;	% rms bunch length used when inp=G or U above [m]
  sigd0     = PARAM.INIT.SIGD0;	% rms relative energy spread used when inp=G or U above [ ]
  
  % 200K sim particles = 100K electrons per sim particle
  Nesim     = PARAM.INIT.NESIM;	% number of particles to generate for simulation when inp=G or U (reasonable: ~1000 to ~100000)
  
  % The Holtzapple skew. Someday they'll name a skew after me. . .
  asym      = PARAM.INIT.ASYM;	% for inp='M' or 'G': sets rise/fall time width (-1<asym<1)
  
  % Our beam has no tail? That's a tall tale! Jesus I hope no one reads this. . .
  tail      = PARAM.INIT.TAIL;  % for inp='M' or 'G': sets rise/fall time width (0<=tail<1)
  cut       = PARAM.INIT.CUT;   % for inp='G': sets rise/fall time width (0.5<=cut<inf)

  % Other stuff
  splots    = 0;	     % if =1, use small plots and show no wakes (for publish size plots)
  plot_frac = 0.05;      % fraction of particles to plot in the delta-z scatter-plots (0 < plot_frac <= 1)
  Ne        = 2.2e10;	 % number of particles initially in bunch
  z0_bar    = 0;	     % axial offset of bunch [m] (used also with file input - mean of file removed first)
  d0_bar    = 0;	     % relative energy offset of bunch [ ]  (used also with file input - mean of file removed first)
  Nbin      = 200;		 % number of bins for z-coordinate (and dE/E for plots)
  gzfit     = 1;		 % if ==1: fit Z-distribution to gaussian (defaults to no-fit if 'gzfit' not provided)
  gdfit     = 0;		 % if ==1: fit dE/E-distribution to gaussian (defaults to no-fit if 'gdfit' not provided)
  contf     = 0;		 % if ==1: get color contour image of z-d space (defaults to scatter plot if not provided)
  
  % S-band wavelength
  lambdaS   = 2.99792458e8/2856e6;
  
% The follwing array of file names, "wake_fn", is the point-charge wakefield filename(s) to be used.  The pointer
% to the used filename appears in the 5th column (wake ON/OFF) of the 'beamline' array below.  A "zero" (i.e. 0)
% means no wake used, and a value of j (e.g. 1,2,...) directs the calculation to use the jth point-charge wakefield
% file (i.e. the jth row of "wake_fn") for that part of the beamline.

% We always use dat slac
  wake_fn=[ ...
    'slac.dat         '; ...
    'slacx.dat        '; ...
    'SlacL.dat        '; ...
    'Slac_cu_rw.dat   '; ...
    'SS_12700um_rw.dat'; ...
    'Al_12700um_rw.dat'; ...
  ]; % name of point-charge wakefield file(s) ["slac.dat" is default if filename not given]

  comment='FACET (sgess''s best guess)'; % text comment which appears at bottom of plots

%==============|==============================================================================================
%	CODES        |		1				2				3				4				5				6
%==============|==============================================================================================
% compressor   |	code= 6           R56/m          T566/m          E_nom/GeV       U5666/m            -
% chicane      |	code= 7           R56/m         E_nom/GeV           -               -               -
% acceleration |	code=10       tot-energy/GeV    phase/deg        lambda/m   wake(ON=1,2**/OFF=0)  Length/m
% acceleration |	code=11     dEacc(phi=0)/GeV    phase/deg        lambda/m   wake(ON=1,2**/OFF=0)  Length/m
% energy fdbk  |	code=13      energy-goal/GeV    voltage/GV      phase1/deg      phase2/deg        lambda/m
% E-spread add |	code=22          rms_dE/E           -               -               -               -
% E-window cut |	code=25         dE/E_window         -               -               -               -
% E-cut limits |	code=26          dE/E_min        dE/E_max           -               -               -
% con-N E-cut  |	code=27            dN/N             -               -               -               -
% notch-coll   |	code=28          dE/E_min        dE/E_max           -               -               -
% Z-cut limits |	code=36            z_min           z_max            -               -               -
% con-N z-cut  |	code=37            dN/N             -               -               -               -
% modulation   |	code=44           mod-amp        lambda/m           -               -               -
% STOP	       |	code=99             -               -               -               -               -
%=============================================================================================================
% CODE<0 makes a plot here, CODE>0 gives no plot here.
%=============================================================================================================

%=============================================================================================================
% FACET
%=============================================================================================================

beamline=[
      
    % RTL: The energy cut is a little tighter than in the file MJH gave me.
     -1     0		   0          0           0        0         % initial particles
    -11     NRTL_ampl  NRTL_phas  lambdaS     1        NRTL_leff % NRTL compressor
     26     NRTL_ELO   NRTL_EHI   0           0        0         % energy spread cut
     -6     NRTL_R56   NRTL_T566  E0          0        0         % NRTL R56, T566
     
    % LI02: SBST PDES ~ 0
    -11     ampl(1)    phas(1)    lambdaS     1        leff(1)   % LI02:KLYS
      
    % LI03: SBST PDES ~ 11
    -11     ampl(2)    phas(2)    lambdaS     1        leff(2)   % LI03:KLYS
     
    % LI04: SBST PDES ~ 22
    -11     ampl(3)    phas(3)    lambdaS     1        leff(3)   % LI04:KLYS
     
    % LI05: SBST PDES ~ 33
    -11     ampl(4)    phas(4)    lambdaS     1        leff(4)   % LI05:KLYS
     
    % LI06: SBST PDES ~ 44
    -11     ampl(5)    phas(5)    lambdaS     1        leff(5)   % LI06:KLYS
     
    % LI07: SECTOR OFF
    -11     ampl(6)    phas(6)    lambdaS     1        leff(6)   % LI07:KLYS
     
    % LI08: SECTOR OFF
    -11     ampl(7)    phas(7)    lambdaS     1        leff(7)   % LI08:KLYS
    
    % LI09: SBST PDES ~ 0
    -11     ampl(8)    phas(8)    lambdaS     1        leff(8)   % LI09:KLYS
     
    % LI10: SBST PDES ~ 55
    -11     ampl(9)    phas(9)    lambdaS     1        leff(9)   % LI10:KLYS
    -13     E1         LONE_ampl  -90         90       lambdaS   % Energy feedback to set 9GeV in chicane
     
    % S10 CHCKN: MJH file has S10 chicken in two sections but only need one
	  7	    LI10_R56   E1	      0           0        0         % LBCC chicane
     22     LI10_ISR   0          0           0        0         % LBCC ISR energy spread
    -37     0.01       1          0           0        0         % Z-cut
    
    % LI11: SBST PDES ~ 0
    -11     ampl(10)   phas(10)   lambdaS     1        leff(10)  % LI11:KLYS
     
    % LI12: SBST PDES ~ 0
    -11     ampl(11)   phas(11)   lambdaS     1        leff(11)  % LI12:KLYS
        
    % LI13: SBST PDES ~ 0
    -11     ampl(12)   phas(12)   lambdaS     1        leff(12)  % LI13:KLYS
        
    % LI14: SBST PDES ~ 0
    -11     ampl(13)   phas(13)   lambdaS     1        leff(13)  % LI14:KLYS
        
    % LI15: SBST PDES ~ 0
    -11     ampl(14)   phas(14)   lambdaS     1        leff(14)  % LI15:KLYS
        
    % LI16: SBST PDES ~ 0
    -11     ampl(15)   phas(15)   lambdaS     1        leff(15)  % LI16:KLYS
        
    % LI17: SBST PDES ~ 0
    -11     ampl(16)   phas(16)   lambdaS     1        leff(16)  % LI17:KLYS
        
    % LI18: SBST PDES ~ 0
    -11     ampl(17)   phas(17)   lambdaS     1        leff(17)  % LI18:KLYS
        
    % LI19: SBST PDES ~ 0
    -11     ampl(18)   phas(18)   lambdaS     1        leff(18)  % LI19:KLYS
    -13		E2	       LTWO_ampl  -90.0	      90.0     lambdaS   % SCAV energy feedback
    
    % LI20: FACET CHICANE
      6     LI20_R56   LI20_T566  E2          0        0         % LI20 R56, T566
     22     LI20_ISR   0          0           0        0         % FFTB ISR energy spread
     37     0.01       1          0           0        0         % Z-cut
    -99		0          0          0           0        0         % end
  ];

% Sign conventions used:
% =====================
%
% phase = 0 is beam on accelerating peak of RF (crest)
% phase < 0 is beam ahead of crest (i.e. bunch-head sees lower RF voltage than tail)
% The bunch head is at smaller values of z than the tail (i.e. head toward z<0)
% With these conventions, the R56 of a chicane is < 0 (and R56>0 for a FODO-arc) - note T566>0 for both
% 
% * = Note, the Energy-window cut (25-code) floats the center of the given FW window in order center it
%     on the most dense region (i.e. maximum number of particles).
% **  1:=1st wake file (e.g. S-band) is used, 2:=2nd wake file (e.g. X-band)