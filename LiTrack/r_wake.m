function [dE,zc,sigz] = r_wake(z,L,Ne,Nbin,RL);

%
%	Function to return the wakefield induced energy profile vs. z for
%	a set of given axial coordinates "z", for a resistive wake with
%   resistance per unit length of R.
%
%  INPUTS:	z:		The internal axial coordinates, within the bunch, of
%					each electron with respect to any fixed point [m]
%			L:		(Optional, DEF=1 m) The length of the linac [m]
%			Ne:		(Optional, DEF=1  ) The number of electrons in the bunch
%			Nbin:   (Optional, DEF=100) The number of bins to use
%			RL:      resistance per unit length [ohms/m]
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			sigz:	rms bunch length (standard deviation) [m]

%=============================================================================

nn   = length(z);
sigz = std(z);

if nn < 100
  disp(' ')
  disp('Probably too few particles in your input "z" array')            
  disp(' ')
end
if nn > 5E6
  disp(' ')
  disp('>5E6 particles in your "z" array - let''s not get carried away now')
  disp(' ')
end
  
if sigz < 1E-6
  disp(' ')
  disp('Bunch length of <1 micron may be too short for this Green''s function')
  disp(' ')
end

if ~exist('L')
  L = 1; 					% default length of S-band linac [m]
end
if ~exist('Ne')
  Ne = 1; 					% default number of e- in bunch      
end
if ~exist('Nbin')
  Nbin = 100;  				% default number simulation particles
end
if ~exist('RL')
  RL = 1. ;  		       % default resistance
end  



[N,zc] = hist(z,Nbin-2);		% add zero height bins to both ends...
dzc = mean(diff(zc));			% ...so particles in 1/2-end bins get...
zc = [zc(1)-dzc zc zc(Nbin-2)+dzc];	% ...linearly interpolated...
N = [0 N 0];				% ...properly


dE  = zeros(Nbin,1);
for j = 1:Nbin				% jth bin is test-bin
  dE(j) = (-4.8066E-17)*abs(RL)*L*N(j)*(Ne/nn)/dzc;		% test-bin's energy loss
end