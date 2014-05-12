function [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z,PF)
% SOLITARYGRAVITYWAVE: Computes the steady irrotational surface solitary gravity 
%    wave solution of the Euler equations (homogeneous, incompressible and 
%    perfect fluids). The wave is defined by its Froude number Fr and the result 
%    is about fifteen digits accurate. The method works for all but the highest 
%    waves, i.e. for all amplitude/depth ratio less than 0.796.
%
% SYNOPSIS:
% [zs,ws,fs,SWP] = SolitaryGravityWave(Fr);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z);
% [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z,1); 
%
% INPUT:
% Fr : Froude number (i.e., c/sqrt(gd), c.f. Note 1 below).
% Z  : Complex abscissa where fields are desired inside the fluid (default Z = []).
%      Z should be strictly below the surface, i.e., -1 <= imag(Z) < eta(real(Z))
%      y = eta(x) being the equation of the free surface.
% PF : Plot Flag. If PF=1 the final results are plotted, if PF~=1 nothing is 
%      plotted (default).
%
% OUTPUT (dimensionless quantities):
% zs   : Complex abscissa at the surface, i.e., x + i*eta.
% ws   : Complex velocity at the surface, i.e., u - i*v.
% fs   : Complex potential at the surface, i.e., phi + i*psi.
% SWP  : Solitary Wave Parameters, i.e.
%        SWP(1) = wave amplitude, max(eta)
%        SWP(2) = wave mass
%        SWP(3) = circulation
%        SWP(4) = impulse
%        SWP(5) = kinetic energy
%        SWP(6) = potential energy
% W    : Complex velocity in the bulk at abscissas Z.
% F    : Complex potential in the bulk at abscissas Z.
% P    : Pressure in the bulk at abscissas Z.
% A    : Complex acceleration in the bulk at abscissas Z (A = dW / dt).
%
% NOTES:
% 1- Computations are performed in dimensionless units such that rho = g = d = 1,
%    where rho is the fluid density, g is the acceleration due to gravity and d 
%    is the constant water depth. It means that all lengths are scaled by d,
%    accelerations by g, speeds by sqrt(g*d), times by sqrt(d/g) and stresses 
%    by rho*g*d.
% 2- The numerical scheme is based on the Petviashvili's iterations of the
%    Babenko equation using a pseudo-spectral method. 
% 3- The solution is obtained in parametric form resulting from a conformal 
%    mapping into a strip.
% 4- The algorithm is efficient for all but the highest waves.
% 5- W, F, P and A are column vectors corresponding to the shape of Z(). The
%    velocities are computed in the frame of reference where the fluid is
%    at rest in the far field.
% 6- For Z very close to the free surface, the fields may be inaccurate.
%
% Have a look at the m-file for more details.

% WARNING: This program was written to illustrate the Petviashvili method 
%          for the Babenko equation. It is not designed to compute the wave
%          of limiting height. The code is freely distributed since it may be 
%          useful to someone, but it is given 'as such' without any guarantee
%          of any kind. Suggestions for improvements and bug reports are 
%          nonetheless welcome.

% Author 1: Denys Dutykh, CNRS & UCD
% E-mail  : Denys.Dutykh@ucd.ie
% URL     : http://www.denys-dutykh.com/
%
% Author 2: Didier Clamond, University of Nice - Sophia Antipolis

if nargin<2, Z=[]; PF=0; end
if nargin<3, PF=0; end

if Fr>1.29421 || Fr<=1,
   error('The Froude number must be in the interval 1 < Fr <= 1.29421');
end
    
% Preliminaries.
F2 = Fr^2;                                        % Froude number squared
kd = fzero(@(kd) sinc(kd)-F2*cos(kd),[0;1.5]);    % trend parameter
L  = 17*log(10)/kd;                               % half-length of the computational domain

% Numerical parameters (can be changed if you know what you are doing).
N   = 16384;                       % number of Fourier modes (must be even)
tol = 1e-15;                       % iterations stopping criterium
gam = 2.0;                         % the exponent parameter in the Petviashvili method
                                   % (needed to compute the stabilizing factor)

% Discrete variables
dxi = 2*L/N;                       % distance between two consecutive points in the mapped space
xi  = (1-N/2:N/2)'*dxi;            % mapped space discretization
k   = [0:N/2-1,-N/2:-1]'*pi/L;     % wavenumbers
COP = k.*coth(k); COP(1) = 1;      % nonlocal C-operator in the Babenko equation
LOP = F2*COP - 1;                  % linear operator of Babenko eq.
LOI = 1./LOP;                      % its inverse
IOP = -1i*coth(k); IOP(1) = 0;     % nonlocal int-C-operator

% Begin numerical resolution.
eta = (F2-1)*sech(0.5*kd*xi).^2;   % initial guess from KdV approximation
err = inf;
Res_err = [];
Dif_err = [];
while err>tol,                     % Petviashvili's iterations
  % Compute intermediate quantities  
  eta_hat  = fft(eta);
  eta2_hat = fft(eta.^2);
  Ceta     = real(ifft(COP.*eta_hat));
  Leta_hat = LOP.*eta_hat;
  Nl_hat   = 0.5*COP.*eta2_hat + fft(eta.*Ceta);

  % Compute weight.
  S1 = eta_hat'*Leta_hat;   % scalar product
  S2 = eta_hat'*Nl_hat;     % scalar product
  S  = (S1/S2)^gam;
  
  % New value and errors.
  eta1 = real(ifft(S*LOI.*Nl_hat));
  err  = norm(eta1-eta,inf);        % errors are measured in l_inf norm
  eta  = eta1;

  % Store errors and residuals.
  Dif_err = [Dif_err; err];
  Res_err = [Res_err; norm(real(ifft(Leta_hat-Nl_hat)),inf)];
end
% End of numerical resolution.

% Post processing.
Ceta   = real(ifft(COP.*fft(eta)));               % C(eta)
dexi   = real(ifft(1i*k.*fft(eta)));              % d eta / d xi
SWP(1) = max(eta);                                % amplitude
SWP(2) = dxi*eta'*(1+Ceta);                       % mass
SWP(3) = Fr*dxi*sum(Ceta);                        % circulation
SWP(4) = Fr*SWP(2);                               % impulse
SWP(5) = 0.5*F2*dxi*eta'*Ceta;                    % kinetic energy
SWP(6) = 0.5*dxi*(eta.^2)'*(1+Ceta);              % potential energy

% Physical variables at the surface.
etaMean = mean(eta);
xs   = (1 + etaMean)*xi + real(ifft(IOP.*fft(eta-etaMean))); 
phis = Fr*(xs - xi);
qs   = (1+Ceta).^2 + dexi.^2;
us   = Fr - Fr*(1+Ceta)./qs;
vs   = -Fr*dexi./qs;

% Output complex variables at the surface.
zs = xs + 1i*eta;
ws = us - 1i*vs;
fs = phis + 1i*eta*Fr;

% Some information about the computations:
Niter = length(Dif_err);
fprintf('+-------------------------------------------------+\n');
fprintf('| Froude Number           = %15.14f      |\n', Fr);
fprintf('| Amplitude               = %15.14f      |\n', SWP(1));
fprintf('| Mass                    = %15.14f      |\n', SWP(2));
fprintf('| Circulation             = %15.14f      |\n', SWP(3));
fprintf('| Impulse                 = %15.14f      |\n', SWP(4));
fprintf('| Kinetic Energy          = %15.14f      |\n', SWP(5));
fprintf('| Potential Energy        = %15.14f      |\n', SWP(6));
fprintf('|                                                 |\n');
fprintf('| Convergence achieved in %05.0f iterations.       |\n', Niter);
fprintf('| Error between two latest iterations: %5.4e |\n', err);
fprintf('| Residual                           : %5.4e |\n', Res_err(end));
fprintf('+-------------------------------------------------+\n');

% Output at desired locations in the bulk.
% The code below is not fully vectorize in order to avoid building a huge 
% matrix when Z is large.
if isempty(Z)==0,
   Z = Z(:);
   LZ = length(Z);
   W = zeros(LZ,1); F=W; P=W; dW=W;
   dzsm1 = Ceta + 1i*dexi;
   for n=1:LZ,
      W(n)  = sum(dzsm1./(zs-Z(n)) - conj(dzsm1)./(conj(zs)-2i-Z(n))); 
      W(n)  = 1i*0.5*Fr/pi*W(n)*dxi;
      dW(n) = sum(dzsm1./(zs-Z(n)).^2 - conj(dzsm1)./(conj(zs)-2i-Z(n)).^2); 
      dW(n) = 1i*0.5*Fr/pi*dW(n)*dxi;
      F(n) = sum(dzsm1.*log((zs+1i)./(zs-Z(n))) - conj(dzsm1.*log((zs+1i)./(zs+2i-conj(Z(n)))))); 
      F(n) = 1i*0.5*Fr/pi*F(n)*dxi;
   end
   P = Fr*real(W) - 0.5*abs(W).^2 - imag(Z); % pressure
   A = dW.*(conj(W) - Fr);                   % acceleration
end

% Plot the results if requested.
if (PF),
   Lc = 2*SWP(2)/SWP(1); % characteristic length
   
   subplot(3,2,1)
   plot(xs, eta, 'b-','LineWidth',1)
   title('Free surface elevation', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$\eta\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 0 1.05*SWP(1)])

   subplot(3,2,2)
   plot(xs, phis, 'b-','LineWidth',1)
   title('Potential at the free surface', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$\phi\ /\ d\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(phis) 1.05*max(phis)])
   
   subplot(3,2,3)
   plot(xs,us, 'b-','LineWidth',1)
   title('Horizontal velocity at the surface', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$u\ /\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 0 1.05*max(us)])
   
   subplot(3,2,4)
   plot(xs,vs, 'b-','LineWidth',1)
   title('Vertical velocity at the surface', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('$x\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$v\ /\ \sqrt{\, g\/d\ }$', 'interpreter', 'latex', 'fontsize', 14)
   axis([-Lc Lc 1.05*min(vs) 1.05*max(vs)])

   subplot(3,2,5)
   semilogy(1:Niter, Dif_err, 'k.-')
   title('Error between two iterations', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('Iterations', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$|| \eta_{n+1} - \eta_n ||_\infty\ /\ d$', 'interpreter', 'latex', 'fontsize', 14)
   axis([0 Niter+1 min(Dif_err) max(Dif_err)])

   subplot(3,2,6)
   semilogy(1:Niter, Res_err, 'k.-')
   title('Residual error', 'interpreter', 'latex', 'fontsize', 14)
   xlabel('Iterations', 'interpreter', 'latex', 'fontsize', 14)
   ylabel('$|| R(\eta_n) ||_\infty$', 'interpreter', 'latex', 'fontsize', 14)
   axis([0 Niter+1 min(Res_err) max(Res_err)])
   set(gcf, 'Color', 'w')   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=sinc(x)
% SINC(X) is SIN(X)/X  with the singularity at zero removed.

%	Drea Thomas	1992
z=(x~=0);
y=x;
y(z) = sin(x(z))./x(z);
y(~z)=ones(sum(sum(~z)),1);