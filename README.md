![Solitary wave](http://www.mathworks.fr/matlabcentral/fileexchange/screenshots/8209/original.jpg)

# Solitary gravity wave computation
======

These Matlab scripts compute the steady irrotational surface solitary gravity wave solution of the Euler equations (homogeneous, incompressible and perfect fluids). The wave is defined by its Froude number *Fr* and the result is about fifteen digits accurate. The method works for all but the highest waves, *i.e.* for all amplitude/depth ratio less than 0.796.

## Synopsis:

* SolitaryGravityWave(Fr, [], 1); % *plot results only*
* [zs, ws, fs, SWP] = SolitaryGravityWave(Fr); % *output results at the surface and parameters*
* [zs, ws, fs, SWP, W, F, P, A] = SolitaryGravityWave(Fr,Z); % *surface and bulk output*
* [zs, ws, fs, SWP, W, F, P, A] = SolitaryGravityWave(Fr,Z,1);

## Input: 
* Fr : Froude number (must be a scalar). 
* Z : Complex abscissa where fields are desired inside the fluid (default Z = []). 
  + Z should be strictly below the surface, *i.e.*, -1 <= imag(Z) < eta(real(Z)) 
  + *y = eta(x)* being the equation of the free surface. 
* PF : Plot Flag. If *PF = 1* the final results are plotted, if *PF ~= 1* nothing is plotted (default).

## Output (dimensionless quantities): 
* zs  : Complex abscissa at the surface, *i.e.*, *x + i\*eta*. 
* ws  : Complex velocity at the surface, *i.e.*, *u - i\*v*. 
* fs  : Complex potential at the surface, *i.e.*, *phi + i\*psi*. 
* SWP : Solitary Wave Parameters, i.e. 
  + SWP(1) = wave amplitude, max(eta) 
  + SWP(2) = wave mass 
  + SWP(3) = circulation 
  + SWP(4) = impulse 
  + SWP(5) = kinetic energy 
  + SWP(6) = potential energy 
* W : Complex velocity in the bulk at abscissas *Z*.
* F : Complex potential in the bulk at abscissas *Z*.
* P : Pressure in the bulk at abscissas *Z*.
* A : Complex acceleration in the bulk at abscissas Z (*A = dW/dt*).

## Example: 
zs = SolitaryGravityWave(1.25);
plot(real(zs), imag(zs));

## Parametrization by the amplitude:

The script *SolitaryGravityAmplitude.m* can be used in a similar way with the only difference is that the first parameter is the wave amplitude/depth ratio which has to be less than 0.796, as before.

## References:

More details on the methods used in these scripts can be found in the following references:

* D. Clamond & D. Dutykh. *Fast accurate computation of the fully nonlinear solitary surface gravity waves*. Computers & Fluids, **84**, 35-38, 2013

* D. Dutykh & D. Clamond. *Efficient computation of steady solitary gravity waves*. Wave Motion, **51**, 86-99, 2014

======

D. Dutykh & D. Clamond  
[www.denys-dutykh.com](http://www.denys-dutykh.com/)  
[math.unice.fr/~didierc/](http://math.unice.fr/~didierc/)