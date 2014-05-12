# Solitary gravity wave computation
===================================

These Matlab scripts compute the steady irrotational surface solitary gravity wave solution of the Euler equations (homogeneous, incompressible and perfect fluids). The wave is defined by its Froude number *Fr* and the result is about fifteen digits accurate. The method works for all but the highest waves, *i.e.* for all amplitude/depth ratio less than 0.796.

## Synopsis:

* SolitaryGravityWave(Fr,[],1); % plot results only 
* [zs,ws,fs,SWP] = SolitaryGravityWave(Fr); % output results at the surface and parameters 
* [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z); % surface and bulk output 
* [zs,ws,fs,SWP,W,F,P,A] = SolitaryGravityWave(Fr,Z,1);

## Input: 
* Fr : Froude number (must be a scalar). 
* Z : Complex abscissa where fields are desired inside the fluid (default Z = []). 
:      Z should be strictly below the surface, i.e., -1 <= imag(Z) < eta(real(Z)) 
:      y = eta(x) being the equation of the free surface. 
* PF : Plot Flag. If PF=1 the final results are plotted, if PF~=1 nothing is plotted (default).

## Output (dimensionless quantities): 
* zs  : Complex abscissa at the surface, i.e., x + i*eta. 
* ws  : Complex velocity at the surface, i.e., u - i*v. 
* fs  : Complex potential at the surface, i.e., phi + i*psi. 
* SWP : Solitary Wave Parameters, i.e. 
        + SWP(1) = wave amplitude, max(eta) 
        + SWP(2) = wave mass 
        + SWP(3) = circulation 
        + SWP(4) = impulse 
        + SWP(5) = kinetic energy 
        + SWP(6) = potential energy 
* W : Complex velocity in the bulk at abscissas Z. 
* F : Complex potential in the bulk at abscissas Z. 
* P : Pressure in the bulk at abscissas Z. 
* A : Complex acceleration in the bulk at abscissas Z (A = dW / dt).

## Example: 
zs = SolitaryGravityWave(1.25);

plot(real(zs),imag(zs))