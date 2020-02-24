function rho = AtmDensityEarth( h )

% Analytic approximation to the 1976 Standard Atmosphere model for density.
% Valid from 0 to 50,000 meters.
%
% Inputs:
%   h     Altitude above sea-level. (meters) 
%
% Outputs
%   rho   Atmospheric density. (kg/m^3)

a = [0.226065450463884; ...
  -0.0938349282722235; ...
  -2.5040823920558e-06; ...
  3.05957469796736e-11 ];

rho = exp(a(1)+a(2)*h/1e3+a(3)*h.^2/1e3+a(4)*h.^3/1e3);