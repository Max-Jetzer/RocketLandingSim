function out = HorizWindsAtAltitude( data, h )

% Given the horizontal wind model across altitude, compute the wind speed
% in the East and North directions at a given altitude "h".
%
% USAGE:
% [wE,wN] = HorizWindsAtAltitude( data, h );
%
% Inputs:
%   data    Data structure obtained from the "LoadNOAAWindData.m" function.
%   h       Altitude (meters)
%
% Outputs:
%   wE      Predicted wind speed in the east direction (m/s)
%   wN      Predicted wind speed in the north direction (m/s)

%output as (wE,wN)
out = [interp1(data.altM,data.windEastMPS,h), interp1(data.altM,data.windNorthMPS,h)];


