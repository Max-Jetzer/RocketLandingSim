%This function is designed to tell you the thrust at a given time 't' for a
%CTI Pro98-6g N1100W commercial solid rocket motor. It takes a solitary
%input 't' and outputs a number telling you the force in newtons. It will
%return NaN if the input is less than zero or greater than 12.19.
function T = BoosterThrust(t)
% 0	0
% 0.1600	2624
% 0.3300	2708
% 0.9100	2055
% 1.2200	1896
% 2.4400	1793
% 3.6600	1625
% 4.8800	1402
% 6.1200	1158
% 7.4100	854
% 9.7700	494
% 12.1800	111.2000
% 12.1900	0
ti=[0; 0.16; 0.33; 0.91; 1.22; 2.44; 3.66; 4.88; 6.12; 7.41; 9.77; 12.18; 12.19];
Ti=[0; 2624; 2708; 2055; 1896; 1793; 1625; 1402; 1158; 854; 494; 111.2; 0];
T=interp1(ti,Ti,t);
end