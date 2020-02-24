function sDot = HighAltAscentRHS(mu,R,rho,theta,stage,t,x,y,z,v,gamma,m)
%HighAltAscentRHS simulates the dynamics of a rocket in a flat Earth model
%   Returns a column vector sDot containing (in order) the rate of change
%   in positions x, y, z (m); velocity v (m/s); flight path angle gamma
%   (rad); and mass m (kg).
%   Other parameters:
%   mu = gravitational constant (km^3/s^2)
%   R = planet radius (km)
%   rho = atmospheric density (kg/m^3)
%   theta = heading angle (rad, CCW from East)
%   t = time since stage start (s)
%   stage = struct giving data about the current rocket stage, containing
%       T = function handle for motor thrust curve
%       delay = time before stage ignition (s)
%       Cd = Drag Coefficient
%       A = Cross-section area (m^2)
%       uE = exhaust velocity (m/s)

%Initialize sDot
sDot = zeros(5,1);

%Calculate g
g = (mu/(R+z/1e3)^2)*1e3;

%Calculate thrust
T = stage.T(t-stage.delay);
if isnan(T), T=0; end

%Time derivatives of positions are just velocities
sDot(1) = v*cos(gamma+randn/10000)*cos(theta+randn/10000);
sDot(2) = v*cos(gamma+randn/10000)*sin(theta+randn/10000);
sDot(3) = v*sin(gamma+randn/10000);

%Time derivative of velocities given by net force divided by mass
sDot(4) = (T-rho*stage.Cd*v^2*stage.A)/m - g*sin(gamma);

%Time derivative of flight path angle is... this, apparently
sDot(5) = -g*cos(gamma)/v+v*cos(gamma)/(R*1e3+z);

%Time derivative of mass given by thrust and exhaust velocity
sDot(6) = -T/stage.uE;
end
