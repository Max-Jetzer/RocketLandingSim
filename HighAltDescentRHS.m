%This function returns the derivatives of x, y, z, vx, vy, and vz. Those 6
%inputs are fairly self explanatory. B is the ballistic coefficient.
%windFun is a matrix of wind data. rhoFun is a matrix of density data. Re
%is the radius of the earth in km. Windrand is a matrix of data
%representing the randomness in the north and east directions of wind.
function sDot = HighAltDescentRHS( x, y, z, vx, vy, vz,  mu, B, windFun, rhoFun, Re, windrand )

%Description of input/output parameters
% x,y,z,vx,vy,vz are obvious (in meters and m/s)
% m is the mass of the descent stage
% rhoFun is a function handle to give density
% B is the current ballistic coefficient B = m/Cd*A

% atmospheric density and winds

rho = rhoFun(z);
if(z<16150)
    vw = [windFun(z)';0].*(1+windrand);                       %NO WIND ABOVE 53 kft aka 16150m
else                          
    vw = [0; 0; 0];
end

%Initial crunching of parameters
r = [x, y, z];    
v = [vx; vy; vz];  %I put these in vectors for now...

vr = v - vw; %Wind-relative velocity
rmag = norm(r);
vrmag = norm(vr);
g = (mu/(Re+z/1e3)^2)*1e3; %Acceleration due to gravity


% position derivative:
rDot = v;

% velocity derivative:
vxDot = -0.5 * rho * B^-1 * vrmag * vr(1);
vyDot = -0.5 * rho * B^-1 * vrmag * vr(2);
vzDot = -0.5 * rho * B^-1 * vrmag * vr(3) - g;

vDot = [vxDot; vyDot; vzDot];

% full state derivative
sDot = [ rDot; vDot ];