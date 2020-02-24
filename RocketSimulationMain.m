%============================================================
%Two-Stage, High-Altitude Rocket Landing Footprint Calculator
%============================================================
%
%This program may or may not eventually outline the probability of landing
%footprints for the 2018 exhibition High-Altitude flight of the Rocket Team
%at the University of Minnesota.
%
%Written by Lord King of Holynesse, Daniel "the Maniel" Toth; 
%Hise Royal Alchemyste, Maxwell Houston "We Have A Problem" Jetzer;
%and the Court Git Wizard, Archmaester Ryan Bowers.
%
%Requires the following files:
%   AtmDensityEarth.m
%   BoosterThrust.m
%   HighAltAscentRHS.m
%   HighAltDescentRHS.m
%   HorizWindsAtAltitude.m
%   LoadNOAAWindData.m
%   plotLandingSite.m
%   SustainerThrust.m
%   SpaceportAmerica.jpg
%
%NOTE: Uncomment lines 180-182 to generate iterative plots of each
%simulated flight trajectory at the expense of slowing the simulation
tic
%% Master Simulation Parameters
numSim = 1000;
footprintx = zeros(1,numSim);
footprinty = zeros(1,numSim);
%% Misc. Constants, Wind Data
Re = 6378.14; %Radius of Earth in km
global ground;
ground = 1400; %Ground elevation above sea level (m)
mu = 398600.4; %Gravitational Constant of Earth in km^3/s^2
windData = LoadNOAAWindData('ELP');
format long;
%% Ascent Parameters
%Booster stage
stage(1).Cd = .65; %Coefficient of Drag
stage(1).A = pi*.25*.102^2; %Cross-section area (m^2)
stage(1).M = 20.627; %Loaded mass
stage(1).Mp = 6.79;  %Propellent mass
stage(1).Isp = 210.3;%Specific impulse
stage(1).uE = 9.81*stage(1).Isp; %exhaust velocity
stage(1).T = @BoosterThrust; %function for thrust curve
stage(1).burntime = 12.19; %motor burntime
stage(1).delay = 0; %time before stage ignites

%Sustainer stage
stage(2).Cd = .22;
stage(2).A = pi*.25*.058^2;
stage(2).M = 5.794;
stage(2).Mp = 2.938;
stage(2).Isp = 186.1;
stage(2).uE = 9.81*stage(2).Isp;
stage(2).T = @SustainerThrust;
stage(2).burntime = 4.079;
stage(2).delay = 10;

%Separate stage for coast (somehow this seems neater)
stage(3) = stage(2);
stage(3).M = stage(2).M-stage(2).Mp; %mass is constant, just structure
stage(3).T = @(t) 0; %No thrust during coast
stage(3).delay = 0; %No delay either
stage(3).burntime = 5*60; %Above max possible time to apogee

%% Initial Ascent Conditions
x0 = 0;
y0 = 0;
z0 = ground; %x-y origin is (0,0); z starts at ground level
v0 = .1;     %Need nonzero initial velocity to avoid NaNs
m0 = stage(1).M; %Initial mass is first stage mass

%% Descent Parameters
D1 = 0.1; % Diameter of drogue parachute in m
D2 = 1.5;  %Diameter of main paracute in m
Bd1 = ((stage(3).M)/(1.4*(pi*(D1/2)^2))); %Ballistic coefficient of drogue stage
Bd2 = ((stage(3).M)/(1.4*(pi*(D2/2)^2))); %Ballistic coefficient of main stage
rhoFun = @(h) AtmDensityEarth(h);

%% Begin simulations
for n=1:numSim
    %Inject randomness into initial flight path angle
    gammaOffset=abs(randn/20000);
    if gammaOffset>.035 %Eliminate outliers
        gammaOffset=abs(randn/100000);
    end
    gamma0 = pi/2-gammaOffset; %Subtract from neutral FPA to get the right direction
    
    windFun = @(h) HorizWindsAtAltitude(windData,h); %INPUT WIND THINGY
    
    %Simulate weathercocking by making rocket travel some random amount
    %into the wind. The distribution is centered about the true wind
    %direction and normally distributed with a standard deviation of 1/25
    windDir = windFun(1410);
    if(windDir(1)>0)
        domainFix = pi;
    else %make sure that arctangent will output true angle
        domainFix = 0;
    end
    
    ang = atan(windDir(2)/windDir(1)) + domainFix;
    
    theta = ang+randn/25;%Assign now-randomized heading angle
    
    %% Ascent Simulation

    %Variables for holding ascent data
    masterT = []; masterS = [];
    %Simulate each stage in sequence
    for i = 1:length(stage)
        rhs = @(t,s) HighAltAscentRHS(mu,Re,AtmDensityEarth(s(3)),theta,stage(i),t,s(1),s(2),s(3),s(4),s(5),s(6));
        opts = odeset();
        if i==1 %Stage one initial conditions are given above
            s0 = [x0;y0;z0;v0;gamma0;m0];
        else %Each successive stage starts at the end of the last
            s0 = masterS(end,:)';
            s0(6) = stage(i).M; %Also drop structural mass when staging
            if i==length(stage) %Check for apogee only on coast
                opts = odeset('Events',@apogeeEvent);
            end
        end
        %ode45 call because non-stiff ODE
        [descentT,sout] = ode45(rhs,[0, stage(i).delay+stage(i).burntime],s0,opts); 
        if i==1
            masterT = descentT; %Assign first value
            masterS = sout;
        else %Concatenate stage results to ascent data
            masterT = [masterT;descentT+masterT(end)];
            masterS = [masterS;sout]; %assign additional values
        end
    end
    
    
    %% Descent setup
    %Initial conditions from ascent
    Vxi = masterS(end,4)*cos(theta);
    Vyi = masterS(end,4)*sin(theta);
    
    %Randomness for North wind and East wind
    windrand = randn(3)*0.25;
    if abs(windrand(1))>0.25 %Again, cut off outliers
        while(abs(windrand(1))>0.25)
            %In this case, reroll outliers until they are not outliers
            windrand(1) = randn*0.25; 
        end
    end
    if abs(windrand(2)) > 0.25
       while(abs(windrand(2))>0.25)
            windrand(2) = randn*0.25;
        end %These will get passed into the RHS dynamic functions
    end
    %% Descent function handles for the "right hand side" dynamics
    rhsd1 = @(t,s) HighAltDescentRHS( s(1), s(2), s(3), s(4), s(5), s(6), mu, Bd1, windFun, rhoFun, Re, windrand ); %STAGE 1
    rhsd2 = @(t,s) HighAltDescentRHS( s(1), s(2), s(3), s(4), s(5), s(6), mu, Bd2, windFun, rhoFun, Re, windrand ); %STAGE 2
    
    %% Descent simulation
    %STAGE 1 (STARTING AT APOGEE)
    
    %Call ode45 using final ascent conditions
    %(ode45 is ok here because the system is not THAT stiff)
    %Check altitude to determine when to open main chute
    opts = odeset('Events',@descentStageEvent);
    [descentT,descentS] = ode45( rhsd1, [0, 40*60], [masterS(end,1); masterS(end,2); masterS(end,3); Vxi; Vyi; 0],opts );
    
    %Stage Two (STARTING AT STAGE 1 END)
    %We used ode15s here because our system is very stiff due to a big
    %parachute
    
    %Stop simulation when rocket hits ground
    opts = odeset('Events',@descentLandingEvent);
    [tout,sout] = ode15s( rhsd2, [0 40*60], descentS(end,:), opts );
    
    %Add main chute data to descent simulation variables
    descentT = [descentT; tout];
    descentS = [descentS; sout];
    
    %UNCOMMENT FOR FUN GRAPHS (Plots all 3d Paths, do it because it's cool)
    figure(5);plot3(masterS(:,1),masterS(:,2),masterS(:,3)-ground, 'b',descentS(:,1),descentS(:,2),descentS(:,3)-ground, 'r'); hold on; grid on; title('All 3D Paths');
    %figure(6);plot(sqrt((masterS(:,1)).^2+(masterS(:,2)).^2),masterS(:,3)-ground, 'b',sqrt((descentS(:,1)).^2+(descentS(:,2)).^2),descentS(:,3)-ground, 'r'); hold on; grid on; title('height vs displacement');
    %figure(7);plot(masterS(:,1),masterS(:,2), 'b',descentS(:,1),descentS(:,2), 'r'); hold on; grid on; title('x vs y positions');
    
    %Add landing location to a new vector for statistics!
    footprintx(n) = descentS(end,1);
    footprinty(n) = descentS(end,2);   
end
% END SIMULATION, BEGIN ANALYSIS AND PLOTTING

%% Statistics 'n Stuff

%First we find the mean
meanx = mean(footprintx);
meany = mean(footprinty);
% %Then the Standard Deviation
% sdx = std(footprintx);
% sdy = std(footprinty);
% 
% %Create an ellipse to represent one standard deviations from the mean
% ellipsextop = linspace(meanx-sdx,meanx+sdx,300);
% ellipsexbottom = linspace(meanx+sdx,meanx-sdx,300);
% ellipseytop = meany + real(sqrt(sdy^2*(1-(ellipsextop-meanx).^2/(sdx^2))));
% ellipseybottom = meany - real(sqrt(sdy^2*(1-(ellipsexbottom-meanx).^2/(sdx^2))));
% ellipsex = [ellipsextop,ellipsexbottom];
% ellipsey = [ellipseytop,ellipseybottom];
% 
% %Create an ellipse to represent two standard deviations from the mean
% ellipsextop2 = linspace(meanx-2*sdx,meanx+2*sdx,300);
% ellipsexbottom2 = linspace(meanx+2*sdx,meanx-2*sdx,300);
% ellipseytop2 = meany + real(sqrt((2*sdy)^2*(1-(ellipsextop2-meanx).^2/((2*sdx)^2))));
% ellipseybottom2 = meany - real(sqrt((2*sdy)^2*(1-(ellipsexbottom2-meanx).^2/((2*sdx)^2))));
% ellipsex2 = [ellipsextop2,ellipsexbottom2];
% ellipsey2 = [ellipseytop2,ellipseybottom2];
% %There is probably a built-in funciton for ellipses but this was more fun
data = [footprintx' footprinty'];

covariance = cov(data);
[eigenvec, eigenval] = eig(covariance);
% Get the index of the largest eigenvector
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% Get the largest eigenvalue
largest_eigenval = max(max(eigenval));

% Get the smallest eigenvector and eigenvalue
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the coordinates of the data mean
avg = mean(data);

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=avg(1);
Y0=avg(2);
a=sqrt(largest_eigenval);
b=sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );
ellipse_xr2 = 2*a*cos(theta_grid);
ellipse_yr2 = 2*b*sin(theta_grid);
%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
r2_ellipse = [ellipse_xr2;ellipse_yr2]'*R;

%% "I watch it for the plots"
%make ALL the plots!
%Ascent graphs
figure(1),clf
plot(masterT,masterS(:,3)-ground); title('Rocket Altitude - Last Flight'); xlabel('Time(s)'); ylabel('Altitude AGL(m)');
figure(2),clf
plot(masterT,(180/pi)*abs(pi/2-masterS(:,5))); title('FPA deviation from vertical - Last Flight'); xlabel('Time(s)'); ylabel('FPA deviation (degrees)');

%Descent Graphs
figure(3),clf
plot(masterS(:,1),masterS(:,2),'b'); title('X-Y Plane Projected Path - Last Flight');hold on;
plot(descentS(:,1),descentS(:,2),'m'); legend('Ascent','Descent','location','best');
xlabel('Downrange East(m)'); ylabel('Downrange North(m)');
%plot3(descentS(:,1),descentS(:,2),descentS(:,3));


%Overlay landing data on image of launch site
figure(4),clf
launchSite = imread('SpaceportAmerica.jpg');
image(launchSite), axis equal, hold on
plotLandingSite(0,0,'rs')
plotLandingSite(footprintx,footprinty,'y.')
plotLandingSite(meanx,meany,'kx')
plotLandingSite(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'b'),
plotLandingSite(r2_ellipse(:,1) + X0,r2_ellipse(:,2) + Y0,'r');
plotLandingSite([-1e4, -1e4+1e3],[8e3 8e3],'g','LineWidth',3)
title('Landing Site Footprint')
legend('Launch Site','Landing Sites','Mean Landing Site','One Std. Dev.','Two Std. Dev.','Scale: 1km','location','best')
set(gca,'xTickLabel',[])
set(gca,'yTickLabel',[])
toc
%% Event functions for ODEs
%Checks for apogee (gamma levels off to 0 or pi)
function [value,isterminal,direction] = apogeeEvent(t,s)
    value = sin(s(5));
    isterminal = 1;
    direction = -1;
end
%Checks for main chute deploy (z = 1 km AGL)
function [value,isterminal,direction] = descentStageEvent(t,s)
    global ground;
    value = s(3)-(ground+1000);
    isterminal = 1;
    direction = -1;
end
%Checks for touchdown (z = 0 AGL)
function [value,isterminal,direction] = descentLandingEvent(t,s)
    global ground;
    value = s(3)-ground;
    isterminal = 1;
    direction = -1;
end