function plotLandingSite(x,y,varargin)
%plotLandingSite plots the given x-y coordinates on the map
%   plotLandingSite(x,y) converts the x,y values in m from launch site to
%   pixels on the mapimage.jpg figure to display their approximate location
%   on a real map. Function assumes current figure already displays the
%   map. Any additional parameters that could be passed to plot(), such as
%   color and line specs, work identically for plotLandingSite.

origin = [1455 1111]; %Pixel count for launch site
runway = [775,228;865,688]; %Endpoints of Spaceport America runway
runway_m = 3657; %Known real distance of runway; used as standard ruler
runway_pix = norm(runway(2,:)-runway(1,:)); %Find runway pixel length
m2pix = runway_pix/runway_m; %conversion factor between m and pixels

%plot converted points on map
plot(origin(1)+x*m2pix,origin(2)-y*m2pix,varargin{:})
end
