%Checks for touchdown (z = 0 AGL)
function [value,isterminal,direction] = descentLandingEvent(t,s)
    global ground;
    value = s(3)-ground;
    isterminal = 1;
    direction = -1;
end