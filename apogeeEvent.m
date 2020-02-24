%Checks for apogee (gamma levels off to 0 or pi)
function [value,isterminal,direction] = apogeeEvent(t,s)
    value = sin(s(5));
    isterminal = 1;
    direction = -1;
end