%Checks for main chute deploy (z = 1 km AGL)
function [value,isterminal,direction] = descentStageEvent(t,s)
    global ground;
    value = s(3)-(ground+1000);
    isterminal = 1;
    direction = -1;
end