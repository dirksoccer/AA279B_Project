function [yout] = differinertial(tout,y,mu)
%differ finding y dot for a propagating orbit calculator

yout(1:3) = y(4:6);% first 3 ydots are the velocities
yout(4:6) = -mu*(y(1:3))./norm(y(1:3))^3;
yout = yout(:);
%columnize that sucker
end

