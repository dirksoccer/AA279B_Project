function [r2] = rotz(r1,t)
%ROTZ Rotate r1 by angle t in z direction
rot = [  cos(t) -sin(t) 0;...
         sin(t) cos(t) 0;...
         0 0 1];
r2 = rot*r1';
end

