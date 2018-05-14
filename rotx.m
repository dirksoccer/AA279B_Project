function [r2] = rotx(r1,t)
%ROTZ Rotate r1 by t in x direction
rot=[1   0      0      ; ...
     0 cos(t) -sin(t)  ;...
     0 sin(t)  cos(t)] ;
r2 = rot*r1';
end



