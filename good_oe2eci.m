function [R,V] = good_oe2eci(a,e,inc,Omega,w,Mo,mu)
% An OE to ECI code that works well

M = mod(Mo,2*pi);
E = M;
d = 1;
while abs(d) > 1e-11
    d = -(E-e*sin(E) - M)/(1-e*cos(E));
    E = E + d;
end
nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));



p = a*(1-e^2);
    
    % CREATING THE R VECTOR IN THE pqw COORDINATE FRAME
    R_pqw(1) = p*cos(nu)/(1 + e*cos(nu));
    R_pqw(2) = p*sin(nu)/(1 + e*cos(nu));
    R_pqw(3) = 0;
    
    % CREATING THE V VECTOR IN THE pqw COORDINATE FRAME    
    V_pqw(1) = -(mu/p)^.5*sin(nu);
    V_pqw(2) =  (mu/p)^.5*(e + cos(nu));
    V_pqw(3) =   0;
 
R1 = [1      0       0;
     0      cos(-inc)  sin(-inc);
     0      -sin(-inc) cos(-inc)];
R3 = [cos(-Omega)  sin(-Omega)     0;
     -sin(-Omega) cos(-Omega)     0;
     0       0          1];
R23 = [cos(-w)  sin(-w)     0;
     -sin(-w) cos(-w)     0;
     0       0          1];
 
    % ROTATING THE pqw VECOTRS INTO THE ECI FRAME (ijk)
size(R1);
size(R3);
size(R23);
size(R_pqw);

    R = R3*R1*R23*R_pqw';
    V = R3*R1*R23*V_pqw';
    
R = R(:);
V = V(:);




end