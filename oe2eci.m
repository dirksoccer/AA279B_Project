function [veci, reci] = oe2eci(mu,a,e, i, RAN, AP, TA)
%oe in parafocal into ECI coords   
p = a-a*e^2;
r = p/(1+e*cos(TA));
v = sqrt(2*mu/r-mu/a);
rpf = [r*cos(TA) r*sin(TA) 0];
fpa = TA + pi/2 - acos(sqrt(mu*p)/r/v);
vpf = [cos(fpa)*v sin(fpa)*v 0];
reci = rotz(rotx(rotz(rpf,AP)',i)',RAN)';
veci = rotz(rotx(rotz(vpf,AP)',i)',RAN)';
end