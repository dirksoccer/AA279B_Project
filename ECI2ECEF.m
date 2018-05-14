function [recef,vecef] = ECI2ECEF(reci,veci,th,we)
%ECEFTOECI This turns ECI coordinates into ECEF at a given input theta (th)
reci = reci(:);
veci = veci(:);
R = [cos(th),sin(th),0;sin(-th),cos(th),0;0,0,1];
recef = R*reci;
vecef = R*veci - (cross([0,0,we],reci))';
recef = recef(:);
vecef = vecef(:);
end

