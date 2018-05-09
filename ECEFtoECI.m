function [reci,veci] = ECEFtoECI(recef,vecef,time,we)
%ECEFTOECI This turns ECEF coordinates into ECI at a given input time
th = 6.7045617/24*2*pi + 360.9856473*pi/180*(time-86400)/86400;
R = [cos(th),sin(th),0;sin(-th),cos(th),0;0,0,1];
reci = inv([R])*recef;
veci = inv([R])*vecef + (cross([0,0,we],recef))';
reci = reci(:);
veci = veci(:);
end

