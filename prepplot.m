function [x,y,z,props,cax] = prepplot(n)


     
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

x = 6378.1*cosphi*cos(theta);
y = 6378.1*cosphi*sintheta;
z = 6378.1*sin(phi)*ones(1,n+1);
cax = [];

load('topo.mat','topo','topomap1');

topo2 = [topo(:,181:360) topo(:,1:180)]; 

props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo2;
end

