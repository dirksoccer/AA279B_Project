


%% Global Orbital Launch Defense
%   Shane Dirks and Nick Folz
%   AA279B Final Project
%   Spring of 2018
%
%   DEPENDENCIES:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% Define Constants
    R_earth = 6378.1; %(km)
    omega_earth = 0.0000729211585530;%(rad/sec) - per HW3
    mu_earth = 398600.4418; %(km^2/s^2)
%     time_sim = 45.4*86400; %seconds (days*86400) since epoch
%     tlength = 1*3600; %seconds of simulation

%     tvec = time_sim:1:time_sim+tlength;
%     thetavec = mod((280.4606 + 360.9856473*tvec/86400)/180*pi,2*pi);
    size_for_things = 300;

    %target things
    target_alt = 1000; %(km)
    target_ratio = (R_earth+target_alt)/R_earth;
    max_delta_v = 4; 

%%  Initialize Satellite Parameters

    % Read Orbital Elements from file
    fileID = fopen('vehicleinfo_oe.txt','r');
    x = textscan(fileID,'%s',1,'delimiter','\n\r');
    header1 = cell2mat(x{1});
    y = textscan(fileID,'%s',1,'delimiter','\n\r');
    header2 = cell2mat(y{1});
    data = cell2mat(textscan(fileID,'%f %f %f %f %f %f %f','delimiter','/n/r'));
    fclose('all');

    numSats = data(end,1);
    oe = data(:,2:7); %pull oe from data
    oe(:,3:6) = oe(:,3:6).*pi/180; %convert appropriate values in oe to radians
    
     
    for index = 1:numSats
        [R_vehicles_ECI(index,1:3),V_vehicles_ECI(index,1:3)] = ...
            oe2eci(mu_earth,oe(index,1),oe(index,2),oe(index,3),oe(index,4),oe(index,5),oe(index,6)); 
        %convert OE for each vehicle to ECI
    end
    
    %solar days between jan 1, 2000 12:00 til jan 1 2018, 00:00 (wolfram alpha)
    clearvars fileID x y data header1 header2

%%  Calculate Satellite Orbits
    % Initialize arrays
    tvecLength = 1001;
    
    % need to have a delay of 0 if you want to "send it" at 0 seconds
 
    for index = 1:numSats
        
        % Pull Position and Velocity
        pos = R_vehicles_ECI(index,:)';
        vel = V_vehicles_ECI(index,:)';
        
        % set this up for 1 complete period of the orbit
        tvec = linspace(0,3*2*pi*sqrt(oe(index,1)^3/mu_earth),tvecLength);

        % Solve for satellite orbit
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        odefun = @(tout,yout) differinertial(tout,yout,mu_earth);
        [tout,yout] = ode113(odefun,tvec,[pos;vel],options);
        satTime(index,:) = tout(:);
        satState(index,:,:) = real(yout(:,:));

        
    end
        
%% Plot Globe, Sat 1 and possibe Trajectories
  
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';


F(3) = struct('cdata',[],'colormap',[]);

for j = 1:size_for_things
        
    theta = j/size_for_things*tvec(end);
    theta = -mod((280.4606 + 360.9856473*theta/86164)/180*pi,2*pi); %Radians

    figure('position',[100 60 1200 745])
    hold on
    load('topo.mat','topo','topomap1');
    [xearth,yearth,zearth,props,cax] = prepplot(size_for_things,theta);
    cax = newplot(cax);
    
    h(1) = surf(xearth,yearth,zearth,props,'parent',cax);
    axis equal
    xlabel(['X [km]'])
    ylabel(['Y [km]'])
    zlabel(['Z [km]'])
    view(127.5,30)
    % %%%%%%%

    cmapsize = 64;  % 64-elements is each colormap
    cvalue1 = [-7473 ,5731];
    C1 = min(cmapsize,round((cmapsize-1)*(props.Cdata-cvalue1(1))/(cvalue1(2)-cvalue1(1)))+1); 
    set(h(1),'CData',C1);

    colormap([topomap1;autumn(64)]);
    clearvars cmapsize C1 cax i topomap1 topo 
 
    h(2) = surf(xearth,yearth,zearth,zearth.*0+128,'faceAlpha',0,'edgecolor','none');
    for satNum = 1:numSats
        h(3) = plot3(satState(satNum,:,1),satState(satNum,:,2),satState(satNum,:,3),'m');
        h(4) = plot3(satState(satNum,round(tvecLength*j/size_for_things),1),...
            satState(satNum,round(tvecLength*j/size_for_things),2),satState(...
            satNum,round(tvecLength*j/size_for_things),3),'or','MarkerFaceColor','r');
    end
    drawnow
    F(j) = getframe(gcf,[390 220 440 340]);
    close

end
v = VideoWriter('Tester.avi');
open(v)
writeVideo(v,F)
close(v)
 


        
        
       
    
   
    
