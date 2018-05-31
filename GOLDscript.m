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
    size_for_things = 100;

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
    days = 6574.5; %Days
    theta_epoch = mod((280.4606 + 360.9856473*days)/180*pi,2*pi); %Radians
    clearvars fileID x y data header1 header2

%%  Calculate Satellite Orbits

    % Initialize arrays
    tvecLength = 1001;
    
    t_start_to_intercept = 25*60; % 25 minute window to intercept
    t_deorbit_delay = [0]*60; % [0 5 10 15 20] variable delays until launch
    

    satState = zeros(numSats,tvecLength,6);
    tvec = linspace(0,2*pi*sqrt(oe(1,1)^3/mu_earth),tvecLength); %one orbit
    
    for index = 1:numSats
        
        % Pull Position and Velocity
        pos = R_vehicles_ECI(index,:)';
        vel = V_vehicles_ECI(index,:)';
       
        % Solve for satellite orbit
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        odefun = @(tout,yout) differinertial(tout,yout,mu_earth);
        [tout,yout] = ode113(odefun,tvec,[pos;vel],options);
        satTime(:) = tout(:);
        satState(index,:,:) = yout(:,:);
        
    end
    clearvars pos vel
%% Find Satellite Intercept Possibilities
    
for timeLoop = 1:1:110
    index_of_deorbit = ceil(timeLoop/110*tvecLength);
    
    %Define search_space
    r = target_ratio*R_earth; %radius of intercept

    [x,y,z] = sphere(size_for_things);
    x = x*r; %psuedo longitude variable
    y = y*r; %psuedo latitude var
    z = z*r; %psuedo radius var
    viableTargets = zeros(numSats,length(x),length(y));

    % Loop through satellite firing solutions
    for satNum = 1:numSats

        suc = x.*NaN; %initialize an array

        %% The first sat

        %iterate through search space
        for i = 1:length(x)
            for j = 1:length(y)

                % Calculate short way lambert solution
                [vSout,~,~] = lambert_v(mu_earth,squeeze(satState(satNum,index_of_deorbit,1:3))',...
                    [x(i,j) y(i,j) z(i,j)],'s',0,t_start_to_intercept);
                % Keep lower velocity solution
                %   Set a dummy variable, success2, as the current step's DV
                %   needed to get to target position
                suc(i,j) = norm(squeeze(satState(satNum,index_of_deorbit,4:6))'-vSout');
                           
                %update the suc matrix with the first set of DV's

            end
        end


        success = suc;
        % Keep only DV values under our threshold
        success(suc>max_delta_v) = NaN;
        viableTargets(satNum,:,:) = success;

        %this value is what amount of the globe is covered by this

        %vehicle. 1 = total coverage, 0 = no coverage
        clearvars success2 suc i j index vout v2out 

    end
        
    
    boolHit = viableTargets;
    boolHit(~isnan(boolHit)) = 1;
    boolHit(isnan(boolHit)) = 0;
    targetOverlap = squeeze(sum(boolHit,1));
    targetOverlap(targetOverlap==0) = nan;
    
    figure
    hold on

    load('topo.mat','topo','topomap1');
    theta = timeLoop/110*tvec(end);
    theta = -mod((280.4606 + 360.9856473*theta/86164)/180*pi,2*pi); %Radians
    
    [xearth,yearth,zearth,props,cax] = prepplot(size_for_things,theta);
    cax = newplot(cax);
    h(1) = surf(xearth,yearth,zearth,props,'parent',cax);
    hold on
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
    %variable1 = jet(64);
    %colormap([topomap1;variable1(:,end:-1:1)]);
    colormap([topomap1;autumn(64)]);
    %clearvars cmapsize C1 cax i j topomap1 topo variable1

    %Plot sat 1 and its possible trajectory
    h(2) = surf(xearth,yearth,zearth,0.*zearth+2*64,'FaceAlpha',0,'EdgeColor','none');
    h(3) = surf(x,y,z,65+(targetOverlap-1)*63/5,'FaceAlpha',.9, 'EdgeColor','none');
    caxis([1 128]);
    hold on

    disp(' ')
    clearvars i

    hold off
    
    M((timeLoop)) = getframe;
    close
    
end