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
for timeLoop = 1:1:110
    % Initialize arrays
    tvecLength = 10001;
    
    % need to have a delay of 0 if you want to "send it" at 0 seconds
    t_simStart_delay = timeLoop*60; % variable delay until simulation start
    t_start_to_intercept = 25*60+t_simStart_delay; % 25 minute window to intercept
    t_deorbit_delay = [0]*60+t_simStart_delay; % [0 5 10 15 20] variable delays until launch
    
    satTime = zeros(numSats,tvecLength);
    satState = zeros(numSats,tvecLength,6);
    satTimeIndex = zeros(numSats,length(t_deorbit_delay));

    for index = 1:numSats
        
        % Pull Position and Velocity
        pos = R_vehicles_ECI(index,:)';
        vel = V_vehicles_ECI(index,:)';
        
        % set this up for 1 complete period of the orbit
        tvec = linspace(0,2*pi*sqrt(oe(index,1)^3/mu_earth),tvecLength);

        % Solve for satellite orbit
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        odefun = @(tout,yout) differinertial(tout,yout,mu_earth);
        [tout,yout] = ode113(odefun,tvec,[pos;vel],options);
        satTime(index,:) = tout(:);
        satState(index,:,:) = yout(:,:);

        timeindex = zeros(1,length(t_deorbit_delay));
        for i = 1:length(t_deorbit_delay)
            [~,timeindex(i)] = min(abs(t_deorbit_delay(i) - tout)); % get index for the needed time step
        end
        
        satTimeIndex(index,:) = timeindex(:);
        
    end
    
%% Find Satellite Intercept Possibilities

        % Define search_space
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

                    %iterate through different delays
                    for index = 1:length(t_deorbit_delay)

                        % Calculate short way lambert solution
                        [vSout,~,~] = lambert_v(mu_earth,squeeze(satState(satNum,satTimeIndex(satNum,index),1:3))',...
                            [x(i,j) y(i,j) z(i,j)],'s',0,t_start_to_intercept-t_deorbit_delay(index));
                        % Calculate long way lambert solution
                        [vLout,~,~] = lambert_v(mu_earth,squeeze(satState(satNum,satTimeIndex(satNum,index),1:3))',...
                            [x(i,j) y(i,j) z(i,j)],'l',0,t_start_to_intercept-t_deorbit_delay(index));
                        
                        % Keep lower velocity solution
                        %   Set a dummy variable, success2, as the current step's DV
                        %   needed to get to target position
                        minDV = min([norm(squeeze(satState(satNum,satTimeIndex(satNum,index),4:6))'-vSout'),...
                                    norm(squeeze(satState(satNum,satTimeIndex(satNum,index),4:6))'-vLout')]);
                        
                        % First delay time, need to initialize numbers in suc somewhere to be able to compare things to it
                        if index == 1   
                            
%                             vout = v2out; %set the DV for that iteration 
                            suc(i,j) = minDV;%update the suc matrix with the first set of DV's

                        % Rest of delay times, keep best solution available
                        elseif(suc(i,j)>minDV)

                            suc(i,j) = minDV;
%                             vout = v2out; %also output the v2out from lambert that gave a better DV
                            if (minDV <= max_delta_v)
                                disp('got better'); %if the improvement is within our limits, 
                                %print that the solver found a better solution
                                %after using a delay
                                %I have played with it a bit and have not got any
                                %increase without reversing the order of the 
                                %t_deorbit_delay variable
                            end
                            
                        end 
                        
                    end
                    
                    %update velocity vector to the best one
%                     v1_out(i,j,1:3) = vout; 
                    
                end
            end


            success = suc;
            % Keep only DV values under our threshold
            success(suc>max_delta_v) = NaN;
            viableTargets(satNum,:,:) = success;
            success_area(satNum) = sum(sum((success>-1).*...
                cos(asin(z/r))))/sum(sum(cos(asin(z/r))));
            

            %this value is what amount of the globe is covered by this

            %vehicle. 1 = total coverage, 0 = no coverage
            clearvars success2 suc i j index vout v2out 
            
        end
        
%% Plot Globe, Sat 1 and possibe Trajectories

%     for satNum = 1:numSats
%         
%         figure(satNum)
%         hold on
% 
%         load('topo.mat','topo','topomap1');
%         [xearth,yearth,zearth,props,cax] = prepplot(size_for_things);
%         cax = newplot(cax);
%         h(1) = surf(xearth,yearth,zearth,props,'parent',cax);
%         hold on
%         axis equal
%         xlabel(['X [km]'])
%         ylabel(['Y [km]'])
%         zlabel(['Z [km]'])
%         view(127.5,30)
%         % %%%%%%%
%         cmapsize = 64;  % 64-elements is each colormap
%         cvalue1 = [-7473 ,5731];
%         C1 = min(cmapsize,round((cmapsize-1)*(topo-cvalue1(1))/(cvalue1(2)-cvalue1(1)))+1); 
%         set(h(1),'CData',C1);
%         colormap([topomap1;autumn(64)]);
%         clearvars cmapsize C1 cax i j topomap1 topo
% 
%         %Plot sat 1 and its possible trajectory
%         h(2) = surf(x,y,z,64+squeeze(viableTargets(satNum,:,:))*(64/max_delta_v),'FaceAlpha',.9, 'EdgeColor','none');
%         h(3) = plot3(satState(satNum,:,1),satState(satNum,:,2),satState(satNum,:,3),'r');
%         h(4) = plot3(satState(satNum,1,1),satState(satNum,1,2),satState(satNum,1,3),'or','MarkerFaceColor','r');
% 
%         hold on
%         for i = 2:length(t_deorbit_delay)
%             h(4+i) = plot3(satState(satNum,timeindex(i),1),satState(satNum,timeindex(i),2),...
                %satState(satNum,timeindex(i),3),'om','MarkerFaceColor','m');
%         end
% 
%         disp(' ')
%         clearvars i tvec
%         
%         hold off
%         
%     end 
    
    boolHit = viableTargets;
    boolHit(~isnan(boolHit)) = 1;
    boolHit(isnan(boolHit)) = 0;
    targetOverlap = squeeze(sum(boolHit,1));
    targetOverlap(targetOverlap==0) = nan;
    
    figure(timeLoop)
    hold on

    load('topo.mat','topo','topomap1');
    [xearth,yearth,zearth,props,cax] = prepplot(size_for_things);
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
    C1 = min(cmapsize,round((cmapsize-1)*(topo-cvalue1(1))/(cvalue1(2)-cvalue1(1)))+1); 
    set(h(1),'CData',C1);
    colormap([topomap1;autumn(64)]);
    clearvars cmapsize C1 cax i j topomap1 topo

    %Plot sat 1 and its possible trajectory
    h(2) = surf(x,y,z,64+targetOverlap*(64/max(max(targetOverlap))),'FaceAlpha',.9, 'EdgeColor','none');

    hold on

    disp(' ')
    clearvars i tvec

    hold off
    
    M((timeLoop)) = getframe;
    
end