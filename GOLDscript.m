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
    time_sim = 45.4*86400; %seconds (days*86400) since epoch
    tlength = 1*3600; %seconds of simulation

    tvec = time_sim:1:time_sim+tlength;
    thetavec = mod((280.4606 + 360.9856473*tvec/86400)/180*pi,2*pi);
    size_for_things = 100;

    %target things
    target_alt = 1200; %(km)
    target_ratio = (R_earth+target_alt)/R_earth;

%% Initialize Satellite Parameters
    fileID = fopen('vehicleinfo_oe.txt','r');

    x = textscan(fileID,'%s',1,'delimiter','\n\r');
    header1 = cell2mat(x{1});
    y = textscan(fileID,'%s',1,'delimiter','\n\r');
    header2 = cell2mat(y{1});

    data = cell2mat(textscan(fileID,'%f %f %f %f %f %f %f',3,'delimiter','/n/r'));
    fclose('all');


    vehicles = data(end,1);
    oe = data(:,2:7); %pull oe from data
    oe(:,3:6) = oe(:,3:6).*pi/180; %convert appropriate values in oe to radians
    
     
    for index = 1:length(vehicles)
        [R_vehicles_ECI(index,1:3),V_vehicles_ECI(index,1:3)] = ...
            oe2eci(mu_earth,oe(index,1),oe(index,2),oe(index,3),oe(index,4),oe(index,5),oe(index,6)); 
        %convert OE for each vehicle to ECI

    end
    
    %solar days between jan 1, 2000 12:00 til jan 1 2018, 00:00 (wolfram alpha)
    days = 6574.5; %Days
    theta_epoch = mod((280.4606 + 360.9856473*days)/180*pi,2*pi); %Radians
    clearvars fileID x y data header1 header2
%% Initial Vehicle position


    r1 = R_vehicles_ECI(1,:)';
    v1 = V_vehicles_ECI(1,:)';
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    tvec = 0:.01:60*5;
    odefun = @(tout,yout) differinertial(tout,yout,mu_earth);
    [tout,yout] = ode113(odefun,tvec,[r1;v1],options);


%% finding sat 1 intercept possibilities

    % Define search_space
    r = target_ratio*R_earth; %radius of intercept

    [x,y,z] = sphere(size_for_things);
    x = x*r; %psuedo longitude variable
    y = y*r; %psuedo latitude var
    z = z*r; %psuedo radius var
    v1_out = x.*NaN; %initialize an array
    suc = v1_out;
    for i = 1:length(x)
        for j = 1:length(y)
            [vout, ~] = lambert_v(mu_earth,r1,[x(i,j) y(i,j) z(i,j)],'s',0,25*60);
            suc(i,j) = norm(v1-vout);
            [v2out,~] = lambert_v(mu_earth,yout(end,1:3),[x(i,j) y(i,j) z(i,j)],'s',0,20*60);
            success2 = norm(yout(end,4:6)-v2out);
            if suc(i,j)>success2
                suc(i,j) = success2;
            end

            v1_out(i,j,1:3) = vout;
        end
    end

    success = suc;
    success(suc>4) = NaN;
clearvars success2 suc
%% Plot Globe, Sat 1 and possibe Trajectories

    figure
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
    h(2) = surf(x,y,z,64+success*16,'FaceAlpha',.9, 'EdgeColor','none');
    h(3) = plot3(r1(1),r1(2),r1(3),'or','MarkerFaceColor','r');
    h(4) = plot3(yout(end,1),yout(end,2),yout(end,3),'og','MarkerFaceColor','g');
    disp(' ')
