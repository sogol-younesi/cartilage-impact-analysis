% Impact Parameters Calculator
% Developed by:
% Hessam Noori
% Version 1.4 - January 2022
% New Features of this Version:
%     - Using derivative of force to find out start and end point of impact
%     - Setting zero deformation to minimum deformation point in the impact
%     - Calculating work for loading cycle only
%     - Scaling graphs based on the parameters defined for each variable
%     - Graph parameters for load and acceleration are controlled automatically
%     - Pressure plot was added
%%

clc;
clear;
close all;

%% Input Parameters

accelscale = (1000)*(1/1.023)*(9.807); % (1000mV / 1V) * (1g / 1.023mV) * (9.807 m/s^2 / 1g)
loadscale = (1000)*(1/0.225); % (1000mV / 1V) * (1N / 0.225mV)   %(In the old calculations, loadscale factor was 0.219 mv/N)

%% Importing Sensor Data

files = dir('*.txt');
[nfiles,~] = size(files);

% Input Data:
% mass = 1412.56;             % Impactor Mass FOR RABBITS STUDY (g)
mass = 358.39;              % Impactor Mass FOR TISSUE STUDY (g)
Diameter = 8.5;             % Sample Diameter (mm)
Thickness = 0.563;            % Sample Thickness (mm), mesured by sample thickness mapping using indentation machine
Imp_diam = 8.5;               % Impactor Diamater, For Calclating Max Stress (mm)

A_impactor = 0.25*pi*(Imp_diam*10^-3)^2;                       % Area, m^2

fprintf('Name,Peak Load (N),Filtered (N),Duration (s),Impulse (N.s),Work (J),Kinetic Energy (J),Max Stress (MPa),Filtered (MPa),Loading Rate (MPa/ms),Filtered (MPa),Displacement (mm),Gravity (m/s^2)\n');

for n=1:nfiles
    
    header_line = 30;                    % number of header lines
    filename = files(n).name;
    datastruct = importdata(filename,'\t',header_line);
    
    % Naming for Output Data and Graphs
    dot = strfind(filename,'.txt');
    Name = filename(1:dot-1);
    
    % Data Analysis
    D = datastruct.data;
    
    t = D(:,1);                           % raw time
    a = D(:,2)*accelscale;                % raw acceleration
    f = D(:,3)*loadscale;                 % raw force
    
    [max_f,i_max_f] = max(f);
    window = 500;                         % The code works well for window times between 2000 and 5000
    
    time = t(i_max_f-window:i_max_f+window);
    A_raw = a(i_max_f-window:i_max_f+window);
    F_raw = f(i_max_f-window:i_max_f+window);
        
    T = (time-time(1))*1000;              % Set t=0 in the begining of the window and changes time unit to mili second;
    
    %% Integrating Acceleration
    
    v = cumtrapz(t,a);                    % Raw velocity, integration was made from t=0 to the end of the timespan
    x = cumtrapz(t,v);                    % raw position, integration was made from t=0 to the end of the timespan
    
    V_raw = v(i_max_f-window:i_max_f+window);
    X_raw_0 = x(i_max_f-window:i_max_f+window);

    %% Filtering and Smoothing Data
    
    % Max of filtered load
    F_prime = smooth(F_raw,0.03,'moving');
    F = smooth(F_raw,0.007,'moving');
    [max_Force,i_max_Force] = max(F);
    
    [n_F,~] = size(F);
    dFdt_raw = zeros(n_F-2,1);
    for i=1:n_F-2
        dFdt_raw(i) = (-3* F_prime(i)+4* F_prime(i+1)- F_prime(i+2))/(2*(T(i+1)-T(i)));
    end
    dFdt = smooth(dFdt_raw,0.015,'moving');
    
    [f1,f2] = butter(3,0.08,'low');
    A = filter(f1,f2,A_raw);
    max_Accel = max(A);
    
    V = smooth(V_raw,0.025,'moving');
    
    X0 = smooth(X_raw_0,0.025,'moving');
    
    %% Impact Start and End Points for Force
    
    % Impact start and end index claculated when derivative of force
    % reaches (p*max(force)) before peak load for start index, or drops to
    % p*max(force) after peak load for end index
    % Raw force data used to calculate max force
    [max_Force_raw,i_max_Force_raw] = max(F_raw);
    fp0 = 0.05;                              % Percentage of max Force derivative to identify the point as impact start and end
    b0 = 100;                                % Dismiss first and last "b0" data in Force and dFdt to remove the effect of noise
    max_dFdt = max(dFdt);
    [min_dFdt,i_min_dFdt] = min(dFdt(i_max_Force_raw:n_F-b0));
    i_min_dFdt = i_min_dFdt+i_max_Force_raw;

    % Impact start index for force using F_prime
    i_start = find(dFdt(b0:i_max_Force_raw)>=max_dFdt*fp0,1)+b0;
    
    % Impact end index for force using F_prime
    i_end = i_min_dFdt+find(dFdt(i_min_dFdt:n_F-b0)>=min_dFdt*fp0,1);
    
    duration = T(i_end)-T(i_start);
    
    gravity = mean(A(1:i_start));
    
    %% Setting starting Point of X to Zero
    
    [min_X0,i_min_X00] = min(X0(i_start:i_end));       % Finding minimum between impact start and end 
    i_min_X0 = i_min_X00+i_start-1;                    % Finding the number of minimum deformation point in the X matrix
    X = (X0-X0(i_min_X0))*1000;                        % Set X=0 to the minimum point of X0 and changing unit of X to milimeter;
    X_raw = (X_raw_0-X_raw_0(i_min_X0))*1000;          % Set X_raw=0 to the starting point of impact, changes the sign and changes X unit to mili meter;
    
    max_def = X(i_start)-X(i_min_X0);
    
    %% Work, Impulse, Max Stress, Impact Energy (Density) and Loading Rate Calculations
    
    % Impulse
    IMP = trapz((T(i_start:i_end)-T(i_start))/1000,F(i_start:i_end));
    
    % Work
    WRK = -trapz(X(i_start:i_max_Force_raw)/1000,F(i_start:i_max_Force_raw));
    
    % Max Stress
    P_raw = F_raw/A_impactor*10^-6;         % Raw pressure, MPa
    [max_P_raw,i_max_P_raw] = max(P_raw);
    
    P = F/A_impactor*10^-6;                  % Smoothed pressure, MPa
    [max_P,i_max_P] = max(P);
    
    %Loading Rate
    LRT = (P(i_max_P)-P(i_start))/(T(i_max_P)-T(i_start));                       % Smooth loading rate, MPa/ms
    LRT_raw = (P_raw(i_max_P_raw)-P_raw(i_start))/(T(i_max_P_raw)-T(i_start));       %Raw loading rate, MPa/ms

    % Impact Energy and Impact Energy Density
    V_imp = V_raw(i_start);
    volume = pi*Diameter^2*Thickness/4*10^-9;                        % Volume, m^3
    
    E_impact = (mass*1e-3)*V_imp^2/2;                                          % Impact Energy, J
    E_impact_density = mass*(10^-3)*V_imp^2/(2*volume)*10^(-6);                % Impact Energy Density, mJ/mm^3 
    
    %% Printing
    fprintf('==============================================================\n');
    fprintf('    \n\n****    Gravity             =     %f    ****\n\n',gravity);
    fprintf('Data for %s (Impactor Mass =  %5.2f, Diameter = %2.1f):\n',Name,mass,Imp_diam);
    fprintf('    Max Load (Filtered)            =     %f (%f) N\n',max_Force_raw,max_Force);
    fprintf('    Max Acceleration               =     %f  m/s^2\n',max_Accel);
    fprintf('    Impact Duration                =     %f  ms\n',duration);
    fprintf('    Impulse                        =     %f  N.s\n',IMP);
    fprintf('    Work                           =     %f  J\n',WRK);
    fprintf('    Impact Energy                  =     %f  J\n',E_impact);
    fprintf('    Max Stress (Filtered)          =     %f  (%f) MPa\n',max_P_raw,max_P);
    fprintf('    Loading Rate (Filtered)        =     %f  (%f) MPa/ms\n',LRT_raw,LRT);
    fprintf('    Impact Energy Density          =     %f  mJ/mm^3\n',E_impact_density);
    fprintf('    Max Deformation                =     %f  mm\n\n',max_def);
    
    % Printing for excel file
    fprintf('%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f', Name,max_Force_raw,max_Force,duration,IMP,WRK,E_impact,max_P_raw,max_P,LRT_raw,LRT,max_def,gravity);
    fprintf('\n\n');
    
    %% Plotting
    
    f_par = ceil(max_Force_raw/500)*500;
    p_par = ceil(max_P_raw/120)*120;
    t_par = 40;
    a_par = ceil(max_Accel/500)*500;
    d_par = ceil(max_def/5)*5;
    v_par = ceil(max(V)/2)*2;
    h_par = 4;
    
    % Force Plot
    figure
    hold on
    xlabel('Time (ms)');
    ylabel('Force (N)');
    title(['Force Plot for ', Name]);
    plot(T, F_raw, 'color', 'b'); % Raw force data in blue
    plot(T, F_prime, 'color', 'g'); % Smoothed force data (F_prime) in green
    plot(T, F, 'color', 'r'); % Smoothed force data (F) in red
    plot(T(i_start), F_prime(i_start), 'x', 'color', 'black', 'linewidth', 2); % Marking start point on F_prime
    plot(T(i_end), F_prime(i_end), 'x', 'color', 'black', 'linewidth', 2); % Marking end point on F_prime
    hold off
    
    % Pressure Plot
    figure
    hold on
    xlabel('Time (ms)');
    ylabel('Pressure (MPa)');
    title(['Pressure Plot for ', Name]);
    plot(T, P_raw, 'color', 'b'); % Raw pressure data in blue
    plot(T, P, 'color', 'r'); % Smoothed pressure data in red
    plot(T(i_start), P(i_start), 'x', 'color', 'black', 'linewidth', 2); % Marking start point on P
    plot(T(i_end), P(i_end), 'x', 'color', 'black', 'linewidth', 2); % Marking end point on P
    hold off
    
    % Acceleration Plot
    figure
    hold on
    xlabel('Time (ms)');
    ylabel('Acceleration (m/s^2)');
    title(['Acceleration Plot for ', Name]);
    plot(T, A_raw, 'color', 'b'); % Raw acceleration data in blue
    plot(T, A, '--', 'color', 'r', 'linewidth', 2); % Smoothed acceleration data in red
    plot(T(i_start), A(i_start), 'x', 'color', 'black', 'linewidth', 2); % Marking start point on A
    plot(T(i_end), A(i_end), 'x', 'color', 'black', 'linewidth', 2); % Marking end point on A
    plot(T(i_min_X0), A(i_min_X0), 'o', 'color', 'black', 'linewidth', 2); % Marking minimum deformation point on A
    hold off
        
    % Velocity Plot
    figure
    hold on
    xlabel('Time (ms)');
    ylabel('Velocity (m/s)');
%    xlim([0 t_par])
%    ylim([-v_par v_par])
    title(['Velocity Plot for ',Name]);
    plot(T,V_raw,'color','b');
    plot(T,V,'--','color','r','linewidth',2);
    plot(T(i_start),V(i_start),'x','color','black','linewidth',2);
    plot(T(i_end),V(i_end),'x','color','black','linewidth',2);
    plot(T(i_min_X0),V(i_min_X0),'o','color','black','linewidth',2);
    hold off
    
    % Deformation Plot
    figure
    hold on
    xlabel('Time (ms)');
    ylabel('Deformation (mm)');
%    xlim([0 t_par])
%    ylim([0 d_par])
    title(['Deformation Plot for ',Name]);
    plot(T,X_raw,'color','b');
    plot(T,X,'--','color','r','linewidth',2);

plot(T(i_start),X(i_start),'x','color','black','linewidth',2);
    plot(T(i_end),X(i_end),'x','color','black','linewidth',2);
    plot(T(i_min_X0),X(i_min_X0),'o','color','black','linewidth',2);
    hold off

    % Hysteresis Curve
    figure
    hold on
    xlabel('Deformation (mm)');
    ylabel('Force (N)');
%    xlim([0 d_par])
%    ylim([0 f_par])
    title(['Hysteresis Curve for ',Name]);
    plot(X(i_start:i_end),F(i_start:i_end),'color','b');
    plot(X(i_start),F(i_start),'*');
    plot(X(i_max_Force_raw),F(i_max_Force_raw),'*');
    hold off
    
end