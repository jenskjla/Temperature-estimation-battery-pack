clear;
syms t s
rng(1);

%Initial parameters values
Cc = 67;
Cs = 3.115;
Rc = 1.83;
Rs = 4.03;

%Sample time
T = 0.5;

%Modify the ambient temperature as needed
Ta = 21;

%number of batteries
nBatt = 7; 

%Directory for estimated internal electrical parameters
baseDir = fullfile('..', 'electrical_model', 'outputGA');

% Preallocate
BatteryID    = (1:nBatt).';
R_estimated    = zeros(nBatt,1);
R1_estimated  = zeros(nBatt,1);
C1_estimated  = zeros(nBatt,1);
C_estimated    = zeros(nBatt,1);
SOCinit      = zeros(nBatt,1);

%Functions
Entrop   = @(x1) 0.001*(14.56547*x1.^5 - 35.3986*x1.^4 + 30.1126*x1.^3 ...
                      - 11.6516*x1.^2 + 2.7279*x1 - 0.3485);
soc_voc  = @(x) 0.824*x + 3.347;

% Pre‑allocate result vectors
R_estimated   = zeros(1, nBatt);
R1_estimated  = zeros(1, nBatt);
C1_estimated  = zeros(1, nBatt);
C_estimated   = zeros(1, nBatt);
SOCinit       = zeros(1, nBatt);

for b = 1:nBatt
    % build the full path to battery<b>_best.csv
    fname = fullfile(baseDir, sprintf('battery%d_best.csv', b));
    if ~isfile(fname)
        error('File not found: %s', fname);
    end
    
    Tb = readtable(fname);
    
    % Dynamic field‐names for this battery
    R_estimated(b)   = Tb.(sprintf('R%d',     b))(1);
    R1_estimated(b) = Tb.(sprintf('R1%d',    b))(1);
    C1_estimated(b) = Tb.(sprintf('C1%d',    b))(1);
    C_estimated(b)   = Tb.(sprintf('C%d',     b))(1);
    SOCinit(b)     = Tb.(sprintf('SOCinit%d',b))(1);
end

% Build a MATLAB table for easy indexing, plotting, and downstream use
bestLossTable = table( R_estimated, R1_estimated, ...
                       C1_estimated, C_estimated, ...
                       SOCinit, ...
    'VariableNames', { ...
      'R_estimated',   ... % Rb
      'R1_estimated', ... % R1b
      'C1_estimated', ... % C1b
      'C_estimated',   ... % Cb
      'SOCinit'        ... % SOCinitb
    });

% Display the assembled table
disp(bestLossTable);


for i = 1:nBatt
    n = num2str(i);    
    Ts_list = ['Ts' n '_list'];
    lists.(Ts_list) = [Ta];

    Tc_list = ['Tc' n '_list'];
    lists.(Tc_list) = [Ta];

    Tm_list = ['Tm' n '_list'];
    lists.(Tm_list) = [Ta];

    soc_list = ['soc' n '_list'];
    lists.(soc_list) = [bestLossTable.SOCinit(i)];
    
    Q_list = ['Q' n '_list'];
    lists.(Q_list) = [0];

    voc_list = ['voc' n '_list'];
    voc = soc_voc(lists.(soc_list)(1));
    lists.(voc_list) = [voc];

    vo_list = ['vo' n '_list'];
    lists.(vo_list) = [0];

    v_list = ['v' n '_list'];
    lists.(v_list) = [0];
end

%While loop parameters
n = 1;
stop_flag = false;  % Stop flag for while loop condition
x = [T];
i_battery = [0];
while ~stop_flag  
    %One sample per while cycle
    x(end+1) = x(n) + T;
    
    % Change this value to your desired current 
    I_nominal = 1; 
    
    % Generate current with ±0.2A error around nominal value
    %I_nominal = (I_nominal - 0.2) + (0.4 * rand()); 
    i_battery(end+1) = I_nominal;
    
    %Current value of i_battery
    Ik = i_battery(n); 
    %For loop to simulate each battery
    for i = 1:nBatt
        k = num2str(i);
        R = bestLossTable.R_estimated(i)
        R1 = bestLossTable.R1_estimated(i)
        C = bestLossTable.C_estimated(i)
        C1 = bestLossTable.C1_estimated(i)
        if eval(['lists.soc' k '_list(n) < 0.95']) %Condition of battery almost fully charged
            %Calculate v and soc
            eval(['v' k ' = exp(-T/(R1*C1)) * lists.v' k '_list(n) + R1* Ik* (1 - exp(-T/(R1*C1)));']);
            eval(['SOC = lists.soc' k '_list(n) + (T/C)*Ik*0.98;']);

            %Function to calculate open circuit voltage with SOC
            voc = soc_voc(SOC);
            
            % Store SOC and voc
            eval(['lists.soc' k '_list(end+1) = SOC;']);
            eval(['lists.voc' k '_list(end+1) = voc;']);
            
            %Output voltage
            eval(['vo' k ' = voc + v' k ' + R* Ik;']);
            %Heat transfer
            eval(['Q = (vo' k ' - voc) * Ik + Ik * (lists.Tm' k '_list(n)+273.15) * Entrop(SOC);']);
            eval(['lists.Q' k '_list(n) = Q;']);
            %Core and surface temperature
            eval(['Tc_new = (lists.Tc' k '_list(end)+273.15) + T*Q/Cc-T/(Cc*Rc)*[lists.Tc' k '_list(end)-lists.Ts' k '_list(end)];']);
            eval(['Ts_new = (lists.Ts' k '_list(end)+273.15) + T/(Cs*Rc)*[lists.Tc' k '_list(end)-lists.Ts' k '_list(end)]-T/(Cs*Rs)*[lists.Ts' k '_list(end)-Ta];']);
            %Medium temperature of battery
            Tm_new = (Tc_new+Ts_new)/2;
            eval(['lists.v' k '_list(end+1) = v' k ';']);
            eval(['lists.vo' k '_list(end+1) = vo' k ';']);
            eval(['lists.Tc' k '_list(end+1) = Tc_new - 273.15;']);
            eval(['lists.Ts' k '_list(end+1) = Ts_new - 273.15;']);
            eval(['lists.Tm' k '_list(end+1) = Tm_new - 273.15;']);
        else
            stop_flag = true;
            n = n+1;
            break
        end
    end

    if stop_flag
        break        % now break out of the while-loop
    end

    n = n+1;
end

%Current graph respective to time
plot(x, i_battery);

%Index modifications needed due to inconsistency for simulation
i_battery(end) = []; 
x(end) = [];

for i = nBatt:-1:1
    k = num2str(i);
    kplus  = num2str( str2double(k) - 1 );
    if str2num(kplus)>0
    % Get the current list previous lengths
    Ts_list_prev = eval(['lists.Ts' k '_list']);
    Tc_list_prev = eval(['lists.Tc' k '_list']);

    Q_list_prev = eval(['lists.Q' k '_list']);

    % Get the current list lengths 
    Ts_list = eval(['lists.Ts' kplus '_list']);
    Tc_list = eval(['lists.Tc' kplus '_list']);

    Q_list = eval(['lists.Q' kplus '_list']);
    
    % Check index inconsistency
    if length(Ts_list_prev) < length(Ts_list)
        eval(['lists.Ts' kplus '_list(end) = [];']);
    end
    if length(Tc_list_prev) < length(Tc_list)
        eval(['lists.Tc' kplus '_list(end) = [];']);
    end
    % Check index inconsistency
    if length(Q_list_prev) < length(Q_list)
        eval(['lists.Q' kplus '_list(end) = [];']);
    end
    end
end


for i = 1:nBatt
    n = num2str(i);
    eval(['figure(' n+1 ');']);
    eval(['plot(x,lists.Ts' n '_list);']);
    hold on
    eval(['plot(x,lists.Tc' n '_list);']);
end


A = [x.', i_battery.', ...
    lists.Ts1_list.', lists.Tc1_list.',...
    lists.Ts2_list.', lists.Tc2_list.',...
    lists.Ts3_list.', lists.Tc3_list.',...
    lists.Ts4_list.', lists.Tc4_list.',...
    lists.Ts5_list.', lists.Tc5_list.',...
    lists.Ts6_list.', lists.Tc6_list.',...
    lists.Ts7_list.', lists.Tc7_list.'];
A = [["t", "ib","Ts1", "Tc1","Ts2", "Tc2",...
    "Ts3", "Tc3","Ts4", "Tc4", ...
    "Ts5", "Tc5","Ts6", "Tc6", ...
    "Ts7", "Tc7",]; A];

x(end)=[];
B = [x.', lists.Q1_list', lists.Q2_list', lists.Q3_list', lists.Q4_list', lists.Q5_list', lists.Q6_list', lists.Q7_list'];
B = [["t", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]; B];

%Output csv files
writematrix(A, "./Simulation_data/targetWave.csv");
writematrix(B, "./Simulation_data/Qvalues.csv");

