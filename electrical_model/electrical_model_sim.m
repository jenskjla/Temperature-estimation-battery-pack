clc 
clearvars
%Sample time
T = 0.5;

%Initial parameters values
%Thermal model
Cc = 67;
Cs = 3.115;
Rc = 1.83;
Rs = 4.03;

%Electrical model
R = 0.05;
Re = 1.3; 
C1 = 500;  
C = 3000;
SOC_init = 0.1;

%Modify the ambient temperature as needed
Ta = 21;

%number of batteries
nBatt = 7; 

%Randomness added to each battery
[R1,R2,R3,R4,R5,R6,R7] = deal(R*0.85,R*1,R*1.15,R*0.9,R*0.85,R*1.1,R*1);
[R11,R12,R13,R14,R15,R16,R17] = deal(Re*1.2,Re*1.1,Re*0.8,Re*0.9,Re*1.0,Re*0.9,Re*1.2);
[C1,C2,C3,C4,C5,C6,C7] = deal(C1*1.1,C1*1.,C1*0.7,C1*0.9,C1*1.2,C1*1.,C1*1.1);

%Functions
Entrop   = @(x1) 0.001*(14.56547*x1.^5 - 35.3986*x1.^4 + 30.1126*x1.^3 ...
                      - 11.6516*x1.^2 + 2.7279*x1 - 0.3485);
soc_voc  = @(x) 0.824*x + 3.347;

vo_list = [];
for i = 1:nBatt
    n = num2str(i);
    v_list = ['v' n '_list'];
    lists.(v_list) = [0];
    
    Ts_list = ['Ts' n '_list'];
    lists.(Ts_list) = [Ta];

    Tc_list = ['Tc' n '_list'];
    lists.(Tc_list) = [Ta];

    Tm_list = ['Tm' n '_list'];
    lists.(Tm_list) = [Ta];

    soc_list = ['soc' n '_list'];
    lists.(soc_list) = [SOC_init];

    voc_list = ['voc' n '_list'];
    voc = soc_voc(lists.(soc_list)(1));
    lists.(voc_list) = [voc];

    vo_list = ['vo' n '_list'];
    lists.(vo_list) = [0];
    
    Q_list = ['Q' n '_list'];
    lists.(Q_list) = [0];
end

%While loop
n = 1;
stop_flag = false;  % Stop flag for while loop condition
x = [T];
i_battery = [0];
while ~stop_flag
    %One sample per iteration
    x(end+1) = x(n) + T;
    
    % Change this value to your desired current 
    I_nominal = 1; 
    i_battery(end+1) = I_nominal;
    
    %Current value of i_battery
    Ik = i_battery(n); 
   
    %For loop to simulate each battery
    for i = 1:nBatt
        k = num2str(i);
        if eval(['lists.soc' k '_list(n) < 0.95']) %Condition of battery almost fully charged 
            % Calculate v and soc
            eval(['v' k ' = exp(-T/(R1' k '*C' k ')) * lists.v' k '_list(n) + R1' k '* Ik* (1 - exp(-T/(R1' k '*C' k ')));']);
            eval(['SOC = lists.soc' k '_list(n) + (T/C)*Ik*0.98;']);
            
            %Function to calculate open circuit voltage with SOC
            voc = soc_voc(SOC);
          
            % Store SOC and voc
            eval(['lists.soc' k '_list(end+1) = SOC;']);
            eval(['lists.voc' k '_list(end+1) = voc;']);
            
            %Output voltage
            eval(['vo' k ' = voc + v' k ' + R' k '* Ik;']);
            %Heat transfer
            eval(['Q = (vo' k ' - voc) * Ik + Ik * (lists.Tm' k '_list(n)+273.15) * Entrop(SOC);']);
            eval(['lists.Q' k '_list(end+1) = Q;']);
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
        end
    end
    n = n + 1;
end

%Current graph respective to time
plot(x, i_battery);

x(end) = [];
i_battery(end) = [];


for i = 1:nBatt
    n = num2str(i);
    %Remove last item from V array b/c of calculus before if statement
    eval(['lists.v' n '_list(end) = [];']);
    eval(['figure(' n+1 ');']);
    eval(['plot(x,lists.Ts' n '_list);']);
    hold on
    eval(['plot(x,lists.Tc' n '_list);']);
    xlabel('samples n');                % X‑axis
    ylabel('Temperature T (°C)');       % Y‑axis
    title('Simulated values of temperature');
end


A = [x.', i_battery.', ...
    lists.vo1_list.', lists.Ts1_list.', lists.Tc1_list.', lists.soc1_list', lists.voc1_list', ...
    lists.vo2_list.', lists.Ts2_list.', lists.Tc2_list.', lists.soc2_list', lists.voc2_list', ...
    lists.vo3_list.', lists.Ts3_list.', lists.Tc3_list.', lists.soc3_list', lists.voc3_list', ...
    lists.vo4_list.', lists.Ts4_list.', lists.Tc4_list.', lists.soc4_list', lists.voc4_list', ...
    lists.vo5_list.', lists.Ts5_list.', lists.Tc5_list.', lists.soc5_list', lists.voc5_list', ...
    lists.vo6_list.', lists.Ts6_list.', lists.Tc6_list.', lists.soc6_list', lists.voc6_list', ...
    lists.vo7_list.', lists.Ts7_list.', lists.Tc7_list.', lists.soc7_list', lists.voc7_list'];
A = [["t", "ib", "vo1", "Ts1", "Tc1", "soc1", "ocv1", "vo2", "Ts2", "Tc2", "soc2", "ocv2",...
    "vo3", "Ts3", "Tc3", "soc3", "ocv3", "vo4", "Ts4", "Tc4", "soc4", "ocv4", ...
    "vo5", "Ts5", "Tc5", "soc5", "ocv5", "vo6", "Ts6", "Tc6", "soc6", "ocv6", ...
    "vo7", "Ts7", "Tc7", "soc7", "ocv7"]; A];


B = [x.', lists.Q1_list', lists.Q2_list', lists.Q3_list', lists.Q4_list', lists.Q5_list', lists.Q6_list', lists.Q7_list'];
B = [["t", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]; B];

%Output csv files
writematrix(A, "./Simulation_data/targetWave_elec.csv");
writematrix(B, "./Simulation_data/Qvalues_elec.csv");

