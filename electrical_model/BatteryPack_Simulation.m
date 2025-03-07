clear;
syms t s
rng(1);

T = 0.1;
x = T:T:300;
[~, len] = size(x);
modu = 2;
carrier = (sawtooth(2*pi*0.1/3*x,1/2)+1)/2*3;
A = [];
for i=1:len
    if carrier(i) <= modu
        A = [A, 0];
    else
        A = [A, 5];
    end
end

B1 = A(1:len-15/T);
B2 = A(15/T+1:end);
i_battery = B1 - B2;
[~, len] = size(i_battery);
x = x(1:len);
% plot(x, i_battery);  

start_index = -1;
end_index = -1;
for i = 51:len
    if i_battery(i)~=0 && i_battery(i-1)==0
        start_index = i;
    elseif i_battery(i)==0 && i_battery(i-1)~=0 && start_index ~= -1
        end_index = i-1;
        a = 0.5;
        b = 2;
        r = (b-a).*rand(1,1) + a;
        i_battery(start_index:end_index) = i_battery(start_index:end_index) * r;
        start_index = -1;
        end_index = -1;
    elseif i_battery(i)==0 && i_battery(i-1)~=0 && start_index == -1
        error("hha")
    end
end

Cc = 67;
Cs = 1*3.115;
Ri = 1.83;
Ro = 4.03;
R = 0.05;
Re1 = 1.3; 
C1 = 500;  
C = 3.0;
SOC_init = 0.9;
Ta = 35;
num = 7; %number of batteries
i_integration = 0;
[R1,R2,R3,R4,R5,R6,R7] = deal(R*0.85,R*1,R*1.15,R*0.9,R*0.85,R*1.1,R*1);
[R11,R12,R13,R14,R15,R16,R17] = deal(Re1*1.2,Re1*1.1,Re1*0.8,Re1*0.9,Re1*1.0,Re1*0.9,Re1*1.2);
[C1,C2,C3,C4,C5,C6,C7] = deal(C1*1.1,C1*1.,C1*0.7,C1*0.9,C1*1.2,C1*1.,C1*1.1);
[SOC1_init,SOC2_init,SOC3_init,SOC4_init,SOC5_init,SOC6_init,SOC7_init] = deal(SOC_init*0.9,SOC_init*0.7,SOC_init*1.,SOC_init*0.9,SOC_init*0.80,SOC_init*1.05,SOC_init*1);


vo_list = [];
for i = 1:num
    n = num2str(i);
    v_list = ['v' n '_list'];
    lists.(v_list) = [];

    Ts_list = ['Ts' n '_list'];
    lists.(Ts_list) = [Ta];

    Tc_list = ['Tc' n '_list'];
    lists.(Tc_list) = [Ta];

    soc_list = ['soc' n '_list'];
    eval(['randnum = SOC' n '_init;']);
%     randnum = 0.7 + (1.1-0.7) *rand();
    lists.(soc_list) = [randnum];

    voc_list = ['voc' n '_list'];
    voc = soc_voc(lists.(soc_list)(1));
    lists.(voc_list) = [voc];

    vo_list = ['vo' n '_list'];
    lists.(vo_list) = [];
    
    Q_list = ['Q' n '_list'];
    lists.(Q_list) = [];

    eval(['v' num2str(i) '_integration = 0;']);
    eval(['SOC' num2str(i) '_integration = 0;']);
    eval(['v' num2str(i) '_init = 0;']);
%     randnum = 0.85 + (1.2 - 0.85) * rand(); 
%     eval(['R' num2str(i) ' = R * randnum;']);
%     eval(['R1' num2str(i) ' = Re1 * randnum;']);
%     eval(['C' num2str(i) ' = C1 * randnum;']);
end


Tc2 = lists.Tc1_list;

for n=1:len
    i_integration = i_integration + i_battery(n)*T;
    for i=1:num
        k = num2str(i);
        eval(['v' k  ' = (i_integration - v' k  '_integration/R1' k ') / C' k ' + v' k '_init;']);
        eval(['v' k  '_integration = v' k '_integration + v' k '*T;']);
        eval(['SOC' k '_integration = SOC' k '_integration + i_battery(n)*T/(C*3600);']);
        eval(['SOC = SOC' k '_integration + lists.soc' k '_list(1);']);
        voc = soc_voc(SOC);
        eval(['lists.soc' k '_list(end+1) = SOC;']);
        eval(['lists.voc' k '_list(end+1) = voc;']);
        eval (['vo' k ' = voc + v' k ' + R' k '* i_battery(n);']);
        eval(['Q = (vo' k ' - voc) * i_battery(n) + i_battery(n) * lists.Tc' k '_list(end) * Entrop(SOC);']);
        eval(['lists.Q' k '_list(end+1) = Q;' ]);
        eval(['Tc_new = T/Cc * (Q - lists.Tc' k '_list(end)/Ri + lists.Ts' k '_list(end)/Ri + Cc/T*lists.Tc' k '_list(end));']);
        eval(['Ts_new = T/Cs * (lists.Tc' k '_list(end)/Ri +(Cs/T-1/Ri-1/Ro)*lists.Ts' k '_list(end) + Ta/Ro);']);
        eval(['lists.v' k '_list(end+1) = v' k ';'])
        eval(['lists.vo' k '_list(end+1) = vo' k ';']);
        eval(['lists.Tc' k '_list(end+1) = Tc_new;']);
        eval(['lists.Ts' k '_list(end+1) = Ts_new;']);
    end
    
    
    if i_battery(n)<0
        ll=1;
    end
    
end

for i = 1:num
    n = num2str(i);
    Tc_list = Tc_list(2:end);
    eval(['lists.Tc' n '_list = lists.Tc' n '_list(2:end);']);
    eval(['lists.Ts' n '_list = lists.Ts' n '_list(2:end);']);
    eval(['lists.soc' n '_list = lists.soc' n '_list(2:end);']);
    eval(['lists.voc' n '_list = lists.voc' n '_list(2:end);']);
    eval(['figure(' n ');']);
    eval(['plot(x,lists.Ts' n '_list);']);
    hold on
    eval(['plot(x,lists.Tc' n '_list);']);
end



% figure
% yyaxis left
% plot(x,vo_list)
% yyaxis right
% plot(x,i_battery)
% ylim([-8, 8]);



A = [x.', i_battery.', ...
    lists.vo1_list.', lists.Ts1_list.', lists.Tc1_list.', lists.soc1_list', lists.voc1_list', ...
    lists.vo2_list.', lists.Ts2_list.', lists.Tc2_list.', lists.soc2_list', lists.voc2_list', ...
    lists.vo3_list.', lists.Ts3_list.', lists.Tc3_list.', lists.soc3_list', lists.voc3_list', ...
    lists.vo4_list.', lists.Ts4_list.', lists.Tc4_list.', lists.soc4_list', lists.voc4_list', ...
    lists.vo5_list.', lists.Ts5_list.', lists.Tc5_list.', lists.soc5_list', lists.voc5_list', ...
    lists.vo6_list.', lists.Ts6_list.', lists.Tc6_list.', lists.soc6_list', lists.voc6_list', ...
    lists.vo7_list.', lists.Ts7_list.', lists.Tc7_list.', lists.soc7_list', lists.voc7_list'];
A = [["t", "ib", "vb1", "Ts1", "Tc1", "soc1", "voc1", "vb2", "Ts2", "Tc2", "soc2", "voc2",...
    "vb3", "Ts3", "Tc3", "soc3", "voc3", "vb4", "Ts4", "Tc4", "soc4", "voc4", ...
    "vb5", "Ts5", "Tc5", "soc5", "voc5", "vb6", "Ts6", "Tc6", "soc6", "voc6", ...
    "vb7", "Ts7", "Tc7", "soc7", "voc7"]; A];


writematrix(A, "./Simulation_data/targetWave.csv");
save('Q_values.mat', "lists")

function y = soc_voc(x)
    y = 0.824*x+3.347; 
end
function y1 = Entrop(x1)
    y1 = 14.56547 * x1.^5 -35.3986 * x1.^4 + 30.1126 * x1.^3 - 11.6516 * x1.^2 + 2.7279 * x1.^1 - 0.3485 * x1.^0; 
    y1 = y1 * 0.001;
end
