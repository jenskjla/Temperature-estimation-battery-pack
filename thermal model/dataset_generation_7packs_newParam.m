T = 0.1;
Ttotal = 300;
Ta = 35;

Rc = 1.83e0; % Ri
Cc = 6.7e1; % Ci
Cs = 3.115e0; % Co
Ru = 4.03e0; % Ro
Rcc = 4e-1;
Re = 0.1;
% R1 = 1.3;
% C1 = 500;
num = 7;

[Rc1,Rc2,Rc3,Rc4,Rc5,Rc6,Rc7] = deal(Rc);
[Re1,Re2,Re3,Re4,Re5,Re6,Re7] = deal(Re*0.85,Re*1,Re*1.15,Re*0.9,Re*0.85,Re*1.1,Re*1); %suppose Re of each battery is different
[Cs1,Cs2,Cs3,Cs4,Cs5,Cs6,Cs7] = deal(Cs);
[Cc1,Cc2,Cc3,Cc4,Cc5,Cc6,Cc7] = deal(Cc);
[Ru1,Ru2,Ru3,Ru4,Ru5,Ru6,Ru7] = deal(Ru);

Ts1_list = [21];
Ts2_list = [21];
Ts3_list = [21];
Ts4_list = [21];
Ts5_list = [21];
Ts6_list = [21];
Ts7_list = [21];
Tc1_list = [21];
Tc2_list = [21];
Tc3_list = [21];
Tc4_list = [21];
Tc5_list = [21];
Tc6_list = [21];
Tc7_list = [21];

i_seq = [];
Time = T:T:Ttotal;
[~, len] = size(Time);
modu = 2;
carrier = (sawtooth(2*pi*0.1/3*Time,1/2)+1)/2*3;
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
Time = Time(1:len);  

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
i_seq = i_battery;
%plot(Time, i_battery); 

listsPath = fullfile('..', 'electrical_model', 'Q_values.mat');
load(listsPath, 'lists');

for i = 1:num
   eval(['Q' num2str(i) ' = lists.Q' num2str(i) '_list';]);
end



A = [-1/(Rc1*Cc1), 1/(Rc1*Cc1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    1/(Rc1*Cs1), -1*(1/Cs1)*(1/Ru1+1/Rc1+1/Rcc), 0, 1/(Rcc*Cs1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, -1/(Rc2*Cc2), 1/(Rc2*Cc2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1/(Rcc*Cs2), 1/(Rc2*Cs2), -1/Cs2*(1/Ru2+1/Rc2+2/Rcc), 0, 1/(Rcc*Cs2), 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, -1/(Rc3*Cc3), 1/(Rc3*Cc3), 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1/(Rcc*Cs3), 1/(Rc3*Cs3), -1/Cs3*(1/Ru3+1/Rc3+2/Rcc), 0, 1/(Rcc*Cs3), 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, -1/(Rc4*Cc4), 1/(Rc4*Cc4), 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1/(Rcc*Cs4), 1/(Rc4*Cs4), -1/Cs4*(1/Ru4+1/Rc4+2/Rcc), 0, 1/(Rcc*Cs4), 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc5*Cc5), 1/(Rc5*Cc5), 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs5), 1/(Rc5*Cs5), -1/Cs5*(1/Ru5+1/Rc5+2/Rcc), 0, 1/(Rcc*Cs5), 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc6*Cc6), 1/(Rc6*Cc6), 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs6), 1/(Rc6*Cs6), -1/Cs6*(1/Ru6+1/Rc6+2/Rcc), 0, 1/(Rcc*Cs6);
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1/(Rc7*Cc7), 1/(Rc7*Cc7);
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/(Rcc*Cs7), 1/(Rc7*Cs7), -1/Cs7*(1/Ru7+1/Rc7+1/Rcc)]; 
B = [1/Cc1, 0;
     0, 1/(Cs1*Ru1);
     1/Cc2, 0;
     0, 1/(Cs2*Ru2);
     1/Cc3, 0;
     0, 1/(Cs3*Ru3);
     1/Cc4, 0;
     0, 1/(Cs4*Ru4);
     1/Cc5, 0;
     0, 1/(Cs5*Ru5);
     1/Cc6, 0;
     0, 1/(Cs6*Ru6);
     1/Cc7, 0;
     0, 1/(Cs7*Ru7)];

for i = 1:len
    x = [Tc1_list(i); Ts1_list(i);Tc2_list(i); Ts2_list(i);Tc3_list(i); Ts3_list(i);Tc4_list(i); Ts4_list(i);
        Tc5_list(i); Ts5_list(i);Tc6_list(i); Ts6_list(i);Tc7_list(i); Ts7_list(i)] 
    for n = 1:num;
        n = num2str(n);
        eval(['u = [Q' n '(i); Ta];'])
        rate = A * x + B * u;
        x_new = x + rate *T;
        eval(['Tc' n '_list = [Tc' n '_list, x_new(1)];'])
        eval(['Ts' n '_list = [Ts' n '_list, x_new(2)];'])
    end
end 
Ts1_list = Ts1_list(2:end);
Ts2_list = Ts2_list(2:end);
Ts3_list = Ts3_list(2:end);
Ts4_list = Ts4_list(2:end);
Ts5_list = Ts5_list(2:end);
Ts6_list = Ts6_list(2:end);
Ts7_list = Ts7_list(2:end);
Tc1_list = Tc1_list(2:end);
Tc2_list = Tc2_list(2:end);
Tc3_list = Tc3_list(2:end);
Tc4_list = Tc4_list(2:end);
Tc5_list = Tc5_list(2:end);
Tc6_list = Tc6_list(2:end);
Tc7_list = Tc7_list(2:end);
plot(Time, Ts1_list)
hold on 
plot(Time, Tc1_list)
title={'t','ib','Ts1','Ts2','Ts3','Ts4','Ts5','Ts6','Ts7','Tc1','Tc2','Tc3','Tc4','Tc5','Tc6','Tc7'};
result_table=table(Time', i_seq', Ts1_list', Ts2_list', Ts3_list', Ts4_list', Ts5_list', Ts6_list', Ts7_list', Tc1_list', Tc2_list', Tc3_list',Tc4_list', Tc5_list', Tc6_list',Tc7_list','VariableNames',title);
writetable(result_table, './Simulation_data/dataset_pack7_NewParam.csv');
title={'t','Q1','Q2','Q3','Q4','Q5','Q6','Q7'};
result_table = table(Time', Q1', Q2', Q3', Q4', Q5', Q6', Q7','VariableNames', title);
writetable(result_table, './Simulation_data/Q_values.csv');


