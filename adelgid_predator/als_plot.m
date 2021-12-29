%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A. tsugae, L. nigrinus, S. tsugae model
% Plot of baseline results- parameters from fixed values
% and parameter estimation results
% Change initial condition to test species scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Baseline- all three species present
% a_0 = 1.46016237638577000; 
% la_0 = 3;
% sa_0 = 0.50009720527582800;

% No L. nigrinus
a_0 = 1.46016237638577000; 
la_0 = 0;
sa_0 = 0.50009720527582800;

% No S. tsugae
% a_0 = 1.46016237638577000; 
% la_0 = 3;
% sa_0 = 0;

% No predators
% a_0 = 1.46016237638577000; 
% la_0 = 0;
% sa_0 = 0;

init_1 = [a_0; la_0; sa_0]; 

% set parameter values
r_a = 6.56761352407444000; 
r_l = 14.99298004937450000;
r_s = 2.22950054246204000;

c_a = 0.04035305800579790; 
c_s = 0.25;

m_a = .0613;
m_aw = .0739; 
m_as = .0010; 
m_s = .0803; 

m_ae = 0;

m_ll = .087;
m_la = .054; 

m_sl = 0.14;
m_sa = 0.051;

f_la = 0.02236604733325600;
f_ll = 0.01952742799642410;
f_sl = 0.00121152119252237;
f_sa = 0.00222354285686278;

g_la = 0.00100429012840320;
g_sl = 0.00112623405252762;
g_sa = 0.00101295697925473;

p_l = 0.89998146858627400; 
p_s = 0.844; 

d_l = 0.75726380107577700;
d_a = 1.25927143640118000; 

b = 64.55902089394340000;
h = 0.10006477793840900; 

Hset = 0.94992632755891300;

%create a vector of the parameters (28 entries)
       %1    2    3    4    5    6     7    8     9     10   11    12    13    14    15    16    17    18    19    20    21    22   23   24   25   26  27  28
pars = [r_a  r_l  r_s  c_a  c_s  m_ae  m_a  m_aw  m_as  m_s  m_ll  m_la  m_sl  m_sa  f_ll  f_la  f_sl  f_sa  g_la  g_sl  g_sa  p_l  p_s  d_l  d_a  b   h   Hset];

% time
% total time in years
end_t_year = 3;

% time point (2 time points per week) at end of each system
% used in saving solutions in model_sol vector
end_t_1 = 12; 
end_t_2 = 32;
end_t_3 = 34;
end_t_4 = 38;
end_t_5 = 46;
end_t_6 = 48;
end_t_7 = 52;
end_t_8 = 56;
end_t_9 = 58;
end_t_10 = 60;
end_t_11 = 62;
end_t_12 = 72;
end_t_13 = 82;
end_t_14 = 98;

% time points (2/week) in each system ((number of weeks)*2+1 entries in each vector)
t_1 = linspace(0, 6, 13);
t_2 = linspace(0, 10, 21);
t_3 = linspace(0, 1, 3);
t_4 = linspace(0, 2, 5);
t_5 = linspace(0, 4, 9);
t_6 = linspace(0, 1, 3);
t_7 = linspace(0, 2, 5);
t_8 = linspace(0, 2, 5);
t_9 = linspace(0, 1, 3);
t_10 = linspace(0, 1, 3);
t_11 = linspace(0, 1, 3);
t_12 = linspace(0, 5, 11);
t_13 = linspace(0, 5, 11);
t_14 = linspace(0, 8, 17);
t_15 = linspace(0, 3, 7);

% total time points in full run of model
t_full = linspace(0, end_t_year * 52, end_t_year * 104 + 1);

% initial conditions
a_0 = x(9); 
la_0 = 3; %from data
sa_0 = x(8);
% create a vector of initial conditions
init_1 = [a_0; la_0; sa_0]; 

% saving space for whole run of state values
model_sol = zeros(104 * end_t_year + 1, 6);

% data
% in time vectors, entries in vector on RHS are in weeks (week numbers with 
% data), entries in vector on LHS are time points in model
% we have data for 5 out of 6 classes in model
ll_data_time = 2 * [21	23	24	73	127	128	129	130	131	132];
ll_data = [3.666666667 1.75 1.666666667 4	3 4.75 7.8 4.6 3 1.5];

la_data_time = 2 * [1	3	3.5	11	18	19	20	21	23	49	51	52	55	56	57	58	62	72	73	74	75	103	104	105	106	107	108	109	111	112	113	114	115	117	118	119	120	121	123	124	125	127	129	130	131	132];
la_data = [1.25	2	1.5	1	1	2.5	1	1.333333333	1.333333333	1	1.4	1	1	1	1	1	1	2.6	2.5	2	1	3.428571429	4.714285714	5.111111111	3.1	3.571428571	4.333333333	4	3.5	3.25	1	1	3.333333333	1	1.833333333	1.714285714	5	2	1	1	2	2	2	1	1	1];

sl_data_time = 2 * [78	79	80	82	84	134	135	136	138	139	140	144];
sl_data = [1	1	1	1	1	2	4.333333333	4	1.5	1.5	2.25	1];

sa_data_time = 2 * [19	21	32	33	35	57	80	83	84	85	88	89	90	94	97	104	106	125	127	130	133	135	136	138	139	140	141	142	143	144];
sa_data = [1	1	1	1	3	1	1	2	2	1	1.5	1.5	1	1	1	1	1	1	2	1	1.25	1	1	1.75	2.6	1.5	1.666666667	1	1	1];

a_data_time = 2 * [14	50];
a_data = [0.519444444	4.4125];

% to solve odes, run through through loop for each year
for i = 1:end_t_year
    
    %save space for individual system solutions
            
    %1 %october, november, december (6 weeks)
    %1 %A, LN adults, ST adults
    y_1 = zeros(length(t_1),3); 
    
    %2 %december, january, february (10 weeks)
    %2 %A, LN adults
    y_2 = zeros(length(t_2),2); 
    
    %3 %february (1 week)
    %3 %A eggs, A, LN adults
    y_3 = zeros(length(t_3),3); 
    
    %4 %february, march (2 weeks)
    %4 %A eggs, A, LN adults
    y_4 = zeros(length(t_4),3); 
    
    %5 %march (4 weeks)
    %5 %A eggs, A, LN larvae, LN adults, ST adults
    y_5 = zeros(length(t_5),5); 
    
    %6 %april (1 week)
    %6 %A eggs, A, LN larvae, LN adults, ST adults
    y_6 = zeros(length(t_6),5); 
    
    %7 %april (2 weeks)
    %7 %A, LN larvae, LN adults, ST larvae, ST adults
    y_7 = zeros(length(t_7),5); 
    
    %8 %april, may (2 weeks)
    %8 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    y_8 = zeros(length(t_8),6); 
    
    %9 %may (1 week)
    %9 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    y_9 = zeros(length(t_9),6); 
    
    %10 %may (1 week)
    %10 %A eggs, A, ST larvae, ST adults
    y_10 = zeros(length(t_10),4); 
    
    %11 %may (1 week)
    %11 %A eggs, A, ST larvae, ST adults
    y_11 = zeros(length(t_11),4); 
    
    %12 %may, june (5 weeks)
    %12 %A, ST larvae, ST adults
    y_12 = zeros(length(t_12),3); 
    
    %13 %july, august (5 weeks)
    %13 %A, ST larvae, ST adults
    y_13 = zeros(length(t_13),3); 
    
    %14 %august, september (8 weeks)
    %14 %A, ST adults
    y_14 = zeros(length(t_14),2); 
    
    %15 %september, october (3 weeks)
    %15 %A, LN adults, ST adults
    y_15 = zeros(length(t_15),3);
      
    
%1 %october, november, december (6 weeks)
%1 %A, LN adults, ST adults
    % solving the first system of odes 
    [s_1,y_1] = ode45(@(s_1,y_1)hwaode_ALS_15system_1(s_1,y_1,pars),t_1,init_1);
    
    % saving initial conditions for next system
            %2 %A, LN adults
    init_2 = [y_1(length(t_1),1);y_1(length(t_1),2)];
    
    % saving results in matrix with full run of all state variables
    for j = 1:(length(t_1)-1)
                                %1 %A, LN adults, ST adults
        model_sol((i-1)*104+j,:) = [0;y_1(j,1);0;y_1(j,2);0;y_1(j,3)];
    end
   
%2 %december, january, february (10 weeks)
%2 %A, LN adults
    [s_2,y_2] = ode45(@(s_2,y_2)hwaode_ALS_15system_2(s_2,y_2,pars),t_2,init_2);
    
            %3 %A eggs, A, LN adults
    init_3 = [0;y_2(length(t_2),1);y_2(length(t_2),2)];
    
    for j = 1:(length(t_2)-1)
                                        %2 %A, LN adults
        model_sol((i-1)*104+j+end_t_1,:) = [0;y_2(j,1);0;y_2(j,2);0;0];
    end
    
%3 %february (1 week)
%3 %A eggs, A, LN adults
    %system 3
    [s_3,y_3] = ode45(@(s_3,y_3)hwaode_ALS_15system_3(s_3,y_3,pars),t_3,init_3);
    
            %4 %A eggs, A, LN adults
    init_4 = [y_3(length(t_3),1);y_3(length(t_3),2);y_3(length(t_3),3)];
    
    for j = 1:(length(t_3)-1)
                                        %3 %A eggs, A, LN adults
        model_sol((i-1)*104+j+end_t_2,:) = [y_3(j,1);y_3(j,2);0;y_3(j,3);0;0];
    end
    
%4 %february, march (2 weeks)
%4 %A eggs, A, LN adults
    %system 4
    [s_4,y_4] = ode45(@(s_4,y_4)hwaode_ALS_15system_4(s_4,y_4,pars),t_4,init_4);
    
            %LN "reproduction"
            %ST overwintering adults become active again
            %5 %A eggs, A, LN larvae, LN adults, ST adults
    init_5 = [y_4(length(t_4),1);y_4(length(t_4),2);pars(2)*y_3(length(t_3),3);y_4(length(t_4),3);pars(23)*y_1(length(t_1),3)];
    
    for j = 1:(length(t_4)-1)
                                        %4 %A eggs, A, LN adults
        model_sol((i-1)*104+j+end_t_3,:) = [y_4(j,1);y_4(j,2);0;y_4(j,3);0;0];
    end
    
%5 %march (4 weeks)
%5 %A eggs, A, LN larvae, LN adults, ST adults
    %system 5
    [s_5,y_5] = ode45(@(s_5,y_5)hwaode_ALS_15system_5(s_5,y_5,pars),t_5,init_5);
    
            %6 %A eggs, A, LN larvae, LN adults, ST adults
    init_6 = [y_5(length(t_5),1);y_5(length(t_5),2);y_5(length(t_5),3);y_5(length(t_5),4);y_5(length(t_5),5)];
    
    for j = 1:(length(t_5)-1)
                                        %5 %A eggs, A, LN larvae, LN adults, ST adults
        model_sol((i-1)*104+j+end_t_4,:) = [y_5(j,1);y_5(j,2);y_5(j,3);y_5(j,4);0;y_5(j,5)];
    end
    
%6 %april (1 week)
%6 %A eggs, A, LN larvae, LN adults, ST adults
    %system 6
    [s_6,y_6] = ode45(@(s_6,y_6)hwaode_ALS_15system_6(s_6,y_6,pars),t_6,init_6);
    
            %7 %A, LN larvae, LN adults, ST larvae, ST adults
    init_7 = [y_6(length(t_6),2);y_6(length(t_6),3);y_6(length(t_6),4);pars(3)*(1+pars(24)*(y_5(7,1)+y_5(7,2))/(1+pars(27)*(y_5(7,1)+y_5(7,2))))*y_5(7,5);y_6(length(t_6),5)];
   
    for j = 1:(length(t_6)-1)
                                         %6 %A eggs, A, LN larvae, LN adults, ST adults
        model_sol((i-1)*104+j+end_t_5,:) = [y_6(j,1);y_6(j,2);y_6(j,3);y_6(j,4);0;y_6(j,5)];
    end
    
%7 %april (2 weeks)
%7 %A, LN larvae, LN adults, ST larvae, ST adults
    %system 7
    [s_7,y_7] = ode45(@(s_7,y_7)hwaode_ALS_15system_7(s_7,y_7,pars),t_7,init_7);
    
            %8 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    init_8 = [0;y_7(length(t_7),1);y_7(length(t_7),2);y_7(length(t_7),3);y_7(length(t_7),4);y_7(length(t_7),5)];
    
    for j = 1:(length(t_7)-1)
                                        %7 %A, LN larvae, LN adults, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_6,:) = [0;y_7(j,1);y_7(j,2);y_7(j,3);y_7(j,4);y_7(j,5)];
    end
    
%8 %april, may (2 weeks)
%8 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    %system 8
    [s_8,y_8] = ode45(@(s_8,y_8)hwaode_ALS_15system_8(s_8,y_8,pars),t_8,init_8);
    
            %9 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    init_9 = [y_8(length(t_8),1);y_8(length(t_8),2);y_8(length(t_8),3);y_8(length(t_8),4);y_8(length(t_8),5);0];
    
    for j = 1:(length(t_8)-1)
                                        %8 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_7,:) = [y_8(j,1);y_8(j,2);y_8(j,3);y_8(j,4);y_8(j,5);y_8(j,6)];
    end
    
%9 %may (1 week)
%9 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
    %system 9
    [s_9,y_9] = ode45(@(s_9,y_9)hwaode_ALS_15system_9(s_9,y_9,pars),t_9,init_9);
    
            %10 %A eggs, A, ST larvae, ST adults
    init_10 = [y_9(length(t_9),1);y_9(length(t_9),2);y_9(length(t_9),5);y_9(length(t_9),6)];
    
    for j = 1:(length(t_9)-1)
                                        %9 %A eggs, A, LN larvae, LN adults, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_8,:) = [y_9(j,1);y_9(j,2);y_9(j,3);y_9(j,4);y_9(j,5);y_9(j,6)];
    end
      
%10 %may (1 week)
%10 %A eggs, A, ST larvae, ST adults
    %system 10
    [s_10,y_10] = ode45(@(s_10,y_10)hwaode_ALS_15system_10(s_10,y_10,pars),t_10,init_10);
    
            %11 %A eggs, A, ST larvae, ST adults
    init_11 = [y_10(length(t_10),1);y_10(length(t_10),2);y_10(length(t_10),3);y_10(length(t_10),4)];
    
    for j = 1:(length(t_10)-1)
                                        %10 %A eggs, A, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_9,:) = [y_10(j,1);y_10(j,2);0;0;y_10(j,3);y_10(j,4)];
    end
    
%11 %may (1 week)
%11 %A eggs, A, ST larvae, ST adults
    %system 11
    [s_11,y_11] = ode45(@(s_11,y_11)hwaode_ALS_15system_11(s_11,y_11,pars),t_11,init_11);
    
            %12 %A, ST larvae, ST adults
    init_12 = [y_11(length(t_11),2);y_11(length(t_11),3);y_11(length(t_11),4)];
    
    for j = 1:(length(t_11)-1)
                                        %11 %A eggs, A, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_10,:) = [y_11(j,1);y_11(j,2);0;0;y_11(j,3);y_11(j,4)];
    end
    
%12 %may, june (5 weeks)
%12 %A, ST larvae, ST adults
    %system 12
    [s_12,y_12] = ode45(@(s_12,y_12)hwaode_ALS_15system_12(s_12,y_12,pars),t_12,init_12);
    
            %13 %A, ST larvae, ST adults
    init_13 = [y_12(length(t_12),1);y_12(length(t_12),2);y_12(length(t_12),3)];
    
    for j = 1:(length(t_12)-1)
                                        %12 %A, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_11,:) = [0;y_12(j,1);0;0;y_12(j,2);y_12(j,3)];
    end
    
%13 %july, august (5 weeks)
%13 %A, ST larvae, ST adults
    %system 13
    [s_13,y_13] = ode45(@(s_13,y_13)hwaode_ALS_15system_13(s_13,y_13,pars),t_13,init_13);
    
            %14 %A, ST adults
    init_14 = [y_13(length(t_13),1);y_13(length(t_13),3)];
    
    for j = 1:(length(t_13)-1)
                                         %13 %A, ST larvae, ST adults
        model_sol((i-1)*104+j+end_t_12,:) = [0;y_13(j,1);0;0;y_13(j,2);y_13(j,3)];
    end
    
%14 %august, september (8 weeks)
%14 %A, ST adults
     %system 14
    [s_14,y_14] = ode45(@(s_14,y_14)hwaode_ALS_15system_14(s_14,y_14,pars),t_14,init_14);
            
            %Laricobius adults emerge from aestivation
            %15 %A, LN adults, ST adults
    init_15 = [y_14(length(t_14),1);pars(22)*y_9(length(t_9),3);y_14(length(t_14),2)];
    
    for j = 1:(length(t_14)-1)
                                         %14 %A, ST adults   
        model_sol((i-1)*104+j+end_t_13,:) = [0;y_14(j,1);0;0;0;y_14(j,2)];
    end
    
%15 %september, october (3 weeks)
%15 %A, LN adults, ST adults
	%system 15
    [s_15,y_15] = ode45(@(s_15,y_15)hwaode_ALS_15system_15(s_15,y_15,pars),t_15,init_15);
    
           %1 %A, LN adults, ST adults
    init_1 = [y_15(length(t_15),1);y_15(length(t_15),2);y_15(length(t_15),3)];
    
    for j = 1:(length(t_15)-1)
                                         %15 %A, LN adults, ST adults  
        model_sol((i-1)*104+j+end_t_14,:) = [0;y_15(j,1);0;y_15(j,2);0;y_15(j,3)];
    end

end

%saving states at final time point
                               %15 %A, LN adults, ST adults
model_sol(104*end_t_year+1,:) = [0;y_15(length(t_15),1);0;y_15(length(t_15),2);0;y_15(length(t_15),3)];

% renaming state solutions
ae_final = model_sol(:,1);
a_final = model_sol(:,2);
ll_final = model_sol(:,3);
la_final = model_sol(:,4);
sl_final = model_sol(:,5);
sa_final = model_sol(:,6);

% renaming state solutions for plotting, excluding 0s
ae_final_plot = ae_final;
ae_final_plot(ae_final_plot == 0) = nan; 
a_final_plot = a_final;
a_final_plot(a_final_plot == 0) = nan; 
ll_final_plot = ll_final;
ll_final_plot(ll_final_plot == 0) = nan; 
la_final_plot = la_final;
la_final_plot(la_final_plot == 0) = nan; 
sl_final_plot = sl_final;
sl_final_plot(sl_final_plot == 0) = nan; 
sa_final_plot = sa_final;
sa_final_plot(sa_final_plot == 0) = nan; 

% plotting all results
figure
hold on
ax = gca;
ax.FontSize = 13;
% adelgid eggs model results in black- dotted line
plot(t_full ./ 52, ae_final_plot, 'k-.', 'LineWidth', 3);
% adelgid model results in back- solid line
plot(t_full ./ 52, a_final_plot, 'k-', 'LineWidth', 3);
% LN larvae model results in orange- dashed line
plot(t_full ./ 52, ll_final_plot, '--', 'Color', [0.8, 0.4, 0], 'LineWidth', 3);
% LN adults model results in orange- solid line 
plot(t_full ./ 52, la_final_plot, '-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3);
% ST larvae model results in blue- dashed line
plot(t_full ./ 52, sl_final_plot, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 3);
% ST adults model results in blue- solid line
plot(t_full ./ 52, sa_final_plot, '-', 'Color', [0, 0.45, 0.7], 'LineWidth', 3);
% adelgid data in black- solid diamonds
plot(a_data_time ./ 104, a_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
% LN larva data in orange- open diamonds
plot(ll_data_time ./ 104, ll_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [1, 1, 1]);
% LN adult data in orange- solid diamonds
plot(la_data_time ./ 104, la_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
% ST larva data in blue- open diamonds
plot(sl_data_time ./ 104, sl_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [1, 1, 1]);
% ST adult data in blue- solid diamonds
plot(sa_data_time ./ 104, sa_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [0, 0.45, 0.7]);
xlabel('Time (years)', 'FontSize', 13)
ylabel('Density', 'FontSize', 13)
legend('A. tsugae eggs, Ae', 'A. tsugae, A', 'L. nigirinus larvae, Ll', 'L. nigrinus adults, La', 'S. tsugae larvae, Sl', 'S. tsugae adults, Sa', 'Adelgid data, A', 'L. nigirinus larvae data, Ll', 'L. nigrinus adult data, La', 'S. tsugae larvae data, Sl', 'S. tsugae adult data, Sa')

% adelgid only
figure
hold on
ax = gca;
ax.FontSize = 13;
% adelgid eggs model results in black- dotted line
plot(t_full ./ 52, ae_final_plot, 'k-.', 'LineWidth', 3);
% adelgid model results in back- solid line
plot(t_full ./ 52, a_final_plot, 'k-', 'LineWidth', 3);
% adelgid data in black- solid diamonds
plot(a_data_time ./ 104, a_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
xlabel('Time (years)', 'FontSize', 13)
ylabel('Density', 'FontSize', 13)
legend('A. tsugae eggs, Ae', 'A. tsugae, A', 'Adelgid data, A')

%plotting only predator results
%plotting all results
figure
hold on
ax = gca;
ax.FontSize = 13;
% LN larvae model results in orange- dashed line
plot(t_full ./ 52, ll_final_plot, '--', 'Color', [0.8, 0.4, 0], 'LineWidth', 3);
% LN adults model results in orange- solid line 
plot(t_full ./ 52, la_final_plot, '-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3);
% ST larvae model results in blue- dashed line
plot(t_full ./ 52, sl_final_plot, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 3);
% ST adults model results in blue- solid line
plot(t_full ./ 52, sa_final_plot, '-', 'Color', [0, 0.45, 0.7], 'LineWidth', 3);
% LN larva data in orange- open diamonds
plot(ll_data_time ./ 104, ll_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [1, 1, 1]);
% LN adult data in orange- solid diamonds
plot(la_data_time ./ 104, la_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
% ST larva data in blue- open diamonds
plot(sl_data_time ./ 104, sl_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [1, 1, 1]);
% ST adult data in blue- solid diamonds
plot(sa_data_time ./ 104, sa_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [0, 0.45, 0.7]);
xlabel('Time (years)', 'FontSize', 13)
ylabel('Density', 'FontSize', 13)
legend('L. nigirinus larvae, Ll','L. nigrinus adults, La','S. tsugae larvae, Sl','S. tsugae adults, Sa','L. nigirinus larvae data, Ll','L. nigrinus adult data, La','S. tsugae larvae data, Sl','S. tsugae adult data, Sa')

%individual species subplots
figure
subplot(3, 1, 1)
hold on
ax = gca;
ax.FontSize = 11;
%adelgid egg in black
plot(t_full ./ 52, ae_final_plot, 'k-.', 'LineWidth', 3);
%adelgid in black
plot(t_full ./ 52, a_final_plot, 'k-', 'LineWidth', 3); 
%adelgid data in black
plot(a_data_time ./ 104, a_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceColor', [0, 0, 0]);
xlabel('Time (years) starting week 43 (October)')
ylabel('A. tsugae density (per cm)')
legend('A. tsugae eggs, Ae', 'A. tsugae, A', 'A. tsugae data, A', 'Location', 'northwest')

subplot(3, 1, 2)
hold on
ax = gca;
ax.FontSize = 11;
%LN larva in orange
plot(t_full ./ 52, ll_final_plot, '--', 'Color', [0.8, 0.4, 0], 'LineWidth', 3);
%LN adult in orange
plot(t_full ./ 52, la_final_plot, '-', 'Color', [0.8, 0.4, 0], 'LineWidth', 3); 
%LN larva data in orange
plot(ll_data_time ./ 104, ll_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [1, 1, 1], 'LineWidth', 2);
%LN adult data in orange
plot(la_data_time ./ 104, la_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.8, 0.4, 0], 'MarkerFaceColor', [0.8, 0.4, 0]);
ylim([0, 9])
xlabel('Time (years) starting week 43 (October)')
ylabel('L. nigrinus density')
legend('L. nigirinus larvae, Ll', 'L. nigrinus adults, La', 'L. nigirinus larvae data, Ll', 'L. nigrinus adults data, La', 'Location', 'northwest')

subplot(3, 1, 3)
hold on
ax = gca;
ax.FontSize = 11;
%ST larva in blue
plot(t_full ./ 52, sl_final_plot, '--', 'Color', [0, 0.45, 0.7], 'LineWidth', 3); 
%ST adult in blue
plot(t_full ./ 52, sa_final_plot, '-', 'Color', [0, 0.45, 0.7], 'LineWidth', 3); 
%ST larva data in blue
plot(sl_data_time ./ 104, sl_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [1, 1, 1], 'LineWidth', 2);
%ST adult data in blue
plot(sa_data_time ./ 104, sa_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.45, 0.7], 'MarkerFaceColor', [0, 0.45, 0.7]); 
ylim([0, 7.5])
xlabel('Time (years) starting week 43 (October)')
ylabel('S. tsugae density')
legend('S. tsugae larvae, Sl', 'S. tsugae adults, Sa', 'S. tsugae larvae data, Sl', 'S. tsugae adults data, Sa', 'Location', 'northwest')

%non egg adelgids model results only
figure

hold on
ax = gca;
ax.FontSize = 13;
%adelgid in black
plot(t_full ./ 52, a_final_plot, 'k-', 'LineWidth', 3); 
xlabel('Time (years)', 'FontSize', 13)
ylabel('A. tsugae Density', 'FontSize', 13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model
% Create a function for each system of ordinary differential
% equations in the model. These functions are called when
% solving the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1 %system 1
function m = hwaode_ALS_15system_1(~, A, pars)

m = zeros(3, 1);

%dA/dt =   -m_a     * A     / Hset     - g_la     * L_a  * A    - g_sa     * S_a  * A
    m(1) = -pars(7) * A(1) ./ pars(28) - pars(19) * A(1) * A(2) - pars(21) * A(1) * A(3); %adelgids%
%dL_a/dt = -m_la     * (1 + d_a      / (A    + 0.1))      * L_a
    m(2) = -pars(12) * (1 + pars(25) * (A(1) + 0.1)^(-1)) * A(2); %LN adults%
%dS_a/dt = -m_sa     * S_a
    m(3) = -pars(14) * A(3); %ST adults%

end

%2 %system 2
function m = hwaode_ALS_15system_2(~, A, pars)

m = zeros(2, 1);

%dA/dt =   -m_aw    * A     / Hset     - g_la     * L_a  * A    
    m(1) = -pars(8) * A(1) ./ pars(28) - pars(19) * A(1) * A(2); %adelgids%
%dL_a/dt = -m_la     * (1 + d_a      / (A    + 0.1))      * L_a
    m(2) = -pars(12) * (1 + pars(25) * (A(1) + 0.1)^(-1)) * A(2); %LN adults%

end

%3 %system 3
function m = hwaode_ALS_15system_3(~, A, pars)

m = zeros(3, 1);

%dA_e/dt = r_a     * A    - m_ae    * A_e  - f_la     * A_e  * L_a
    m(1) = pars(1) * A(2) - pars(6) * A(1) - pars(16) * A(1) * A(3); %adelgid eggs%
%dA/dt =   -m_aw    * A     / Hset     - g_la     * L_a  * A    
    m(2) = -pars(8) * A(2) ./ pars(28) - pars(19) * A(2) * A(3); %adelgids%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(3) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(3); %LN adults%

end

%4 %system 4
function m = hwaode_ALS_15system_4(~, A, pars)

m = zeros(3, 1);

%dA_e/dt = r_a     * A    - c_a     * A_e  - m_ae    * A_e  - f_la     * A_e  * L_a
    m(1) = pars(1) * A(2) - pars(4) * A(1) - pars(6) * A(1) - pars(16) * A(1) * A(3); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_aw    * A     / Hset     - g_la     * L_a  * A    
    m(2) = pars(4) * A(1) - pars(8) * A(2) ./ pars(28) - pars(19) * A(2) * A(3); %adelgids%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(3) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(3); %LN adults%

end

%5 %system 5
function m = hwaode_ALS_15system_5(~, A, pars)

m = zeros(5, 1);

%dA_e/dt = r_a     * A    -c_a      * A_e  - m_ae    * A_e  - f_ll     * A_e  * L_l  - f_la     * A_e  * L_a  - f_sa     * A_e  * S_a
    m(1) = pars(1) * A(2) - pars(4) * A(1) - pars(6) * A(1) - pars(15) * A(1) * A(3) - pars(16) * A(1) * A(4) - pars(18) * A(1) * A(5); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - g_la     * A    * L_a  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(19) * A(2) * A(4) - pars(21) * A(2) * A(5); %adelgids%
%dL_l/dt = -m_ll     * (1 + d_l      / (A_e  + 0.1))      * L_l    
    m(3) = -pars(11) * (1 + pars(24) * (A(1) + 0.1)^(-1)) * A(3); %LN larvae%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(4) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(4); %LN adults%
%dS_a/dt = -m_sa     * S_a    
    m(5) = -pars(14) * A(5); %ST adults%
    
end

%6 %system 6
function m = hwaode_ALS_15system_6(~, A, pars)

m = zeros(5, 1);

%dA_e/dt = -c_a     * A_e  - m_ae    * A_e  - f_ll     * A_e  * L_l  - f_la     * A_e  * L_a  - f_sa     * A_e  * S_a
    m(1) = -pars(4) * A(1) - pars(6) * A(1) - pars(15) * A(1) * A(3) - pars(16) * A(1) * A(4) - pars(18) * A(1) * A(5); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - g_la     * A    * L_a  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(19) * A(2) * A(4) - pars(21) * A(2) * A(5); %adelgids%
%dL_l/dt = -m_ll     * (1 + d_l      / (A_e  + 0.1))      * L_l    
    m(3) = -pars(11) * (1 + pars(24) * (A(1) + 0.1)^(-1)) * A(3); %LN larvae%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(4) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(4); %LN adults%
%dS_a/dt = -m_sa     * S_a    
    m(5) = -pars(14) * A(5); %ST adults%

end

%7 %system 7
function m = hwaode_ALS_15system_7(~, A, pars)

m = zeros(5, 1);

%dA/dt =   -m_a     * A     / Hset     - g_la     * A    * L_a  - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(1) = -pars(7) * A(1) ./ pars(28) - pars(19) * A(1) * A(3) - pars(20) * A(1) * A(4) - pars(21) * A(1) * A(5); %adelgids%
%dL_l/dt = -m_ll     * L_l    
    m(2) = -pars(11) * A(2); %LN larvae%
%dL_a/dt = -m_la     * (1 + d_a      / (A    + 0.1))      * L_a
    m(3) = -pars(12) * (1 + pars(25) * (A(1) + 0.1)^(-1)) * A(3); %LN adults%
%dS_l/dt = -m_sl     * S_l    
    m(4) = -pars(13) * A(4); %ST larvae%
%dS_a/dt = -m_sa     * S_a    
    m(5) = -pars(14) * A(5); %ST adults%
    
end

%8 %system 8
function m = hwaode_ALS_15system_8(~, A, pars)

m = zeros(6, 1);

%dA_e/dt = r_a     * A(2) - c_a     * A_e  - m_ae    * A_e  - f_ll     * A_e  * L_l  - f_la     * A_e  * L_a  - f_sl     * A_e  * S_l  - f_sa     * A_e  * S_a
    m(1) = pars(1) * A(2) - pars(4) * A(1) - pars(6) * A(1) - pars(15) * A(1) * A(3) - pars(16) * A(1) * A(4) - pars(17) * A(1) * A(5) - pars(18) * A(1) * A(6); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - m_s      * A^2    - g_la     * A    * L_a  - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(10) * A(2)^2 - pars(19) * A(2) * A(4) - pars(20) * A(2) * A(5) - pars(21) * A(2) * A(6); %adelgids%
%dL_l/dt = -m_ll     * (1 + d_l      / (A_e  + 0.1))      * L_l    
    m(3) = -pars(11) * (1 + pars(24) * (A(1) + 0.1)^(-1)) * A(3); %LN larvae%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(4) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(4); %LN adults%
%dS_l/dt = -m_sl     * S_l    
    m(5) = -pars(13) * A(5); %ST larvae%
%dS_a/dt = -m_sa     * S_a    
    m(6) = -pars(14) * A(6); %ST adults%
    
end

%9 %system 9
function m = hwaode_ALS_15system_9(~, A, pars)

m = zeros(6, 1);

%dA_e/dt = r_a     * A(2) - c_a     * A_e  - m_ae    * A_e  - f_ll     * A_e  * L_l  - f_la     * A_e  * L_a  - f_sl     * A_e  * S_l  - f_sa     * A_e  * S_a
    m(1) = pars(1) * A(2) - pars(4) * A(1) - pars(6) * A(1) - pars(15) * A(1) * A(3) - pars(16) * A(1) * A(4) - pars(17) * A(1) * A(5) - pars(18) * A(1) * A(6); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - m_s      * A^2    - g_la     * A    * L_a  - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(10) * A(2)^2 - pars(19) * A(2) * A(4) - pars(20) * A(2) * A(5) - pars(21) * A(2) * A(6); %adelgids%
%dL_l/dt = -m_ll     * (1 + d_l      / (A_e  + 0.1))      * L_l    
    m(3) = -pars(11) * (1 + pars(24) * (A(1) + 0.1)^(-1)) * A(3); %LN larvae%
%dL_a/dt = -m_la     * (1 + d_a      / (A_e  + A    + 0.1))      * L_a
    m(4) = -pars(12) * (1 + pars(25) * (A(1) + A(2) + 0.1)^(-1)) * A(4); %LN adults%
%dS_l/dt = -c_s     * S_l  - m_sl     * S_l    
    m(5) = -pars(5) * A(5) - pars(13) * A(5); %ST larvae%
%dS_a/dt = c_s     * S_l  - m_sa     * S_a    
    m(6) = pars(5) * A(5) - pars(14) * A(6); %ST adults%
    
end

%10 %system 10
function m = hwaode_ALS_15system_10(~, A, pars)

m = zeros(4, 1);

%dA_e/dt = r_a     * A(2) -c_a      * A_e  - m_ae    * A_e  - f_sl     * A_e  * S_l  - f_sa     * A_e  * S_a
    m(1) = pars(1) * A(2) - pars(4) * A(1) - pars(6) * A(1) - pars(17) * A(1) * A(3) - pars(18) * A(1) * A(4); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - m_s      * A^2    - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(10) * A(2)^2 - pars(20) * A(2) * A(3) - pars(21) * A(2) * A(4); %adelgids%
%dS_l/dt = -c_s     * S_l  - m_sl     * S_l    
    m(3) = -pars(5) * A(3) - pars(13) * A(3); %ST larvae%
%dS_a/dt = c_s     * S_l  - m_sa     * S_a    
    m(4) = pars(5) * A(3) - pars(14) * A(4); %ST adults%

end

%11 %system 11
function m = hwaode_ALS_15system_11(~, A, pars)

m = zeros(4, 1);

%dA_e/dt = -c_a     * A_e  - m_ae    * A_e  - f_sl     * A_e  * S_l  - f_sa     * A_e  * S_a
    m(1) = -pars(4) * A(1) - pars(6) * A(1) - pars(17) * A(1) * A(3) - pars(18) * A(1) * A(4); %adelgid eggs%
%dA/dt =   c_a     * A_e  - m_a     * A     / Hset     - m_s      * A^2    - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(2) = pars(4) * A(1) - pars(7) * A(2) ./ pars(28) - pars(10) * A(2)^2 - pars(20) * A(2) * A(3) - pars(21) * A(2) * A(4); %adelgids%
%dS_l/dt = -c_s     * S_l  - m_sl     * S_l    
    m(3) = -pars(5) * A(3) - pars(13) * A(3); %ST larvae%
%dS_a/dt = c_s     * S_l  - m_sa     * S_a    
    m(4) = pars(5) * A(3) - pars(14) * A(4); %ST adults%

end

%12 %system 12
function m = hwaode_ALS_15system_12(~, A, pars)

m = zeros(3, 1);

%dA/dt =   -m_a     * A     / Hset     - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(1) = -pars(7) * A(1) ./ pars(28) - pars(20) * A(1) * A(2) - pars(21) * A(1) * A(3); %adelgids%
%dS_l/dt = -c_s     * S_l  - m_sl     * S_l    
    m(2) = -pars(5) * A(2) - pars(13) * A(2); %ST larvae%
%dS_a/dt = c_s     * S_l  - m_sa     * S_a    
    m(3) = pars(5) * A(2) - pars(14) * A(3); %ST adults%

end

%13 %system 13
function m = hwaode_ALS_15system_13(~, A, pars)

m = zeros(3, 1);

%dA/dt =   -m_as    * A     / Hset     - g_sl     * A    * S_l  - g_sa     * A    * S_a
    m(1) = -pars(9) * A(1) ./ pars(28) - pars(20) * A(1) * A(2) - pars(21) * A(1) * A(3); %adelgids%
%dS_l/dt = -c_s     * S_l  - m_sl     * S_l    
    m(2) = -pars(5) * A(2) - pars(13) * A(2); %ST larvae%
%dS_a/dt = c_s     * S_l  - m_sa     * S_a    
    m(3) = pars(5) * A(2) - pars(14) * A(3); %ST adults%
    
end

%14 %system 14
function m = hwaode_ALS_15system_14(~, A, pars)

m = zeros(2, 1);

%dA/dt =   -m_as    * A     / Hset     - g_sa     * A    * S_a
    m(1) = -pars(9) * A(1) ./ pars(28) - pars(21) * A(1) * A(2); %adelgids%
%dS_a/dt = -m_sa     * S_a    
    m(2) = -pars(14) * A(2); %ST adults%
    
end

%15 %system 15
function m = hwaode_ALS_15system_15(~, A, pars)

m = zeros(3, 1);

%dA/dt =   -m_as    * A     / Hset     - g_la     * A    * L_a  -g_sa      * A    * S_a
    m(1) = -pars(9) * A(1) ./ pars(28) - pars(19) * A(1) * A(2) - pars(21) * A(1) * A(3); %adelgids%
%dL_a/dt = -m_la     * (1 + d_a      /  (A + 0.1))        * L_a
    m(2) = -pars(12) * (1 + pars(25) * (A(1) + 0.1)^(-1)) * A(2); %LN adults%    
%dS_a/dt = -m_sa     * S_a    
    m(3) = -pars(14) * A(3); %ST adults%
    
end
