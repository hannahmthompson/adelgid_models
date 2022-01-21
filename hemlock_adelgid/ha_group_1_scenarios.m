%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T. canadensis and A. tsugae model
% Testing varied initial condition scenarios
% Group I parameter estimation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% controls steepness of transition from positive to negative tips alive growth
b_1 = 9.9974;
% controls threshold of between positive and negative tips alive growth
b_2 = 2.4262;
% related to threshold and symmetry of recovery and decay
l = 0.0915;
% tips alive growth rate (value relative to l controls symmetry of recovery and decay)
g_h = 0.1500;
% adelgid growth rate
g_a = 0.27;
% adelgid death rate due to sexuparae
m_s = 0.0399;
% background per capita adelgid death rate 
m_a = 0.0005;
% winter per capita adelgid death rate
m_aw = 0.0010;
% summer per capita adelgid death rate
m_as = 0.0412;
% tips alive carrying capacity
k = 0.8789;

% create vector of parameters
%       1    2    3  4    5    6    7     8     9    10   
pars = [b_1  b_2  l  g_h  g_a  m_a  m_aw  m_as  m_s  k];

% set times
% total number of years
end_t_year = 15;
% set time vectors for individual systems with (number of weeks) * 2 + 1 entries
t_1 = linspace(0, 3, 7);
t_2 = linspace(0, 1, 3);
t_3 = linspace(0, 5, 11);
t_4 = linspace(0, 4, 9);
t_5 = linspace(0, 1, 3);
t_6 = linspace(0, 16, 33);
t_7 = linspace(0, 5, 11);
t_8 = linspace(0, 14, 29);
t_9 = linspace(0, 1, 3);
t_10 = linspace(0, 2, 5);
t_full = linspace(0, end_t_year * 52, end_t_year * 104 + 1);

% initalize matrix to store full mode=l run solutions 
% for all scenarios (total of 1215 scenarios tested)
all_model_sol = zeros(104 * end_t_year + 1, 2 * 1215);

% save scenario results in matrices where each row
% represents one of 1215 scenario: 
% A(0) values from 0 to 8 hwa/cm in increments of 0.1
% H(0) values from 0.15 to 0.85 proportion tips alive
% in increments of 0.05
% if model results reach 0.1 proportion tips alive
% (I) 1 will replace 0 for that entry in deadorno matrix
deadorno = zeros(1215, 3);
% (II) first time model results reach 0.1 proportion tips 
% alive will replace 0 for that entry in deadtime matrix
deadtime = zeros(1215, 3);

% loop through each value of A(0)
% within loop, loop through each value of H(0)
% save results in all_model_sol and appropriate entry of
% deadorno and deadtime
for Acount = 1 : 81
    
    % calculate adelgid initial condition (in hwa/cm)
    a_0 = Acount / 10 - 0.1;
    
    % loop through each value of H(0)
    for Hcount = 1 : 15

        % calculate tips alive initial condition
        h_0 = Hcount / 20 + 0.1; 

        % create a vector of initial conditions
        init_1 = [h_0; a_0];

        % initalize matrix to store full model run solutions 
        % for this scenario
        model_sol = zeros(104 * end_t_year + 1, 2);
        
        % yearly loop
        for i = 1 : end_t_year
        
            % initialize/reset matrices to store solutions for each system 
            y_1 = zeros(length(t_1), 2);
            y_2 = zeros(length(t_2), 2);
            y_3 = zeros(length(t_3), 2);
            y_4 = zeros(length(t_4), 2);
            y_5 = zeros(length(t_5), 2);
            % dH/dt=0 for 6th system, so only adelgid solution is stored
            y_6 = zeros(length(t_6), 1);
            y_7 = zeros(length(t_7), 2);
            y_8 = zeros(length(t_8), 2);
            y_9 = zeros(length(t_9), 2);
            y_10 = zeros(length(t_10), 2);
            
            % solve system 1 and store results 
            [s_1, y_1] = ode45(@(s_1, y_1)hwaode_HA_tipsalive_1(s_1, y_1, pars), t_1, init_1);
            
            % set initial conditions for next system 
            init_2 = [y_1(length(t_1), 1); y_1(length(t_1), 2)];
            
            % store results in full model run solution matrix
            for j = 1 : (length(t_1) - 1)
                model_sol((i - 1) * 104 + j, :) = y_1(j, :);
            end
            
            % system 2
            [s_2, y_2] = ode45(@(s_2, y_2)hwaode_HA_tipsalive_2(s_2, y_2, pars), t_2, init_2);
            
            init_3 = [y_2(length(t_2), 1); y_2(length(t_2), 2)];
            
            for j = 1 : (length(t_2) - 1)
                model_sol((i - 1) * 104 + j + 6, :) = y_2(j, :);
            end
            
            % system 3
            [s_3, y_3] = ode45(@(s_3, y_3)hwaode_HA_tipsalive_3(s_3, y_3, pars), t_3, init_3);
            
            init_4 = [y_3(length(t_3), 1); y_3(length(t_3), 2)];
            
            for j = 1 : (length(t_3) - 1)
                model_sol((i - 1) * 104 + j + 8, :) = y_3(j, :);
            end
            
            % system 4
            [s_4, y_4] = ode45(@(s_4, y_4)hwaode_HA_tipsalive_1(s_4, y_4, pars), t_4, init_4);
            
            init_5 = [y_4(length(t_4), 1); y_4(length(t_4), 2)];
            
            for j = 1 :(length(t_4) - 1)
                model_sol((i - 1) * 104 + j + 18, :) = y_4(j, :);
            end
            
            % system 5
            [s_5, y_5] = ode45(@(s_5, y_5)hwaode_HA_tipsalive_2(s_5, y_5, pars), t_5, init_5);
            
            init_6 = y_5(length(t_5), 2);
            % dH/dt = 0 for system 6 so H = hem_6
            hem_6 = y_5(length(t_5), 1);
            
            for j = 1 : (length(t_5) - 1)
                model_sol((i - 1) * 104 + j + 26, :) = y_5(j, :);
            end
            
            % system 6
            [s_6, y_6] = ode45(@(s_6, y_6)hwaode_HA_tipsalive_6(s_6, y_6, pars, hem_6), t_6, init_6);
            
            init_7 = [hem_6; y_6(length(t_6), 1)];
            
            for j = 1 : (length(t_6) - 1)
                model_sol((i - 1) * 104 + j + 28, :) = [hem_6; y_6(j, :)];
            end
            
            % system 7
            [s_7, y_7] = ode45(@(s_7, y_7)hwaode_HA_tipsalive_2(s_7, y_7, pars), t_7, init_7);
            
            init_8 = [y_7(length(t_7), 1); y_7(length(t_7), 2)];
            
            for j = 1 : (length(t_7) - 1)
                model_sol((i - 1) * 104 + j + 60, :) = y_7(j, :);
            end
            
            % system 8
            [s_8, y_8] = ode45(@(s_8, y_8)hwaode_HA_tipsalive_8(s_8, y_8, pars), t_8, init_8);
            
            init_9 = [y_8(length(t_8), 1); y_8(length(t_8), 2)];
            
            for j = 1 : (length(t_8) - 1)
                model_sol((i - 1) * 104 + j + 70, :) = y_8(j, :);
            end
            
            % system 9
            [s_9,y_9] = ode45(@(s_9, y_9)hwaode_HA_tipsalive_2(s_9, y_9, pars), t_9, init_9);
            
            init_10 = [y_9(length(t_9), 1); y_9(length(t_9), 2)];
            
            for j = 1 : (length(t_9) - 1)
                model_sol((i - 1) * 104 + j + 98, :) = y_9(j, :);
            end
            
            % system 10
            [s_10, y_10] = ode45(@(s_10, y_10)hwaode_HA_tipsalive_1(s_10, y_10, pars), t_10,  init_10);
            
            init_1=[y_10(length(t_10), 1 ); y_10(length(t_10), 2)];
            
            for j=1 : (length(t_10) - 1)
                model_sol((i - 1) * 104 + j + 100, :)=y_10(j, :);
            end
        
        end
        
        % save solution for last time point 
        model_sol(104 * end_t_year + 1, :) = y_10(length(t_10), :);
        
        % renaming state solutions
        h_final = model_sol(:, 1);
        a_final = model_sol(:, 2);
        
        % save initial conditions in result matrices
        deadorno(15 * (Acount - 1) + Hcount, 1) = h_0;
        deadorno(15 * (Acount - 1) + Hcount, 2) = a_0;
        deadtime(15 * (Acount - 1) + Hcount, 1) = h_0;
        deadtime(15 * (Acount - 1) + Hcount, 2) = a_0;

        % if model solution ever reaches 0.1 proportion tips alive
        if any(h_final <= 0.1)
            % change entry of deadorno from 0 to 1
            deadorno(15 * (Acount - 1) + Hcount, 3) = 1;
            % change entry of deadtime to first time model reachs 0.1
            % proportion tips alive
            deadtime(15 * (Acount - 1) + Hcount, 3) = t_full(min(find(h_final <= 0.1))) ./ 52;
        end
        
        %saving model solutions for these initial conditions
        all_model_sol(:, 30 * (Acount - 1) + 2 * Hcount - 1) = model_sol(:, 1);
        all_model_sol(:, 30 * (Acount - 1) + 2 * Hcount) = model_sol(:, 2);

    end    
    
end

% save result matrices as csv files 
writematrix(deadorno, 'ha_group_1_deadorno.csv') 
writematrix(deadtime, 'ha_group_1_deadtime.csv')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model
% Create a function for each system of ordinary differential
% equations in the model. 
% These functions are called when solving the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function m = hwaode_HA_tipsalive_1(t, A, pars)

    m = zeros(2, 1);
    
 % dH/dt = g_h       (-1  / (1 + e^{ -b_1      (A     - b_2)}      +  l)         H     (1 - H     / k)
    m(1) = pars(4) * (-1 ./ (1 + exp(-pars(1) * (A(2) - pars(2)))) + pars(3)) * A(1) * (1 - A(1) ./ pars(10));
 % dA/dt = g_a       A    - m_a       A     / H
    m(2) = pars(5) * A(2) - pars(6) * A(2) ./ A(1);

end

function m = hwaode_HA_tipsalive_2(t, A, pars)
    
    m = zeros(2, 1);
    
 % dH/dt = g_h       (-1  / (1 + e^{ -b_1       (A    - b_2)}      +  l)        H      (1 - H     / k)
    m(1) = pars(4) * (-1 ./ (1 + exp(-pars(1) * (A(2) - pars(2)))) + pars(3)) * A(1) * (1 - A(1) ./ pars(10));
 % dA/dt = - m_a      A     / H
    m(2) = -pars(6) * A(2) ./ A(1);

end

function m = hwaode_HA_tipsalive_3(t, A, pars)

    m = zeros(2, 1);
 
 % dH/dt = g_h       (-1  / (1 + e^{  -b_1        (A    - b_2)}      +  l)        H      (1 - H     / k)
    m(1) = pars(4) * (-1 ./ (1 + exp( - pars(1) * (A(2) - pars(2)))) + pars(3)) * A(1) * (1 - A(1) ./ pars(10));
 % dA/dt = - m_a      A     / H    - m_s       A^2
    m(2) = -pars(6) * A(2) ./ A(1) - pars(9) * A(2)^2;

end

function m = hwaode_HA_tipsalive_6(t, A, pars, hem_6)

    m = zeros(1, 1);
    
 % dH/dt = 0
 % dA/dt = - m_as     A    / H
    m(1) = -pars(8) * A(1) / hem_6;

end

function m = hwaode_HA_tipsalive_8(t, A, pars)

    m = zeros(2, 1);
   
 % dH/dt = g_h       (-1  / (1 + e^{ -b_1       (A    - b_2)}      + l)         H      (1 - H     / k)
    m(1) = pars(4) * (-1 ./ (1 + exp(-pars(1) * (A(2) - pars(2)))) + pars(3)) * A(1) * (1 - A(1) ./ pars(10));
 % dA/dt = - m_aw     A     / H
    m(2) = -pars(7) * A(2) ./ A(1);

end