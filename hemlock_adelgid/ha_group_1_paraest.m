%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T. canadensis and A. tsugae model
% Parameter estimation
% Group I data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up optimization problem- set bounds on parameters, 
% starting values of parameters, optimization problem, 
% multistart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear figures and workspace
clf;
clear all;

% set bounds for parameters to be estimated
%               b_1   b_2   l     g_h     g_a   m_s    m_a     m_aw   m_as   k    
lower_bounds = [0.01  0.01  0.01  0.0001  0.27  0.005  0.0005  0.001  0.001  0.775];
upper_bounds = [10    8     1     0.15    0.27  4      0.2     0.3    0.3    0.95];

% number of starting points
num_start_points = 10; 

% set one starting point
x_start = 0.5 * (lower_bounds + upper_bounds);

% create optimization problem structure and set options
problem = createOptimProblem('fmincon', 'objective', @HA_tipsalive_10system, 'x0', x_start, 'lb', lower_bounds, 'ub', upper_bounds);
problem.options = optimoptions(problem.options, 'MaxFunctionEvaluations', 99999, 'MaxIterations', 99999);
% lower these values for testing/to make problem move on if it gets stuck

% create a solver for optimization problem 
% print output from each starting point
mult_start = MultiStart('Display', 'iter');
% if using parallel computing:
% mult_start=MultiStart('UseParallel',true,'Display', 'iter');
% parpool

% run optimization problem and save results 
[x, yval, exitflag, output, solutions] = run(mult_start, problem, num_start_points);

% if using parallel computing:
% delete(gcp)

% initialize matrix of results from all starting points
r = size(solutions);
results = zeros(r(1, 2),  length(upper_bounds) + 3);

% save results- exit flag, first order optimality, objective function value, parameter values
for i = 1 : r(1, 2)
   results(i, :) = [solutions(1,  i).Exitflag,  solutions(1,  i).Output.bestfeasible.firstorderopt,  solutions(1, i).Fval solutions(1, i).X];
end

% displaying result- x contains the parameter values, 
% yval is the ojective function vlaue
x
yval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save model results that correspond to optimal parameter 
% choices found 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign entries of x to corresponding parameters
% controls steepness of transition from positive to negative tips alive growth
b_1 = x(1);
% controls threshold of between positive and negative tips alive growth
b_2 = x(2);
% related to threshold and symmetry of recovery and decay
l = x(3);
% tips alive growth rate (value relative to l controls symmetry of recovery and decay)
g_h = x(4);
% adelgid growth rate
g_a = x(5);
% adelgid death rate due to sexuparae
m_s = x(6);
% background per capita adelgid death rate 
m_a = x(7);
% winter per capita adelgid death rate
m_aw = x(8);
% summer per capita adelgid death rate
m_as = x(9);
% tips alive carrying capacity
k = x(10);

% create vector of parameters
%       1    2    3  4    5    6    7     8     9    10   
pars = [b_1  b_2  l  g_h  g_a  m_a  m_aw  m_as  m_s  k];

% set times
% total number of years
end_t_year = 13;
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
t_full = linspace(0, end_t_year * 52 ,end_t_year * 104 + 1);

% set initial conditions
% hemlock initial condition, proportion of tips alive
h_0 = 0.775; 
% other adelgid initial condition, density in hwa/cm
a_0 = 2.844166667;  
% create a vector of initial conditions
init_1 = [h_0; a_0]; 

% initalize matrix to store full model run solutions
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

% rename state solutions
h_final = model_sol(:, 1);
a_final = model_sol(:, 2);

% Group I data
% time in calendar weeks
time_entries_a_data = [59	77	119	129	160	181	216	234	277	285	323	338	375	389]';
a_data = [1.9125	3.845	0.050833333	0.005	0.02	0.011111111	0.258333333	0.016666667	0	0.041666667	0.079166667	0.429166667	0.133333333	2.070833333];
time_entries_h_data = [59	77	119	160	181	234	285	338	389	675]';
h_data = [0.54	0.57	0.117	0.202	0.12	0.525	0.55	0.525	0.55	0.575];

% test data- same time points as group 1 data
% comment out either real data or this test data

% test- estimated parameter should agree with parameters
% used to generate test data.

% pars value used to generate test data:
% [5.6 0.13 0.53 0.11 0.4 0.071 0.0005 0.05 0.05 0.9]
% a_data = [1.51787754722614 2.95828807074122 0.0736788762675378 0.329409756029034 0.00659696162793616 0.0157221343699808 0.00163464685975224 0.0168772199769714 0.0129197552243110 0.0283872307099463 0.0119080488438500 0.176290302258266 0.0586720244564968 0.526021161227084];
% h_data = [0.540505188299063 0.339935316522196 0.251571515875157 0.248795645473372 0.343073352151633 0.527374975741411 0.678578768417368 0.768655293520170 0.757395091245183 0.654430982698062];

% to test that model solution here matches model solution in 
% calculation of objective function value function,
% calculate objective function value here
a_diff = zeros(length(a_data), 1);
h_diff = zeros(length(h_data), 1);
    
for k = 1 : length(a_data)
    a_diff(k) = a_final(2 * (time_entries_a_data(k) - 15) + 1) - a_data(k);
end

for j = 1 : length(h_data)
    h_diff(j) = h_final(2 * (time_entries_h_data(j) - 15) + 1) - h_data(j);  
end

% calculate and output objective function value
norm(h_diff) / norm(h_data) + norm(a_diff) / norm(a_data)

% plot results
figure() 

subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 16;
% tips alive in green
% model results as solid line, data as diamonds
plot(t_full ./ 52, h_final, '-', 'Color', [0, 0.6, 0.5], 'LineWidth', 3)
plot((time_entries_h_data - 15) ./ 52, h_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0, 0.7, 0.6], 'MarkerFaceColor', [0, 0.7, 0.6])
xlabel('Time (years) starting week 15 (April)', 'FontSize', 16)
ylabel('Proportion tips alive', 'FontSize', 16)

subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 16;
% adelgid density in black
% model results as solid line, data as diamonds
plot(t_full ./ 52, a_final, 'k-', 'LineWidth',  3)
plot((time_entries_a_data - 15) ./ 52, a_data, 'd', 'MarkerSize', 8, 'MarkerEdgeColor', [0.1, 0.1, 0.1], 'MarkerFaceColor', [0.1, 0.1, 0.1])
xlabel('Time (years) starting week 15 (April)', 'FontSize', 16)
ylabel('A. tsugae density (per cm)', 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function value function
% Solve the model and calculate and return the objective 
% function value for a given set of parameter values z. 
% The values of z are determined by MultiStart and fmincon.
% System solved is similar to the above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = HA_tipsalive_10system(z)

    b_1 = z(1);
    b_2 = z(2);
    l = z(3);
    g_h = z(4);
    g_a = z(5);
    m_s = z(6);
    m_a = z(7);
    m_aw = z(8);
    m_as = z(9);
    k = z(10);

    pars = [b_1  b_2  l  g_h  g_a  m_a  m_aw  m_as  m_s  k];
    
    end_t_year = 13;
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
    
    h_0 = 0.775; 
    a_0 = 2.844166667;  
    init_1 = [h_0; a_0]; 
    
    model_sol = zeros(104 * end_t_year + 1, 2);
    
    for i = 1 : end_t_year
        
        y_1 = zeros(length(t_1), 2);
        y_2 = zeros(length(t_2), 2);
        y_3 = zeros(length(t_3), 2);
        y_4 = zeros(length(t_4), 2);
        y_5 = zeros(length(t_5), 2);
        y_6 = zeros(length(t_6), 1);
        y_7 = zeros(length(t_7), 2);
        y_8 = zeros(length(t_8), 2);
        y_9 = zeros(length(t_9), 2);
        y_10 = zeros(length(t_10), 2);
        
        [s_1, y_1] = ode45(@(s_1, y_1)hwaode_HA_tipsalive_1(s_1, y_1, pars), t_1, init_1);
        
        init_2 = [y_1(length(t_1), 1); y_1(length(t_1), 2)];
        
        for j = 1 : (length(t_1) - 1)
            model_sol((i - 1) * 104 + j, :) = y_1(j, :);
        end
        
        [s_2, y_2] = ode45(@(s_2, y_2)hwaode_HA_tipsalive_2(s_2, y_2, pars), t_2, init_2);
        
        init_3 = [y_2(length(t_2), 1); y_2(length(t_2), 2)];
        
        for j = 1 : (length(t_2) - 1)
            model_sol((i - 1) * 104 + j + 6, :) = y_2(j, :);
        end
        
        [s_3, y_3] = ode45(@(s_3, y_3)hwaode_HA_tipsalive_3(s_3, y_3, pars), t_3, init_3);
        
        init_4 = [y_3(length(t_3), 1); y_3(length(t_3), 2)];
        
        for j = 1 : (length(t_3) - 1)
            model_sol((i - 1) * 104 + j + 8, :) = y_3(j, :);
        end
        
        [s_4, y_4] = ode45(@(s_4, y_4)hwaode_HA_tipsalive_1(s_4, y_4, pars), t_4, init_4);
        
        init_5 = [y_4(length(t_4), 1); y_4(length(t_4), 2)];
        
        for j = 1 : (length(t_4) - 1)
            model_sol((i - 1) * 104 + j + 18, :) = y_4(j, :);
        end
        
        [s_5, y_5] = ode45(@(s_5, y_5)hwaode_HA_tipsalive_2(s_5, y_5, pars), t_5, init_5);
        
        init_6 = y_5(length(t_5), 2);
        hem_6 = y_5(length(t_5), 1);
        
        for j = 1 : (length(t_5) - 1)
            model_sol((i - 1) * 104 + j + 26, :) = y_5(j, :);
        end
        
        [s_6, y_6] = ode45(@(s_6, y_6)hwaode_HA_tipsalive_6(s_6, y_6, pars, hem_6), t_6, init_6);
        
        init_7 = [hem_6; y_6(length(t_6), 1)];
        
        for j = 1 : (length(t_6) - 1)
            model_sol((i - 1) * 104 + j + 28, :) = [hem_6; y_6(j, :)];
        end
        
        [s_7, y_7] = ode45(@(s_7, y_7)hwaode_HA_tipsalive_2(s_7, y_7, pars), t_7, init_7);
        
        init_8 = [y_7(length(t_7), 1); y_7(length(t_7), 2)];
        
        for j = 1 : (length(t_7) - 1)
            model_sol((i - 1) * 104 + j + 60, :) = y_7(j, :);
        end
        
        [s_8, y_8] = ode45(@(s_8, y_8)hwaode_HA_tipsalive_8(s_8, y_8, pars), t_8, init_8);
        
        init_9 = [y_8(length(t_8), 1); y_8(length(t_8), 2)];
        
        for j = 1 : (length(t_8) - 1)
            model_sol((i - 1) * 104 + j + 70, :) = y_8(j, :);
        end
        
        [s_9,y_9] = ode45(@(s_9, y_9)hwaode_HA_tipsalive_2(s_9, y_9, pars), t_9, init_9);
        
        init_10 = [y_9(length(t_9), 1); y_9(length(t_9), 2)];
        
        for j = 1 : (length(t_9) - 1)
            model_sol((i - 1) * 104 + j + 98, :) = y_9(j, :);
        end
        
        [s_10, y_10] = ode45(@(s_10, y_10)hwaode_HA_tipsalive_1(s_10, y_10, pars), t_10,  init_10);
        
        init_1=[y_10(length(t_10), 1 ); y_10(length(t_10), 2)];
        
        for j=1 : (length(t_10) - 1)
            model_sol((i - 1) * 104 + j + 100, :)=y_10(j, :);
        end

    end
    
    model_sol(104 * end_t_year + 1, :) = y_10(length(t_10), :);
    
    h_final = model_sol(:, 1);
    a_final = model_sol(:, 2);

    % reset objective function value
    value = 0;
    
    % Group I data
    time_entries_a_data = [59	77	119	129	160	181	216	234	277	285	323	338	375	389]';
    a_data = [1.9125	3.845	0.050833333	0.005	0.02	0.011111111	0.258333333	0.016666667	0	0.041666667	0.079166667	0.429166667	0.133333333	2.070833333];
    time_entries_h_data = [59	77	119	160	181	234	285	338	389	675]';
    h_data = [0.54	0.57	0.117	0.202	0.12	0.525	0.55	0.525	0.55	0.575];

    % test data- same time points as group I data
    % a_data = [1.51787754722614 2.95828807074122 0.0736788762675378 0.329409756029034 0.00659696162793616 0.0157221343699808 0.00163464685975224 0.0168772199769714 0.0129197552243110 0.0283872307099463 0.0119080488438500 0.176290302258266 0.0586720244564968 0.526021161227084];
    % h_data = [0.540505188299063 0.339935316522196 0.251571515875157 0.248795645473372 0.343073352151633 0.527374975741411 0.678578768417368 0.768655293520170 0.757395091245183 0.654430982698062];

    % initialize vectors to store error at each time point with data
    a_diff = zeros(length(a_data), 1);
    h_diff = zeros(length(h_data), 1);
    
    % calculate error at each time point with data
    % for adelgid data 
    for k = 1 : length(a_data)
        a_diff(k) = a_final(2 * (time_entries_a_data(k) - 15) + 1) - a_data(k);
    end
    
    % for hemlock data
    for j = 1 : length(h_data)
        h_diff(j) = h_final(2 * (time_entries_h_data(j) - 15) + 1) - h_data(j);  
    end
    
    % calculate objective function value- sum of squares of relative errors
    value = norm(h_diff) / norm(h_data) + norm(a_diff) / norm(a_data);

    if isnan(value)
        value = 0;
    end

end

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