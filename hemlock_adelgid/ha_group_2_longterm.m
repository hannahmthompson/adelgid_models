%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T. canadensis and A. tsugae model
% Plot of long term behavior
% Group I parameter estimation results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save model results that correspond to optimal parameter 
% choices found by parameter estimation 
% Plotted for 200 years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% controls steepness of transition from positive to negative tips alive growth
b_1 = 9.9999;
% controls threshold of between positive and negative tips alive growth
b_2 = 0.3617;
% related to threshold and symmetry of recovery and decay
l = 0.1346;
% tips alive growth rate (value relative to l controls symmetry of recovery and decay)
g_h = 0.1500;
% adelgid growth rate
g_a = 0.6000;
% adelgid death rate due to sexuparae
m_s = 0.0803;
% background per capita adelgid death rate 
m_a = 0.0613;
% winter per capita adelgid death rate
m_aw = 0.0739;
% summer per capita adelgid death rate
m_as = 0.0010;
% tips alive carrying capacity
k = 0.9398;

% create vector of parameters
%       1    2    3  4    5    6    7     8     9    10   
pars = [b_1  b_2  l  g_h  g_a  m_a  m_aw  m_as  m_s  k];

% set times
% total number of years
end_t_year = 200;
% set time vectors for individual systems with (number of weeks) * 2 + 1 entries
t_1 = linspace(0, 1, 3);
t_2 = linspace(0, 1, 3);
t_3 = linspace(0, 16, 33);
t_4 = linspace(0, 5, 11);
t_5 = linspace(0, 14, 29);
t_6 = linspace(0, 1, 3);
t_7 = linspace(0, 5, 11);
t_8 = linspace(0, 1, 3);
t_9 = linspace(0, 5, 11);
t_10 = linspace(0, 3, 7);
t_full = linspace(0, end_t_year * 52, end_t_year * 104 + 1);

% set initial conditions
% hemlock initial condition, proportion of tips alive
h_0 = 0.661111111;
% other adelgid initial condition, density in hwa/cm
a_0 = 2.481481481;
% create a vector of initial conditions
init_1 = [h_0; a_0];

% initalize matrix to store full model run solutions
model_sol = zeros(104 * end_t_year + 1, 2);

% yearly loop
for i = 1 : end_t_year

    % initialize/reset matrices to store solutions for each system 
    y_1 = zeros(length(t_1), 2);
    y_2 = zeros(length(t_2), 2);
    % dH/dt=0 for 3rd system, so only adelgid solution is stored
    y_3 = zeros(length(t_3), 1);
    y_4 = zeros(length(t_4), 2);
    y_5 = zeros(length(t_5), 2);
    y_6 = zeros(length(t_6), 2);
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
    
    init_3 = y_2(length(t_2), 2);
    % dH/dt = 0 for system 3 so H = hem_3
    hem_3 = y_2(length(t_2), 1);
    
    for j = 1 : (length(t_2) - 1)
        model_sol((i - 1) * 104 + j + 2, :) = y_2(j, :);
    end
    
    % system 3
    [s_3, y_3] = ode45(@(s_3, y_3)hwaode_HA_tipsalive_6(s_3, y_3, pars, hem_3), t_3, init_3);
    
    init_4 = [hem_3; y_3(length(t_3), 1)];
    
    for j = 1 : (length(t_3) - 1)
        model_sol((i - 1) * 104 + j + 4, :) = [hem_3; y_3(j, :)];
    end
    
    % system 4
    [s_4, y_4] = ode45(@(s_4, y_4)hwaode_HA_tipsalive_2(s_4, y_4, pars), t_4, init_4);
    
    init_5 = [y_4(length(t_4), 1); y_4(length(t_4), 2)];
    
    for j = 1 : (length(t_4) - 1)
        model_sol((i - 1) * 104 + j + 36, :) = y_4(j, :);
    end
    
    % system 5
    [s_5, y_5] = ode45(@(s_5, y_5)hwaode_HA_tipsalive_8(s_5, y_5, pars), t_5, init_5);
    
    init_6 = [y_5(length(t_5), 1); y_5(length(t_5), 2)];
    
    for j = 1 : (length(t_5) - 1)
        model_sol((i - 1) * 104 + j + 46, :) = y_5(j, :);
    end
    
    % system 6
    [s_6, y_6] = ode45(@(s_6, y_6)hwaode_HA_tipsalive_2(s_6, y_6, pars), t_6, init_6);
    
    init_7 = [y_6(length(t_6), 1); y_6(length(t_6), 2)];
    
    for j = 1 : (length(t_6) - 1)
        model_sol((i - 1) * 104 + j + 74, :) = y_6(j, :);
    end
    
    % system 7
    [s_7, y_7] = ode45(@(s_7, y_7)hwaode_HA_tipsalive_1(s_7, y_7, pars), t_7, init_7);
    
    init_8 = [y_7(length(t_7), 1); y_7(length(t_7), 2)];
    
    for j = 1 : (length(t_7) - 1)
        model_sol((i - 1) * 104 + j + 76, :) = y_7(j, :);
    end
    
    % system 8
    [s_8, y_8] = ode45(@(s_8, y_8)hwaode_HA_tipsalive_2(s_8, y_8, pars), t_8, init_8);
    
    init_9 = [y_8(length(t_8), 1); y_8(length(t_8), 2)];
    
    for j = 1 : (length(t_8) - 1)
        model_sol((i - 1) * 104 + j + 86, :) = y_8(j, :);
    end
    
    % system 9
    [s_9, y_9] = ode45(@(s_9, y_9)hwaode_HA_tipsalive_3(s_9, y_9, pars), t_9, init_9);
    
    init_10 = [y_9(length(t_9), 1); y_9(length(t_9), 2)];
    
    for j = 1 : (length(t_9) - 1)
        model_sol((i - 1) * 104 + j + 88, :) = y_9(j, :);
    end
    
    % system 10
    [s_10, y_10] = ode45(@(s_10, y_10)hwaode_HA_tipsalive_1(s_10, y_10, pars), t_10, init_10);
    
    init_1 = [y_10(length(t_10), 1); y_10(length(t_10), 2)];
    
    for j = 1 : (length(t_10) - 1)
        model_sol((i - 1) * 104 + j + 98, :) = y_10(j, :);
    end

end

model_sol(104 * end_t_year + 1, :) = y_10(length(t_10), :);

%renaming state solutions
h_final = model_sol(:, 1);
a_final = model_sol(:, 2);

% plot results
figure() 

subplot(2, 1, 1)
hold on
ax = gca;
ax.FontSize = 16;
% tips alive in green
% model results in solid line, data as diamonds
plot(t_full ./ 52, h_final, '-', 'Color', [0, 0.6, 0.5], 'LineWidth', 1.15)
ylim([0, 1])
xlabel('Time (years) starting week 27 (July)', 'FontSize', 16)
ylabel('Proportion tips alive', 'FontSize', 16)

subplot(2, 1, 2)
hold on
ax = gca;
ax.FontSize = 16;
% adelgid density in black
% model results in solid line, data as diamonds
plot(t_full ./ 52, a_final, 'k-', 'LineWidth', 1.15)
xlabel('Time (years) starting week 27 (July)', 'FontSize', 16)
ylabel('A. tsugae density (per cm)', 'FontSize', 16)

% save figure as .jpg file
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 14, 6], 'PaperUnits', 'Inches', 'PaperSize', [14, 6])
f = gcf;
exportgraphics(f, 'ha_group_2_longterm.jpg', 'Resolution', 600)

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