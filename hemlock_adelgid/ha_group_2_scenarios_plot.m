b_1=9.999905788;%controls steepness of transition from positive to negative tips alive growth
b_2=0.36167139;%controls threshold of between positive and negative tips alive growth
l=0.134569601;%related to threshold and symmetry of recovery and decay
g_hp=0.149999345;%tips alive growth rate (value relative to l controls symmetry of recovery and decay)
g_ap=0.599999043;%adelgid growth rate
m_ap=0.06133784;%background per capita adelgid death rate 
m_awp=0.07389841;%winter per capita adelgid death rate
m_asp=0.001000503;%summer per capita adelgid death rate
m_sp=0.080323418;%adelgid death rate due to sexuparae
k=0.939770364;


%       1  2  3   4    5   6     7     8    9   10   
pars=[b_1 b_2 l g_hp g_ap m_ap m_awp m_asp m_sp k];

end_t_year=15;
t_1=linspace(0,1,3);
t_2=linspace(0,1,3);
t_3=linspace(0,16,33);
t_4=linspace(0,5,11);
t_5=linspace(0,14,29);
t_6=linspace(0,1,3);
t_7=linspace(0,5,11);
t_8=linspace(0,1,3);
t_9=linspace(0,5,11);
t_10=linspace(0,3,7);
t_year=linspace(0,52,105);
t_full=linspace(0,end_t_year*52,end_t_year*104+1);

init_test=[0.15  0.5 
           0.15  1.4
           0.15  1.5
           0.15  3.0
           0.20  2.1
           0.20  2.2
           0.20  2.3
           0.20  3.2
           0.20  3.7
           0.20  3.8
           0.20  5.0
           0.25  1.0
           0.25  3.0
           0.25  3.1
           0.25  3.2  
           0.25  4.0
           0.25  4.4
           0.25  4.5
           0.25  4.6
           0.25  4.7
           0.25  5.5
           0.30  0.5
           0.30  4.0
           0.30  5.6
           0.30  5.7
           0.30  5.8
           0.30  5.9
           0.30  7.0
           0.35  1.0
           0.35  3.0
           0.35  5.0
           0.35  7.0
           0.40  1.0
           0.40  3.0
           0.40  5.0
           0.40  7.0
           0.45  1.0
           0.45  3.0
           0.45  5.0
           0.45  7.0
           0.50  1.0
           0.50  3.0
           0.50  5.0
           0.50  7.0
           0.55  1.0
           0.55  3.0
           0.55  5.0
           0.55  7.0
           0.60  1.0
           0.60  3.0
           0.60  5.0
           0.60  7.0
           0.65  1.0
           0.65  3.0
           0.65  5.0
           0.65  7.0
           0.70  1.0
           0.70  3.0
           0.70  5.0
           0.70  7.0
           0.75  0.5
           0.75  3.0
           0.75  6.0
           0.75  7.8
           0.75  7.9
           0.75  8.0
           0.80  2.0
           0.80  5.2
           0.80  5.3
           0.80  5.4
           0.80  7.0
           0.85  0.5
           0.85  1.2
           0.85  1.5
           0.85  1.6
           0.85  1.7
           0.85  2.0
           0.85  2.3
           0.85  2.4
           0.85  2.9
           0.85  3.1
           0.85  3.2
           0.85  5.0
           0.85  6.8
           0.85  7.4
           0.85  7.5
           0.85  8.0
           0.90  1.0
           0.90  3.8
           0.90  3.9
           0.90  5.0
           0.90  7.0
           0.90  8.0];

for q=1:length(init_test)
    
    h_0=init_test(q,1);
    a_0=init_test(q,2);

        init_1=[h_0;a_0]; %create a vector of initial conditions

        model_sol=zeros(104*end_t_year+1,2);

        %weekly timepoints in tspan
        for i=1:end_t_year

            %set solution   
            y_1=zeros(length(t_1),2);
            y_2=zeros(length(t_2),2);
            y_3=zeros(length(t_3),1);
            y_4=zeros(length(t_4),2);
            y_5=zeros(length(t_5),2);
            y_6=zeros(length(t_6),2);
            y_7=zeros(length(t_7),2);
            y_8=zeros(length(t_8),2);
            y_9=zeros(length(t_9),2);
            y_10=zeros(length(t_10),2);



            [s_1,y_1] = ode45(@(s_1,y_1)hwaode_HA_tipsalive_1(s_1,y_1,pars),t_1,init_1);

            init_2=[y_1(length(t_1),1);y_1(length(t_1),2)];

            for j=1:(length(t_1)-1)
                model_sol((i-1)*104+j,:)=y_1(j,:);
            end



            [s_2,y_2] = ode45(@(s_2,y_2)hwaode_HA_tipsalive_2(s_2,y_2,pars),t_2,init_2);

            init_3=y_2(length(t_2),2);
            hem_3=y_2(length(t_2),1);

            for j=1:(length(t_2)-1)
                model_sol((i-1)*104+j+2,:)=y_2(j,:);
            end




            [s_3,y_3] = ode45(@(s_3,y_3)hwaode_HA_tipsalive_6(s_3,y_3,pars, hem_3),t_3,init_3);

            init_4=[hem_3;y_3(length(t_3),1)];

            for j=1:(length(t_3)-1)
                model_sol((i-1)*104+j+4,:)=[hem_3;y_3(j,:)];
            end



            [s_4,y_4] = ode45(@(s_4,y_4)hwaode_HA_tipsalive_2(s_4,y_4,pars),t_4,init_4);

            init_5=[y_4(length(t_4),1);y_4(length(t_4),2)];

            for j=1:(length(t_4)-1)
                model_sol((i-1)*104+j+36,:)=y_4(j,:);
            end



            [s_5,y_5] = ode45(@(s_5,y_5)hwaode_HA_tipsalive_8(s_5,y_5,pars),t_5,init_5);

            init_6=[y_5(length(t_5),1);y_5(length(t_5),2)];

            for j=1:(length(t_5)-1)
                model_sol((i-1)*104+j+46,:)=y_5(j,:);
            end



            [s_6,y_6] = ode45(@(s_6,y_6)hwaode_HA_tipsalive_2(s_6,y_6,pars),t_6,init_6);

            init_7=[y_6(length(t_6),1);y_6(length(t_6),2)];

            for j=1:(length(t_6)-1)
                model_sol((i-1)*104+j+74,:)=y_6(j,:);
            end



            [s_7,y_7] = ode45(@(s_7,y_7)hwaode_HA_tipsalive_1(s_7,y_7,pars),t_7,init_7);

            init_8=[y_7(length(t_7),1);y_7(length(t_7),2)];

            for j=1:(length(t_7)-1)
                model_sol((i-1)*104+j+76,:)=y_7(j,:);
            end



            [s_8,y_8] = ode45(@(s_8,y_8)hwaode_HA_tipsalive_2(s_8,y_8,pars),t_8,init_8);

            init_9=[y_8(length(t_8),1);y_8(length(t_8),2)];

            for j=1:(length(t_8)-1)
                model_sol((i-1)*104+j+86,:)=y_8(j,:);
            end



            [s_9,y_9] = ode45(@(s_9,y_9)hwaode_HA_tipsalive_3(s_9,y_9,pars),t_9,init_9);

            init_10=[y_9(length(t_9),1);y_9(length(t_9),2)];

            for j=1:(length(t_9)-1)
                model_sol((i-1)*104+j+88,:)=y_9(j,:);
            end



            [s_10,y_10] = ode45(@(s_10,y_10)hwaode_HA_tipsalive_1(s_10,y_10,pars),t_10,init_10);

            init_1=[y_10(length(t_10),1);y_10(length(t_10),2)];

            for j=1:(length(t_10)-1)
                model_sol((i-1)*104+j+98,:)=y_10(j,:);
            end
            
        end

        model_sol(104*end_t_year+1,:)=y_10(length(t_10),:);
  

        %renaming state solutions
        h_final=model_sol(:,1);
        a_final=model_sol(:,2);
        
        
        death = ones(1,length(t_full))*0.1;
      
        figure
        subplot(2,1,1)
        hold on
        ax = gca;
        ax.FontSize = 16;
        plot(t_full./52,h_final,'-','Color',[0, 0.6, .5], 'LineWidth', 3); %tips alive in green
        plot(t_full./52,death,'-','Color',[0, 0.4, .3], 'LineWidth', 2)
        xlabel('Time (years) starting week 27 (July)', 'FontSize',16)
        ylabel('Proportion tips alive','FontSize',16)
        ylim([0,1])
        text(.1, 0.92, sprintf('H(0)= %1.2f', h_0), 'FontSize',14) 
        text(3.1, 0.92, sprintf('A(0)= %1.1f', a_0),'FontSize',14)

        subplot(2,1,2)
        hold on
        ax = gca;
        ax.FontSize = 16;
        plot(t_full./52,a_final,'k-','LineWidth', 3); %aldegid density in black
        xlabel('Time (years) starting week 27 (July)', 'FontSize',16)
        ylabel('A. tsugae density (per cm)','FontSize',16)
        
        hold off
    
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