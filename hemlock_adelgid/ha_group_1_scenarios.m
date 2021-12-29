
b_1=9.997352069;%controls steepness of transition from positive to negative tips alive growth
b_2=2.426159708;%controls threshold of between positive and negative tips alive growth
l=0.091471924;%related to threshold and symmetry of recovery and decay
g_hp=0.149997195;%tips alive growth rate (value relative to l controls symmetry of recovery and decay)
g_ap=0.27;%adelgid growth rate
m_ap=0.000500969;%background per capita adelgid death rate 
m_awp=0.001002748;%winter per capita adelgid death rate
m_asp=0.041234489;%summer per capita adelgid death rate
m_sp=0.039882022;%adelgid death rate due to sexuparae
k=0.878886384;

%       1  2  3   4    5   6     7     8    9   10   
pars=[b_1 b_2 l g_hp g_ap m_ap m_awp m_asp m_sp k];

end_t_year=30;
t_1=linspace(0,3,7);
t_2=linspace(0,1,3);
t_3=linspace(0,5,11);
t_4=linspace(0,4,9);
t_5=linspace(0,1,3);
t_6=linspace(0,16,33);
t_7=linspace(0,5,11);
t_8=linspace(0,14,29);
t_9=linspace(0,1,3);
t_10=linspace(0,2,5);
t_year=linspace(0,52,105);
t_full=linspace(0,end_t_year*52,end_t_year*104+1);

all_model_sol=zeros(104*end_t_year+1, 930);
deadorno=zeros(162,15);
deadtime=zeros(162,15);

for Acount=1:162
    
    a_0=Acount/10-0.1;  %other adelgid initial condition, density in hwa/cm

    
    for Hcount=1:15

        %initial conditions
        h_0=Hcount/20+.1; %hemlock initial condition, proportion of tips alive  

        init_1=[h_0;a_0]; %create a vector of initial conditions

        model_sol=zeros(104*end_t_year+1,2);

        %weekly timepoints in tspan
        for i=1:end_t_year

            %set solution   
            y_1=zeros(length(t_1),2);
            y_2=zeros(length(t_2),2);
            y_3=zeros(length(t_3),2);
            y_4=zeros(length(t_4),2);
            y_5=zeros(length(t_5),2);
            y_6=zeros(length(t_6),1);
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

            init_3=[y_2(length(t_2),1);y_2(length(t_2),2)];

            for j=1:(length(t_2)-1)
                model_sol((i-1)*104+j+6,:)=y_2(j,:);
            end




            [s_3,y_3] = ode45(@(s_3,y_3)hwaode_HA_tipsalive_3(s_3,y_3,pars),t_3,init_3);

            init_4=[y_3(length(t_3),1);y_3(length(t_3),2)];

            for j=1:(length(t_3)-1)
                model_sol((i-1)*104+j+8,:)=y_3(j,:);
            end



            [s_4,y_4] = ode45(@(s_4,y_4)hwaode_HA_tipsalive_1(s_4,y_4,pars),t_4,init_4);

            init_5=[y_4(length(t_4),1);y_4(length(t_4),2)];

            for j=1:(length(t_4)-1)
                model_sol((i-1)*104+j+18,:)=y_4(j,:);
            end



            [s_5,y_5] = ode45(@(s_5,y_5)hwaode_HA_tipsalive_2(s_5,y_5,pars),t_5,init_5);

            init_6=y_5(length(t_5),2);
            hem_6=y_5(length(t_5),1);

            for j=1:(length(t_5)-1)
                model_sol((i-1)*104+j+26,:)=y_5(j,:);
            end



            [s_6,y_6] = ode45(@(s_6,y_6)hwaode_HA_tipsalive_6(s_6,y_6,pars, hem_6),t_6,init_6);

            init_7=[hem_6;y_6(length(t_6),1)];

            for j=1:(length(t_6)-1)
                model_sol((i-1)*104+j+28,:)=[hem_6;y_6(j,:)];
            end



            [s_7,y_7] = ode45(@(s_7,y_7)hwaode_HA_tipsalive_2(s_7,y_7,pars),t_7,init_7);

            init_8=[y_7(length(t_7),1);y_7(length(t_7),2)];

            for j=1:(length(t_7)-1)
                model_sol((i-1)*104+j+60,:)=y_7(j,:);
            end



            [s_8,y_8] = ode45(@(s_8,y_8)hwaode_HA_tipsalive_8(s_8,y_8,pars),t_8,init_8);

            init_9=[y_8(length(t_8),1);y_8(length(t_8),2)];

            for j=1:(length(t_8)-1)
                model_sol((i-1)*104+j+70,:)=y_8(j,:);
            end



            [s_9,y_9] = ode45(@(s_9,y_9)hwaode_HA_tipsalive_2(s_9,y_9,pars),t_9,init_9);

            init_10=[y_9(length(t_9),1);y_9(length(t_9),2)];

            for j=1:(length(t_9)-1)
                model_sol((i-1)*104+j+98,:)=y_9(j,:);
            end



            [s_10,y_10] = ode45(@(s_10,y_10)hwaode_HA_tipsalive_1(s_10,y_10,pars),t_10,init_10);

            init_1=[y_10(length(t_10),1);y_10(length(t_10),2)];

            for j=1:(length(t_10)-1)
                model_sol((i-1)*104+j+100,:)=y_10(j,:);
            end

        end

        model_sol(104*end_t_year+1,:)=y_10(length(t_10),:);
  

        %renaming state solutions
        h_final=model_sol(:,1);
        a_final=model_sol(:,2);
        
        if any(h_final<=0.1)
            deadorno(Acount,Hcount)=1;
            deadtime(Acount,Hcount)=t_full(min(find(h_final<=0.1)))./52;
        end
        
        %saving model solutions for these initial conditions
        all_model_sol(:,30*(Acount-1)+2*Hcount-1)=model_sol(:,1);
        all_model_sol(:,30*(Acount-1)+2*Hcount)=model_sol(:,2);

%         figure
%         subplot(2,1,1)
%         hold on
%         ax = gca;
%         ax.FontSize = 16;
%         plot(t_full./52,h_final,'-','Color',[0, 0.6, .5], 'LineWidth', 3); %tips alive in green
%         xlabel('Time (years) starting week 15 (April)', 'FontSize',16)
%         ylabel('Proportion tips alive','FontSize',16)
%         ylim([0,1])
%         text(.1, .1, sprintf('H0= %1.1f', h_0), 'FontSize',16) 
%         text(2.1, .1, sprintf('A0= %1.1f', a_0),'FontSize',16)
% 
%         subplot(2,1,2)
%         hold on
%         ax = gca;
%         ax.FontSize = 16;
%         plot(t_full./52,a_final,'k-','LineWidth', 3); %aldegid density in black
%         xlabel('Time (years) starting week 15 (April)', 'FontSize',16)
%         ylabel('A. tsugae density (per cm)','FontSize',16)
%         
%         hold off
%         
    end
    
    
end



function m=hwaode_HA_tipsalive_1(t, A, pars)

m=zeros(2,1);

m(1)=pars(4)*(-1./(1+exp(-pars(1)*(A(2)-pars(2))))+pars(3))*A(1)*(1-A(1)./pars(10)); %hemlock health
m(2)=pars(5)*A(2)-pars(6)*A(2)./A(1); %adelgid density

end

function m=hwaode_HA_tipsalive_2(t, A, pars)

m=zeros(2,1);

m(1)=pars(4)*(-1./(1+exp(-pars(1)*(A(2)-pars(2))))+pars(3))*A(1)*(1-A(1)./pars(10)); %hemlock health
m(2)=-pars(6)*A(2)./A(1); %adelgid density

end

function m=hwaode_HA_tipsalive_3(t, A, pars)

m=zeros(2,1);

m(1)=pars(4)*(-1./(1+exp(-pars(1)*(A(2)-pars(2))))+pars(3))*A(1)*(1-A(1)./pars(10)); %hemlock health
m(2)=-pars(6)*A(2)./A(1)-pars(9)*A(2)^2; %adelgid density

end

function m=hwaode_HA_tipsalive_6(t, A, pars, hem_6)

m=zeros(1,1);

%m(1)=0; %hemlock health
m(1)=-pars(8)*A(1)/hem_6; %adelgid density

end

function m=hwaode_HA_tipsalive_8(t, A, pars)

m=zeros(2,1);

m(1)=pars(4)*(-1./(1+exp(-pars(1)*(A(2)-pars(2))))+pars(3))*A(1)*(1-A(1)./pars(10)); %hemlock health
m(2)=-pars(7)*A(2)./A(1); %adelgid density

end