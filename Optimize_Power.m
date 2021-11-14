clear all, close all
clc
%yalmip('clear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Parameters - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = 100; % height
Qe = 40; % radius of Eve
Qcu = 20; % radius of cognitive user
Qpu = 20; % radius of primary user
x0 = -100; % initial location - x
y0 = 200; % initial location - y
ht0 = H; % initial location -ht
xF = 500; % final location - x
yF = 200; % final location - y
htF = H; % final location -ht
x_cs = 0; % cognitive source location - x
y_cs = 0; % cognitive source location - y
ht_cs = 0; % cognitive source location - ht
x_cu = 300; % cognitive user location -x
y_cu = 0; % cognitive user location -y
ht_cu = 0; % cognitive user location -ht
x_e = 150; % eavesdropper location -x
y_e = 250; % evaesdropper location -y
ht_e = 0; % eveasdropper location -ht
x_pu = 0; % primary user location -x
y_pu = 250; % primary user location -y
ht_pu = 0; % primary user location -ht
Pcs_average = 10^(40/10-3)/2; % the average transmit power
Pcs_max = 10^(40/10-3); % the peak transmit power
Pu_average = (10^(4/10-3))/2; % the average jamming signal power
Pu_max = (10^(4/10-3)); % the peak jamming signal power
V = 10; % the max velocity of the UAV
l = 0; % iteration index
noise = 10^(-70/10-3); % the power spectral density of the AWGN
beta = 10^(10/10-3); % the channel power gain at the reference distance d = 1m
gamma = 10^8; % gamma = beta/noise = 90dB
phi = 3; % the path loss exponent
index = 0.5772156;
d_cs_cu = (sqrt((x_cs-x_cu)^2+(y_cs-y_cu)^2) + Qcu)^(-phi); % the distance between the cognitive source and the cognitive user
d_cs_e = (sqrt((x_cs-x_e)^2+(y_cs-y_e)^2) - Qe)^(-phi); % the distance between the cognitive source and cognitive eavesdropper
d_cs_pu = (sqrt((x_cs-x_pu)^2+(y_cs-y_pu)^2) - Qpu)^(-phi); % the distance between the cognitive source and the primary user

%%%%%%%%%%%%%%%%%%%%%%%%% - SecrecyRate_Iteration - %%%%%%%%%%%%%%%%%%%%%%%
fid_SercrecyRate_Threshold = fopen('..\40m\Output\Threshold\SecrecyRate_Iteration\Optimize_Power\Secrecy_Rate.txt','w');
fid_Iteration_Threshold = fopen('..\40m\Output\Threshold\SecrecyRate_Iteration\Optimize_Power\Iteration.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%% - SecrecyRate_Time - %%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_SecrecyRate = fopen('..\40m\Output\SecrecyRate_Time\Optimize_Power\Secrecy_Rate.txt','w');%
fid_Time = fopen('..\40m\Output\SecrecyRate_Time\Optimize_Power\Time.txt','w');

fid_x_1 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_1.txt','w');
fid_y_1 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_1.txt','w');
fid_ht_1 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_1.txt','w');
fid_x_2 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_2.txt','w');
fid_y_2 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_2.txt','w');
fid_ht_2 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_2.txt','w');
fid_x_3 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_3.txt','w');
fid_y_3 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_3.txt','w');
fid_ht_3 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_3.txt','w');
fid_x_4 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_4.txt','w');
fid_y_4 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_4.txt','w');
fid_ht_4 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_4.txt','w');
fid_x_5 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_5.txt','w');
fid_y_5 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_5.txt','w');
fid_ht_5 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_5.txt','w');
fid_x_6 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_6.txt','w');
fid_y_6 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_6.txt','w');
fid_ht_6 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_6.txt','w');
fid_x_7 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_7.txt','w');
fid_y_7 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_7.txt','w');
fid_ht_7 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_7.txt','w');
fid_x_8 = fopen('..\40m\Output\Trajectory\Optimize_Power\x_8.txt','w');
fid_y_8 = fopen('..\40m\Output\Trajectory\Optimize_Power\y_8.txt','w');
fid_ht_8 = fopen('..\40m\Output\Trajectory\Optimize_Power\ht_8.txt','w');

fid_PowerSource_1 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_1.txt','w');
fid_TimeSlot1_1 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot1.txt','w');
fid_PowerSource_2 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_2.txt','w');
fid_TimeSlot1_2 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot2.txt','w');
fid_PowerSource_3 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_3.txt','w');
fid_TimeSlot1_3 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot3.txt','w');
fid_PowerSource_4 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_4.txt','w');
fid_TimeSlot1_4 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot4.txt','w');
fid_PowerSource_5 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_5.txt','w');
fid_TimeSlot1_5 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot5.txt','w');
fid_PowerSource_6 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_6.txt','w');
fid_TimeSlot1_6 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot6.txt','w');
fid_PowerSource_7 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_7.txt','w');
fid_TimeSlot1_7 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot7.txt','w');
fid_PowerSource_8 = fopen('..\40m\Output\Power_Source\Optimize_Power\Power_Source_8.txt','w');
fid_TimeSlot1_8 = fopen('..\40m\Output\Power_Source\Optimize_Power\Time_Slot8.txt','w');

fid_PowerUAV_1 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_1.txt','w');
fid_TimeSlot2_1 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot1.txt','w');
fid_PowerUAV_2 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_2.txt','w');
fid_TimeSlot2_2 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot2.txt','w');
fid_PowerUAV_3 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_3.txt','w');
fid_TimeSlot2_3 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot3.txt','w');
fid_PowerUAV_4 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_4.txt','w');
fid_TimeSlot2_4 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot4.txt','w');
fid_PowerUAV_5 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_5.txt','w');
fid_TimeSlot2_5 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot5.txt','w');
fid_PowerUAV_6 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_6.txt','w');
fid_TimeSlot2_6 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot6.txt','w');
fid_PowerUAV_7 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_7.txt','w');
fid_TimeSlot2_7 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot7.txt','w');
fid_PowerUAV_8 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Power_UAV_8.txt','w');
fid_TimeSlot2_8 = fopen('..\40m\Output\Power_UAV\Optimize_Power\Time_Slot8.txt','w');


fid_Rate_1 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_1.txt','w');
fid_TimeSlot3_1 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot1.txt','w');
fid_Rate_2 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_2.txt','w');
fid_TimeSlot3_2 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot2.txt','w');
fid_Rate_3 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_3.txt','w');
fid_TimeSlot3_3 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot3.txt','w');
fid_Rate_4 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_4.txt','w');
fid_TimeSlot3_4 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot4.txt','w');
fid_Rate_5 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_5.txt','w');
fid_TimeSlot3_5 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot5.txt','w');
fid_Rate_6 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_6.txt','w');
fid_TimeSlot3_6 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot6.txt','w');
fid_Rate_7 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_7.txt','w');
fid_TimeSlot3_7 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot7.txt','w');
fid_Rate_8 = fopen('..\40m\Output\Rate\Optimize_Power\Rate_8.txt','w');
fid_TimeSlot3_8 = fopen('..\40m\Output\Rate\Optimize_Power\Time_Slot8.txt','w');

fid_Distance_1 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_1.txt','w');
fid_TimeSlot4_1 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot1.txt','w');
fid_Distance_2 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_2.txt','w');
fid_TimeSlot4_2 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot2.txt','w');
fid_Distance_3 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_3.txt','w');
fid_TimeSlot4_3 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot3.txt','w');
fid_Distance_4 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_4.txt','w');
fid_TimeSlot4_4 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot4.txt','w');
fid_Distance_5 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_5.txt','w');
fid_TimeSlot4_5 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot5.txt','w');
fid_Distance_6 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_6.txt','w');
fid_TimeSlot4_6 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot6.txt','w');
fid_Distance_7 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_7.txt','w');
fid_TimeSlot4_7 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot7.txt','w');
fid_Distance_8 = fopen('..\40m\Output\Distance\Optimize_Power\Distance_8.txt','w');
fid_TimeSlot4_8 = fopen('..\40m\Output\Distance\Optimize_Power\Time_Slot8.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_flight = 80; % the UAV's trajectory versus period
%N = 2*T_flight; % the number of time slots
N = 500;
k = 1;
while (1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Initial points - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x_initial = zeros(1,N); % the initial jamming trajectory - x                         
    y_initial = zeros(1,N); % the initial  jamming trajectory - y
    ht_initial = zeros(1,N); % the initial jamming trajectory - ht
    for i = 1:1:N
        x_initial(i) = x0 + (i)*(xF-x0)/(N);
    end
    for i = 1:1:N
        y_initial(i) = y0;
    end
    for i = 1:1:N
        ht_initial(i) = H;
    end
    Pcs_initial_1 = zeros(1,N);
    Pu_initial_1 = zeros(1,N);
    Pcs_initial = zeros(1,N);
    Pu_initial = zeros(1,N);
    for i = 1:1:N
        Pcs_initial(1,i) = Pcs_average;
        Pu_initial(1,i) = Pu_average;
    end

    % - the initial slack variables - %
    z_initial = (x_cu-x_initial).^2 + (y_initial).^2 + ht_initial.^2;

    t_initial = exp(-index)*gamma*d_cs_cu*Pcs_initial.*z_initial./(gamma*Pu_initial+z_initial);

    k_initial = (x_e-x_initial).^2 + (y_e-y_initial).^2 + ht_initial.^2;

    u_initial = gamma*d_cs_e.*Pcs_initial.*k_initial./(gamma*Pu_initial+k_initial);

    t_0_initial = 1./t_initial;

    v_initial = gamma*Pu_initial./k_initial;

    w_initial = 1./(v_initial+1);

    q_initial = (x_pu-x_initial).^2 + (y_pu-y_initial).^2 + ht_initial.^2;

    h_initial = 1./q_initial;
    
    a = zeros(1,N);
    a(1)=1;
    for i = 2:1:N
        a(i)= a(i-1) + 1;
    end
    %R_new=double(sum(d(1,1:N)-e(1,1:N)));
    cvx_solver SDPT3
    scale_factor = 10^5;
    R_cvx_update = zeros(1,N);
    %while (1)
    for i = 1:1:15
        A = (gamma*Pu_initial+z_initial)./t_initial; 
        B = k_initial./v_initial;
        C = (Pcs_initial)./w_initial;
        D = Pu_initial./h_initial;
        %R_old=R_new;
        cvx_begin sdp
            variables R_cvx(1,N) Pu_cvx(1,N) Pcs_cvx(1,N) x_cvx(1,N);  
            variables epsilon_cvx_e(1,N) schur_cvx_e(1,N) epsilon_cvx_cu(1,N) schur_cvx_cu(1,N) epsilon_cvx_pu(1,N) schur_cvx_pu(1,N);
            variables z_cvx(1,N) t_cvx(1,N) d_cvx(1,N) t_0_cvx(1,N) e_cvx(1,N) k_cvx(1,N) u_cvx(1,N) v_cvx(1,N) w_cvx(1,N) c_cvx(1,N) q_cvx(1,N) h_cvx(1,N);    
            maximize(sum(d_cvx(1,1:N)-e_cvx(1,1:N))/N) %
            subject to
                %%% - Trajectory constraints - %%%
                (x_cvx(1)-x0) <= T_flight*V/N;
                (x_cvx(1)-x0) >= -T_flight*V/N;
                for n = 1:1:N-1
                    (x_cvx(n+1)-x_cvx(n)) <= T_flight*V/N;
                    (x_cvx(n+1)-x_cvx(n)) >= -T_flight*V/N;
                    x_cvx(n) <= xF;
                    x_cvx(n) >= x0;
                end
                x_cvx(N) == xF;
      
                %%% - Power constraints - %%%
                (1/N)*sum(Pu_cvx(1,1:N)) <= Pu_average;
                (1/N)*sum(Pcs_cvx(1,1:N)) <= Pcs_average;
            
                %%% - Cognitive radio constraint - %%%
                % scale_factor*((1/N)*sum(beta*d_cs_pu*Pcs_cvx(1,1:N)+c_cvx(1,1:N))) <= scale_factor*10^(-5); % 5 6 7 standard 5
            
                for n = 1:1:N
                    %%% - Power constraints - %%%
                    Pu_cvx(n) >= 0;
                    Pu_cvx(n) <= Pu_max;
                    Pcs_cvx(n) >= 0;
                    Pcs_cvx(n) <= Pcs_max;
                
                    %%% - Rd constraints - %%%
                    [(epsilon_cvx_cu(n)+1) 0 (x_cu-x_cvx(n));...
                     0 (epsilon_cvx_cu(n)+1) (y_cu-y_initial(n));... 
                     (x_cu-x_cvx(n)) (y_cu-y_initial(n)) (-Qcu^2*epsilon_cvx_cu(n)-schur_cvx_cu(n))] >= 0;
                    epsilon_cvx_cu(n) >= 0;
                    z_cvx(n) - schur_cvx_cu(n) <= (x_cu-x_initial(n))^2 + 2*(x_initial(n)-x_cu)*(x_cvx(n)-x_initial(n)) + (y_initial(n))^2 + ht_initial(n)^2;
                    % scale_factor*z_cvx(n) <= scale_factor*((x_cu-x_initial(n))^2+2*(x_initial(n)-x_cu)*(x_cvx(n)-x_initial(n))+(y_initial(n))^2+(ht_initial(n))^2);
                    scale_factor*norm([t_cvx(n)*sqrt(0.5*A(n)) (gamma*Pu_cvx(n)+z_cvx(n))*sqrt(0.5/A(n)) (Pcs_cvx(n)-z_cvx(n))*(sqrt(exp(-index)*gamma*d_cs_cu))/2])... 
                    <= scale_factor*(Pcs_cvx(n)+z_cvx(n))*(sqrt(exp(-index)*gamma*d_cs_cu))/2;%%%Add Euler constant
                    scale_factor*d_cvx(n) <= scale_factor*((log2(1+1/t_0_initial(n)))-(t_0_cvx(n)-t_0_initial(n))/((t_0_initial(n)*(t_0_initial(n)+1)*log(2))));
                    scale_factor*norm([1 0.5*(t_cvx(n)-t_0_cvx(n))]) <= scale_factor*0.5*(t_cvx(n)+t_0_cvx(n));
                
                    %%% - Re constraints - %%%
                    [(epsilon_cvx_e(n)-1) 0 (x_cvx(n)-x_e);...
                     0 (epsilon_cvx_e(n)-1) (y_initial(n)-y_e);... 
                     (x_cvx(n)-x_e) (y_initial(n)-y_e) (-Qe^2*epsilon_cvx_e(n)-schur_cvx_e(n))] >= 0;
                    epsilon_cvx_e(n) >= 0;
                    norm([(x_cvx(n)-x_e) (y_initial(n)-y_e) (ht_initial(n)-ht_e) (1/2)*(k_cvx(n)+schur_cvx_e(n)-1)]) <= (1/2)*(k_cvx(n)+schur_cvx_e(n)+1);
                    %scale_factor*norm([(-x_e+x_cvx(n)) (-y_e+y_initial(n)) ht_initial(n) (1/2)*(k_cvx(n)-1)]) <= scale_factor*(1/2)*(k_cvx(n)+1); % k(n) >= (x_e-x(n))^2+(y_e-y(n))^2+H^2,
                    scale_factor*norm([w_cvx(n)*sqrt(gamma*0.5*d_cs_e*C(n)) (Pcs_cvx(n))*sqrt(gamma*0.5*d_cs_e/C(n)) 0.5*(u_cvx(n)-1)]) <= scale_factor*0.5*(u_cvx(n)+1);
                    scale_factor*norm([1 0.5*(w_cvx(n)-v_cvx(n)-1)]) <= scale_factor*0.5*(w_cvx(n)+v_cvx(n)+1);
                    scale_factor*norm([v_cvx(n)*sqrt(0.5*B(n)) k_cvx(n)*sqrt(0.5/B(n)) 0.5*(Pu_cvx(n)-1)*sqrt(gamma)]) <= scale_factor*0.5*(Pu_cvx(n)+1)*sqrt(gamma);
                    scale_factor*e_cvx(n) >= scale_factor*((log2(1+u_initial(n))+(u_cvx(n)-u_initial(n))/((1+u_initial(n))*log(2))));
                
                    %%% - Cognitive radio constraints - %%%
                    [(epsilon_cvx_pu(n)+1) 0 (x_pu-x_cvx(n));...
                     0 (epsilon_cvx_pu(n)+1) (y_pu-y_initial(n));... 
                     (x_pu-x_cvx(n)) (y_pu-y_initial(n)) (-Qpu^2*epsilon_cvx_pu(n)-schur_cvx_pu(n))] >= 0;
                    epsilon_cvx_pu(n) >= 0;
                    q_cvx(n) - schur_cvx_pu(n) <= (x_pu-x_initial(n))^2 + 2*(x_initial(n)-x_pu)*(x_cvx(n)-x_initial(n)) + (y_pu-y_initial(n))^2 + (ht_pu-ht_initial(n))^2;
                    % scale_factor*q_cvx(n) <= scale_factor*((x_pu-x_initial(n))^2+2*(x_initial(n)-x_pu)*(x_cvx(n)-x_initial(n))+(y_pu-y_initial(n))^2+(ht_pu-ht_initial(n))^2); 
                    scale_factor*norm([h_cvx(n)*sqrt(0.5*beta*D(n)) Pu_cvx(n)*sqrt(0.5*beta/D(n)) 0.5*(c_cvx(n)-1)]) <= scale_factor*0.5*(c_cvx(n)+1);
                    scale_factor*norm([1 0.5*(h_cvx(n)-q_cvx(n))]) <= scale_factor*0.5*(h_cvx(n)+q_cvx(n));
                    scale_factor*(beta*d_cs_pu*Pcs_cvx(n) + c_cvx(n)) <= scale_factor*10^(-5);
                
                    %%% - Implicit constraints - %%%
                    Pcs_cvx(n)+z_cvx(n) >= 0; % Rd constraint
                    t_0_cvx(n) >= 0; %  Rd constraint
                    t_cvx(n)+t_0_cvx(n) >= 0; % Rd constraint
                    k_cvx(n) >= H^2; % Re constraint
                    v_cvx(n) >= 0; % Re constraint
                    w_cvx(n) >= 0; % Re constraint
                    w_cvx(n)+v_cvx(n)+1 >= 0; % Re constraint
                    Pu_cvx(n) >= 0; % 
                    h_cvx(n)+q_cvx(n) >= 0; % 
                    c_cvx(n) >= 0; % 
                
                end
        cvx_end
        %app = log2((gamma*d_cs_e*Pcs_cvx./(1+gamma*Pu_cvx./k_cvx))+1);
        %disp(e_cvx - log2(u_cvx + 1));
    
        %%% - Plot - %%%
        figure(1)
        hold on
        plot(x_initial,y_initial);
        title('Trajectories of the UAV')
        xlabel('x(m)')
        ylabel('y(m)')
    
        %%% - Update - %%%
        delta_cvx = x_cvx - x_initial;
        x_initial = x_cvx;
        t_0_initial = t_0_cvx;
        t_initial = t_cvx;
        Pu_initial = Pu_cvx;
        z_initial = z_cvx;
        k_initial = k_cvx;
        u_initial = u_cvx;
        v_initial = v_cvx;
        w_initial = w_cvx;
        h_initial = h_cvx;
        Pcs_initial = Pcs_cvx;
        
        x_ke = 129.4201698;
        y_ke = (5/3)*x_ke;
        ht_ke = ht_e;
        
        x_kcu = x_cu + Qcu;
        y_kcu = y_cu;
        ht_kcu = ht_cu;
        x_kpu = x_pu;
        y_kpu = y_pu - Qpu;
        ht_kpu = ht_pu;
        
        uav_destination = (x_kcu-x_initial).^2 + (y_kcu-y_initial).^2 + (ht_kcu-ht_initial).^2;
        uav_eavesdropper = (x_ke-x_initial).^2 + (y_ke-y_initial).^2 + (ht_ke-ht_initial).^2;
        uav_primary = (x_kpu-x_initial).^2 + (y_kpu-y_initial).^2 + (ht_kpu-ht_initial).^2;
        
        %%%%%%%%%%%%%%%%%% - Original Value - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        rate_d_original = log2(1+exp(-index)*gamma*d_cs_cu*Pcs_initial./(gamma*Pu_initial./uav_destination+1));
        rate_e_original = log2(1+gamma*d_cs_e*Pcs_initial./(gamma*Pu_initial./uav_eavesdropper+1));
        rate_s_original = rate_d_original-rate_e_original;
        Original_Value = (1/N)*sum(rate_s_original(1,1:N));
        
        %%%%%%%%%%%%%%%%% - Secrecy rate average - %%%%%%%%%%%%%%%%%%%%%%%%
        AA = exp(-index)*gamma*d_cs_cu./(gamma*Pu_initial./uav_destination+1);
        BB = gamma*d_cs_e./(gamma*Pu_initial./uav_eavesdropper+1);
        for n = 1:1:length(AA)
            if (AA(n) > BB(n))
                Pcs_initial_1(n) = Pcs_cvx(n);
                Pu_initial_1(n) = Pu_cvx(n);
            else
                Pcs_initial_1(n) = 0;
                Pu_initial_1(n) = 0;
            end
        end
        
        uav_destination = (x_kcu-x_initial).^2 + (y_kcu-y_initial).^2 + (ht_kcu-ht_initial).^2;
        uav_eavesdropper = (x_ke-x_initial).^2 + (y_ke-y_initial).^2 + (ht_ke-ht_initial).^2;
        rate_d = log2(1+exp(-index)*gamma*d_cs_cu*Pcs_initial_1./(gamma*Pu_initial_1./uav_destination+1));
        rate_e = log2(1+gamma*d_cs_e*Pcs_initial_1./(gamma*Pu_initial_1./uav_eavesdropper+1));
        rate_s = rate_d-rate_e;
        SceceryRate_Average = (1/N)*sum(rate_s(1,1:N));
        
        if (k == 8)
            fprintf(fid_SercrecyRate_Threshold, '%.32f \n',SceceryRate_Average);
            fprintf(fid_Iteration_Threshold, '%.32f \n',i);
        end
    
        figure(2)
        hold on
        subplot(3,1,1)
        plot(a(1,1:N),Pu_initial(1,1:N));
        title('Power of the UAV')
        xlabel('slot')
        ylabel('Power')
    
        %{
        subplot(3,1,2)
        plot(a(1,1:N),((x_initial-x_k).^2 + (y_initial-y_e).^2 + (ht_initial-ht_e).^2).^(1/2));
        title('distance between the UAV an eavesdropper')
        xlabel('slot')
        ylabel('Distance')
        %}
        
        subplot(3,1,3)
        plot(a(1,1:N),rate_s);%d_cvx-e_cvx
        title('Rate')
        xlabel('slot')
        ylabel('Rate')
    
        disp('Approximated value');
        disp(sum(d_cvx(1,1:N)-e_cvx(1,1:N))/N); %disp(double(Objective))
        disp('Original value');
        disp(Original_Value);
        disp('Secrecy rate average');
        disp(SceceryRate_Average);
        
        disp('Iteration k');
        disp(k);
        disp('Iteration i');
        disp(i);
        %l = l + 1;
    
        R_cvx = d_cvx-e_cvx;
    
        disp('Step');
        disp(sum(R_cvx(1,1:N))/N-sum(R_cvx_update(1,1:N))/N);
    
        disp('Interfere Power');
        disp((1/N)*sum(beta*d_cs_pu*Pcs_cvx(1,1:N)+c_cvx(1,1:N)));
        %{
        if (abs(sum(R_cvx(1,1:N))/N-sum(R_cvx_update(1,1:N))/N)<=10^-6)
            break;
        end
        %}
        R_cvx_update = R_cvx;
    end
    
    %%%%%%%%%%%%%%%%%%%%%% - Scerecy Rate Output - %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid_SecrecyRate, '%.32f \n',SceceryRate_Average);
    fprintf(fid_Time, '%.32f \n',T_flight);
    
    %%%%%%%%%%%%%%%%%%%%%% - Trajectory Output - %%%%%%%%%%%%%%%%%%%%%%%%%%
    if (k == 1)
        fprintf(fid_x_1, '%.32f \n',x_initial);
        fprintf(fid_y_1, '%.32f \n',y_initial);
        fprintf(fid_ht_1, '%.32f \n',ht_initial);
    end
    
    if (k == 2)
        fprintf(fid_x_2, '%.32f \n',x_initial);
        fprintf(fid_y_2, '%.32f \n',y_initial);
        fprintf(fid_ht_2, '%.32f \n',ht_initial);
    end
    
    if (k == 3)
        fprintf(fid_x_3, '%.32f \n',x_initial);
        fprintf(fid_y_3, '%.32f \n',y_initial);
        fprintf(fid_ht_3, '%.32f \n',ht_initial);
    end
    
    if (k == 4)
        fprintf(fid_x_4, '%.32f \n',x_initial);
        fprintf(fid_y_4, '%.32f \n',y_initial);
        fprintf(fid_ht_4, '%.32f \n',ht_initial);
    end
    
    if (k == 5)
        fprintf(fid_x_5, '%.32f \n',x_initial);
        fprintf(fid_y_5, '%.32f \n',y_initial);
        fprintf(fid_ht_5, '%.32f \n',ht_initial);
    end
    
    if (k == 6)
        fprintf(fid_x_6, '%.32f \n',x_initial);
        fprintf(fid_y_6, '%.32f \n',y_initial);
        fprintf(fid_ht_6, '%.32f \n',ht_initial);
    end
    
    if (k == 7)
        fprintf(fid_x_7, '%.32f \n',x_initial);
        fprintf(fid_y_7, '%.32f \n',y_initial);
        fprintf(fid_ht_7, '%.32f \n',ht_initial);
    end
    
    if (k == 8)
        fprintf(fid_x_8, '%.32f \n',x_initial);
        fprintf(fid_y_8, '%.32f \n',y_initial);
        fprintf(fid_ht_8, '%.32f \n',ht_initial);
    end
    
    if (k == 1)
        fprintf(fid_PowerSource_1, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_1, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_1, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_1, '%.32f \n',a);
        
        fprintf(fid_Rate_1, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_1, '%.32f \n',a);
        
        %fprintf(fid_Distance_1, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_1, '%.32f \n',a);
    end
    
    if (k == 2)
        fprintf(fid_PowerSource_2, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_2, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_2, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_2, '%.32f \n',a);
        
        fprintf(fid_Rate_2, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_2, '%.32f \n',a);
        
        %fprintf(fid_Distance_2, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_2, '%.32f \n',a);
    end
    
    if (k == 3)
        fprintf(fid_PowerSource_3, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_3, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_3, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_3, '%.32f \n',a);
        
        fprintf(fid_Rate_3, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_3, '%.32f \n',a);
        
        %fprintf(fid_Distance_3, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_3, '%.32f \n',a);
    end
    
    if (k == 4)
        fprintf(fid_PowerSource_4, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_4, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_4, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_4, '%.32f \n',a);
        
        fprintf(fid_Rate_4, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_4, '%.32f \n',a);
        
        %fprintf(fid_Distance_4, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_4, '%.32f \n',a);
    end
    
    if (k == 5)
        fprintf(fid_PowerSource_5, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_5, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_5, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_5, '%.32f \n',a);
        
        fprintf(fid_Rate_5, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_5, '%.32f \n',a);
        
        %fprintf(fid_Distance_5, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_5, '%.32f \n',a);
    end
    
    if (k == 6)
        fprintf(fid_PowerSource_6, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_6, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_6, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_6, '%.32f \n',a);
        
        fprintf(fid_Rate_6, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_6, '%.32f \n',a);
        
        %fprintf(fid_Distance_6, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_6, '%.32f \n',a);
    end
    
    if (k == 7)
        fprintf(fid_PowerSource_7, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_7, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_7, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_7, '%.32f \n',a);
        
        fprintf(fid_Rate_7, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_7, '%.32f \n',a);
        
        %fprintf(fid_Distance_7, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_7, '%.32f \n',a);
    end
    
    if (k == 8)
        fprintf(fid_PowerSource_8, '%.32f \n',Pcs_initial_1);
        fprintf(fid_TimeSlot1_8, '%.32f \n',a);
        
        fprintf(fid_PowerUAV_8, '%.32f \n',Pu_initial_1);
        fprintf(fid_TimeSlot2_8, '%.32f \n',a);
        
        fprintf(fid_Rate_8, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot3_8, '%.32f \n',a);
        
        %fprintf(fid_Distance_8, '%.32f \n',((x_initial-x_k).^2+(y_initial-y_k).^2+(ht_initial-ht_k).^2).^(1/2));
        %fprintf(fid_TimeSlot4_8, '%.32f \n',a);
    end
    
    k = k+1;
    
    if (T_flight >= 500)
        break;
    end
    T_flight = T_flight+60;
    %N = 2*T_flight;
end











