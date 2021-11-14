clear all, close all
clc
%yalmip('clear');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - Parameters - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_cs = 0; % cognitive source location - x
y_cs = 0; % cognitive source location - y
x_cu = 300; % cognitive user location -x
y_cu = 0; % cognitive user location -y
x_e = 150; % eavesdropper location -x
y_e = 250; % evaesdropper location -y
x_pu = 300-300; % primary user location -x
y_pu = 250; % primary user location -y
Pcs_average = 10^(40/10-3)/2; % the average transmit power
Pcs_max = 10^(40/10-3); % the peak transmit power

l = 0; % iteration index
noise = 10^(-70/10-3); % the power spectral density of the AWGN -104
beta = 10^(10/10-3); % the channel power gain at the reference distance d = 1m -14
gamma = 10^8; % gamma = beta/noise = 90dB
phi = 3; % the path loss exponent
index = 0.5772156;
d_cs_cu = (sqrt((x_cs-x_cu)^2+(y_cs-y_cu)^2))^(-phi); % the distance between the cognitive source and the cognitive user
d_cs_e = (sqrt((x_cs-x_e)^2+(y_cs-y_e)^2))^(-phi); % the distance between the cognitive source and cognitive eavesdropper
d_cs_pu = (sqrt((x_cs-x_pu)^2+(y_cs-y_pu)^2))^(-phi); % the distance between the cognitive source and the primary user

%%%%%%%%%%%%%%%%%%%%%%%%% - SecrecyRate_Iteration - %%%%%%%%%%%%%%%%%%%%%%%
fid_SercrecyRate_Threshold = fopen('..\UAV\Output\Threshold\SecrecyRate_Iteration\No_Jamming\Secrecy_Rate.txt','w');
fid_Iteration_Threshold = fopen('..\UAV\Output\Threshold\SecrecyRate_Iteration\No_Jamming\Iteration.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%% - SecrecyRate_Time - %%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_SecrecyRate = fopen('..\UAV\Output\SecrecyRate_Time\No_Jamming\Secrecy_Rate.txt','w');
fid_Time = fopen('..\UAV\Output\SecrecyRate_Time\No_Jamming\Time.txt','w');

fid_PowerSource_1 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_1.txt','w');
fid_TimeSlot1_1 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot1.txt','w');
fid_PowerSource_2 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_2.txt','w');
fid_TimeSlot1_2 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot2.txt','w');
fid_PowerSource_3 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_3.txt','w');
fid_TimeSlot1_3 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot3.txt','w');
fid_PowerSource_4 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_4.txt','w');
fid_TimeSlot1_4 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot4.txt','w');
fid_PowerSource_5 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_5.txt','w');
fid_TimeSlot1_5 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot5.txt','w');
fid_PowerSource_6 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_6.txt','w');
fid_TimeSlot1_6 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot6.txt','w');
fid_PowerSource_7 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_7.txt','w');
fid_TimeSlot1_7 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot7.txt','w');
fid_PowerSource_8 = fopen('..\UAV\Output\Power_Source\No_Jamming\Power_Source_8.txt','w');
fid_TimeSlot1_8 = fopen('..\UAV\Output\Power_Source\No_Jamming\Time_Slot8.txt','w');

fid_Rate_1 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_1.txt','w');
fid_TimeSlot2_1 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot1.txt','w');
fid_Rate_2 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_2.txt','w');
fid_TimeSlot2_2 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot2.txt','w');
fid_Rate_3 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_3.txt','w');
fid_TimeSlot2_3 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot3.txt','w');
fid_Rate_4 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_4.txt','w');
fid_TimeSlot2_4 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot4.txt','w');
fid_Rate_5 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_5.txt','w');
fid_TimeSlot2_5 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot5.txt','w');
fid_Rate_6 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_6.txt','w');
fid_TimeSlot2_6 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot6.txt','w');
fid_Rate_7 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_7.txt','w');
fid_TimeSlot2_7 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot7.txt','w');
fid_Rate_8 = fopen('..\UAV\Output\Rate\No_Jamming\Rate_8.txt','w');
fid_TimeSlot2_8 = fopen('..\UAV\Output\Rate\No_Jamming\Time_Slot8.txt','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 1;
T_flight = 80; % the UAV's trajectory versus period
%N = 2*T_flight; % the number of time slots
N = 500;
while (1)
    Pcs_initial = zeros(1,N);
    for i = 1:1:N
        Pcs_initial(1,i) = Pcs_average/2;%2.01
    end
    % - the initial slack variables - %
    t_0_initial = 1./(exp(-index)*gamma*d_cs_cu*Pcs_initial);
    
    a = zeros(1,N);
    a(1)=1;
    for i = 2:1:N
        a(i)= a(i-1)+1;
    end
    %R_new=double(sum(d(1,1:N)-e(1,1:N)));
    cvx_solver SDPT3
    R_cvx_update = zeros(1,N);
    %while (1)
    for i = 1:1:2
        cvx_begin
            variables Pcs_cvx(1,N) 
            variables d_cvx(1,N) t_0_cvx(1,N) e_cvx(1,N) k_cvx(1,N);    
            maximize(sum(d_cvx(1,1:N)-e_cvx(1,1:N))/N) %d_cvx(1,1:N)-e_cvx(1,1:N)
            subject to    
                %%% - Power constraints - %%%
                (1/N)*sum(Pcs_cvx(1,1:N)) <= Pcs_average;
            
                %%% - Cognitive radio constraint - %%%
                (1/N)*sum(beta*d_cs_pu*Pcs_cvx(1,1:N)) <= 10^(-5);% 5 6 7
            
                for n = 1:1:N
                    %%% - Power constraints - %%%
                    Pcs_cvx(n) >= 0;
                    Pcs_cvx(n) <= Pcs_max;
                
                    %%% - Rd constraints - %%%
                    d_cvx(n) <= log2(1+1/t_0_initial(n))-(t_0_cvx(n)-t_0_initial(n))/(t_0_initial(n)*(t_0_initial(n)+1)*log(2));
                    norm([1 0.5*sqrt(exp(-index)*gamma*d_cs_cu)*(Pcs_cvx(n)-t_0_cvx(n))]) <= 0.5*sqrt(exp(-index)*gamma*d_cs_cu)*(Pcs_cvx(n)+t_0_cvx(n));
                
                    %%% - Re constraints - %%%
                    e_cvx(n) >= log2(1+gamma*d_cs_e*Pcs_initial(n))+gamma*d_cs_e*(Pcs_cvx(n)-Pcs_initial(n))/((1+gamma*d_cs_e*Pcs_initial(n))*log(2));
                
                    %%% - Implicit constraints - %%%
                    Pcs_cvx(n)+t_0_cvx(n) >= 0;
                    t_0_cvx(n) >= 0;   
                end
        cvx_end
        %disp(d_cvx - (log2(1 + 1./t_0_initial) - (t_0_cvx - t_0_initial)./(t_0_initial.*(t_0_initial + 1)*log(2))));
        %disp(Pcs_cvx.*t_0_cvx);
        %disp(e_cvx - (log2(1 + gamma*d_cs_e*Pcs_initial) + gamma*d_cs_e*(Pcs_cvx - Pcs_initial)./((1 + gamma*d_cs_e*Pcs_initial)*log(2))));
        %%% - Update - %%%
        t_0_initial = t_0_cvx;
        %Pcs_initial = Pcs_cvx;
        AA = exp(-index)*gamma*d_cs_cu;
        BB = gamma*d_cs_e;
        for n = 1:1:length(Pcs_cvx)
            if (AA > BB)
                Pcs_initial(n) = Pcs_cvx(n);   
            else
                Pcs_initial(n) = 0;
            end
        end
    
        rate_d = log2(1+exp(-index)*gamma*d_cs_cu*Pcs_initial);
        rate_e = log2(1+gamma*d_cs_e*Pcs_initial);
        rate_s = rate_d-rate_e;
        SecrecyRate_Average = (1/N)*sum(rate_s(1,1:N));
        %fprintf(fid_SercrecyRate, '%.32f \n',rate_average);
        %fprintf(fid_Iteration, '%.32f \n',i);
        if (k == 8)
            fprintf(fid_SercrecyRate_Threshold, '%.32f \n',SecrecyRate_Average);
            fprintf(fid_Iteration_Threshold, '%.32f \n',i);
        end
        
        disp('Approximated value');
        disp(sum(d_cvx(1,1:N)-e_cvx(1,1:N))/N); %disp(double(Objective))
        disp('Secrecy rate average');
        disp(SecrecyRate_Average);
  
        disp('Iteration k');
        disp(k);
        disp('Iteration i');
        disp(i);
        %l = l + 1;
    
        R_cvx = d_cvx-e_cvx;
    
        disp('Step');
        disp(sum(R_cvx(1,1:N))/N-sum(R_cvx_update(1,1:N))/N);
    
        disp('Interfere Power');
        disp((1/N)*sum(beta*d_cs_pu*Pcs_cvx(1,1:N)));
        %{
        if (abs(sum(R_cvx(1,1:N))/N-sum(R_cvx_update(1,1:N))/N) <= 10^-6)
            break;
        end
        %}
        R_cvx_update = R_cvx;
    end
    
    %%%%%%%%%%%%%%%%%%%%%% - Scerecy Rate Output - %%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid_SecrecyRate, '%.32f \n',SecrecyRate_Average);
    fprintf(fid_Time, '%.32f \n',T_flight);
    
    if (k == 1)
        fprintf(fid_PowerSource_1, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_1, '%.32f \n',a);
        
        fprintf(fid_Rate_1, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_1, '%.32f \n',a);
    end
    
    if (k == 2)
        fprintf(fid_PowerSource_2, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_2, '%.32f \n',a);
        
        fprintf(fid_Rate_2, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_2, '%.32f \n',a);
    end
    
    if (k == 3)
        fprintf(fid_PowerSource_3, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_3, '%.32f \n',a);
        
        fprintf(fid_Rate_3, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_3, '%.32f \n',a);
    end
    
    if (k == 4)
        fprintf(fid_PowerSource_4, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_4, '%.32f \n',a);
              
        fprintf(fid_Rate_4, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_4, '%.32f \n',a);
    end
    
    if (k == 5)
        fprintf(fid_PowerSource_5, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_5, '%.32f \n',a);
        
        fprintf(fid_Rate_5, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_5, '%.32f \n',a);
    end
    
    if (k == 6)
        fprintf(fid_PowerSource_6, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_6, '%.32f \n',a);
        
        fprintf(fid_Rate_6, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_6, '%.32f \n',a);
    end
    
    if (k == 7)
        fprintf(fid_PowerSource_7, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_7, '%.32f \n',a);
        
        fprintf(fid_Rate_7, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_7, '%.32f \n',a);
    end
    
    if (k == 8)
        fprintf(fid_PowerSource_8, '%.32f \n',Pcs_initial);
        fprintf(fid_TimeSlot1_8, '%.32f \n',a);
        
        fprintf(fid_Rate_8, '%.32f \n',rate_s);
        fprintf(fid_TimeSlot2_8, '%.32f \n',a);
    end
    
    k = k+1;
    
    if (T_flight >= 500)
        break;
    end
    T_flight = T_flight+60;
    %N = 2*T_flight;
end











