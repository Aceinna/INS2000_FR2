clear;
clc;

datadir = '../../data/294/';

%% 

% ===== PVA
fname_gnss = '../femtomes_ins330-2020_10_23_10_15_58-gnss.txt';
fname_sol =  '../femtomes_ins330-2020_10_23_10_15_58-ins.txt';
% fname_fil = '../femtomes_ins330-2020_10_19_09_48_26-process';
fname_ref =  '../novatel_CPT7-2020_10_23_10_15_46-ins.txt'; % NovAtel 
% fname_imu = '../femtomes_ins330-2020_10_19_09_48_26-imu.txt';
% t_process = [460024,463192];
% t_cut = [460024,463192];

% ====== dual antenna heading
fname_ant_heading = '../femtomes_ins330-2020_10_23_10_15_58-heading.txt';
fname_ant_heading_ref =  '../novatel_CPT7-2020_10_23_10_15_46-heading2.txt';


%% ---------------------------Plot raw imu data-------------------------
% rawimu = load(fname_imu);
% % rawodo = load(fname_odo);
% [N,~]=size(rawimu);
% t_shift=floor(rawimu(1,2)/1000)*1000;     % to aviod to display the large GPS time.
% 
% figure,
% ax1(1)=subplot(411); plot(rawimu(1:N,2)-t_shift,rawimu(1:N,4:6)), ylabel('Acc (m/s^2)'), 
%     grid on; legend('Forward','Right','Down'),
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% ax1(2)=subplot(412); plot(rawimu(1:N,2)-t_shift,rawimu(1:N,7:9)), ylabel('Gyro (deg/s)'), 
%     grid on; legend('Forward','Right','Down'),
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% 
% ax1(3)=subplot(413); plot(rawimu(2:N,2)-t_shift,diff(rawimu(1:N,2))), ylabel('IMU Time Diff (s)'), 
% grid on;
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% ax1(4)=subplot(414); plot(rawodo(:,2)-t_shift,rawodo(:,4)), ylabel('IMU Time Diff (s)'), 
% grid on;
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% linkaxes(ax1,'x');
% clear rawimu;


%% 


rawgnss = load(fname_gnss);
sol_X = load(fname_sol);
sol_0 = sol_X(:,2:end);
% fil_X = load(fname_fil);
% fil_0 = fil_X(:,2:end);

true_0 = load(fname_ref);

ant_heading_x = load(fname_ant_heading);
ant_heading_ref = load(fname_ant_heading_ref);

% t_process = [sol_X(1,2), sol_X(end, 2)];
t_process = [440313, 441003];
t_start = t_process(1);
t_end = t_process(2);


%%
% figure,
% 
% ax2(1)=subplot(211);
% plot(rawodo(:,2)-t_shift,rawodo(:,4)),
% hold on ,plot(sol_0(:,1)-t_shift,sqrt(sol_0(:,5).*sol_0(:,5)+sol_0(:,6).*sol_0(:,6)+sol_0(:,7).*sol_0(:,7)));
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% legend('ODO','INS');
% ax2(2)=subplot(212);
% plot(sol_0(:,1)-t_shift,sol_0(:,10)*180/pi),
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% linkaxes(ax2,'x');

%%
WGS84=GetWGS84;
M=RC_Meridian(WGS84.a, WGS84.e2, rawgnss(1,3)*pi/180);
N=RC_PrimeVertical(WGS84.a, WGS84.e2, rawgnss(1,3)*pi/180);



rawgnss(:,3:4) = rawgnss(:,3:4)*pi/180 ;
t_shift=floor(rawgnss(1,2)/10000)*10000;     % to aviod to display the large GPS time.
% 
% figure,
% ax(1) = subplot(211),plot(rawgnss(2:end,2)-t_shift,diff(rawgnss(:,2))), ylabel('GNSS Time Diff (s)'),
%     grid on; 
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% ax(2) = subplot(212),plot(rawgnss(2:end,2)-t_shift,rawgnss(2:end,9)), ylabel('GNSS POSTION TYPE (s)'),
%     grid on; 
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end

%N_gps = rawgnss(:,2)>=t_start & rawgnss(:,2)<=t_end;
gps=rawgnss(:,2:end);

%%

N_true = true_0(:,2)>=t_start & true_0(:,2)<=t_end;
true_1 = true_0(N_true,2:end);
true_1(:,2:3) = true_1(:,2:3) * pi/180;
true_1(:,8:10) = true_1(:,8:10) * pi/180;

sol_0(:,2:3) = sol_0(:,2:3) * pi/180;
sol_0(:,8:10) = sol_0(:,8:10) * pi/180;

[N_U,N_2] =unique(gps(:,1));
gps_2 =gps(N_2 , :); 
gps_1 = interp1(gps_2(:,1),gps_2,true_1(:,1));

[v0, idx0] = unique(sol_0(:, 1));
sol_1 = interp1(sol_0(idx0,1),sol_0(idx0, :),true_1(:,1),'pchip');
sol_1(:,10)=interp1_Azimuth(sol_0(idx0,1),sol_0(idx0,10)-2*pi,true_1(:,1));
% fil= interp1(fil_0(:,1),fil_0,true_1(:,1));

sol = sol_1;
true = true_1;

% ant heading error
h_active_row = ant_heading_ref(:,2)>=t_start & ant_heading_ref(:,2)<=t_end;
ant_heading_ref_1 = ant_heading_ref(h_active_row,[2 6]);
[a0, aidx0] = unique(ant_heading_x(:, 2));
ant_heading_x_1 = interp1(ant_heading_x(aidx0,2),ant_heading_x(aidx0,[2 6]),ant_heading_ref_1(:,1),'pchip');
 
error_h = abs(ant_heading_x_1(:,2)-ant_heading_ref_1(:,2)) - 180;
t_disp2 = ant_heading_ref_1(:,1)-t_shift;

%% apply lever arm
offset = [0.28;-1.43;-1];

for i = 1:length(true_1)
    att_aa = [true_1(i,8), true_1(i,9), true_1(i,10)];
    C_vn = euler2dcm(true_1(i,8), true_1(i,9), true_1(i,10));
    tmp = true_1(i,2:4)';
    la_r = C_vn * offset;
    d_l(1) = la_r(2) / (M);
    d_l(2) = la_r(1) / (N * cos(true_1(i,2))) ;
    d_l(3) = la_r(3);
    for j = 2:4
        true(i,j) = true_1(i,j) + d_l(j-1);
    end
%         v_ins = [true_0(i,5); true_0(i,6); true_0(i,7)];
%         v_gnss = C_bn*v_ins;
%         for j =5:7
%         true_1(i,j) = v_gnss(j-4);
%         end
%     C_vn = C_bn*C_vb;
%     [att(1),att(2),att(3)] = dcm2euler(C_vn);
%     for j = 8:10
%   %  true(i,j) = att(j-7);
%     end       
end

%% plot trajectory

r_N_gps=(gps(:,2)-gps(1,2))*(M+gps(1,4));
r_E_gps=(gps(:,3)-gps(1,3))*(N+gps(1,4))*cos(gps(1,2));

r_N_gps_1=(gps_1(:,2)-gps(1,2))*(M+gps(1,4));
r_E_gps_1=(gps_1(:,3)-gps(1,3))*(N+gps(1,4))*cos(gps(1,2));

r_N_sol=(sol(:,2)-gps(1,2))*(M+gps(1,4));
r_E_sol=(sol(:,3)-gps(1,3))*(N+gps(1,4))*cos(gps(1,2));

r_N_true=(true(:,2)-gps(1,2))*(M+gps(1,4));
r_E_true=(true(:,3)-gps(1,3))*(N+gps(1,4))*cos(gps(1,2));

% [N_gps_nub,~]=size(gps_1);
% 
% %ͨ�^GPSλ�����GPS�ٶȺͺ���
% for i =2:N_gps_nub 
%     v_gnss(i,1) = ( r_N_gps_1(i) - r_N_gps_1(i-1) )/(gps_1(i,1)-gps_1(i-1,1));
%     v_gnss(i,2) = ( r_E_gps_1(i) - r_E_gps_1(i-1) )/(gps_1(i,1)-gps_1(i-1,1));
%     v_gnss(i,3) = -( gps_1(i,4) - gps_1(i-1,4) )/(gps_1(i,1)-gps_1(i-1,1));
% end
% for i =2:N_gps_nub 
%     att_gnss(i,3) = atan2(r_E_gps_1(i) - r_E_gps_1(i-1), r_N_gps_1(i) - r_N_gps_1(i-1));
% end
% v_gnss(1,1) = v_gnss(2,1);
% v_gnss(1,2) = v_gnss(2,2);
% v_gnss(1,3) = v_gnss(2,3);
% att_gnss(1,3) = att_gnss(2,3);
% att_gnss(:,1:2) = 0;
% error_v=sol(:,5:7)-v_gnss(:,1:3);
% 
% az_limit=pi;
% azimuth_error=sol(:,10)-(att_gnss(:,3));
% I1=find(azimuth_error<-(2*pi-az_limit));
% azimuth_error(I1)=azimuth_error(I1)+2*pi;
% I2=find(azimuth_error>+(2*pi-az_limit));
% azimuth_error(I2)=azimuth_error(I2)-2*pi;
% error_A=[sol(:,8:9)-att_gnss(:,1:2) azimuth_error] *180/pi;
% 
% %������?
t_disp=true(:,1)-t_shift;
t_gps=gps(:,1)-t_shift;
t_sol = sol(:,1)-t_shift;

figure, 
plot3( r_E_true, r_N_true, t_disp,'g*-', 'markersize',5,'linewidth', 0.5), xlabel('East (m)'), ylabel('North (m)'),grid on;  view(0,90);
hold on,plot3( r_E_gps, r_N_gps, t_gps,'m*-', 'markersize',5,'linewidth', 0.5),
hold on,plot3( r_E_sol, r_N_sol, t_sol,'k*-', 'markersize',5,'linewidth', 0.5),
hold on,plot3( r_E_gps(1), r_N_gps(1), t_gps(1), 'r*', 'markersize',20),
legend('True','GNSS','GNSS/INS','Start');


error_r=[r_N_sol - r_N_true, r_E_sol - r_E_true, sol(:,4)-true(:,4)];
error_r1=[r_N_sol - r_N_gps_1, r_E_sol - r_E_gps_1, sol(:,4)-gps_1(:,4)];
error_r2=[r_N_true - r_N_gps_1, r_E_true - r_E_gps_1, true(:,4)-gps_1(:,4)];

t_cut = [sol(1,1), sol(end,1)];

true(:,7) = -true(:,7);
error_v = sol(:,5:7) - true(:,5:7);
az_limit = pi;
azimuth_error=sol(:,10)-(true(:,10));
I1=find(azimuth_error<-(2*pi-az_limit));
azimuth_error(I1)=azimuth_error(I1)+2*pi;
I2=find(azimuth_error>+(2*pi-az_limit));
azimuth_error(I2) = azimuth_error(I2)-2*pi;
error_A = [sol(:,8:9) - true(:,8:9) azimuth_error] *180/pi;

nn_cut=find(true(:,1)>t_cut(1) & true(:,1)<t_cut(2));

sum_error=[ sum_statistic(error_r(nn_cut,1));sum_statistic(error_r(nn_cut,2));sum_statistic(error_r(nn_cut,3));...
            sum_statistic(error_v(nn_cut,1));sum_statistic(error_v(nn_cut,2));sum_statistic(error_v(nn_cut,3));...
            sum_statistic(error_A(nn_cut,1));sum_statistic(error_A(nn_cut,2));sum_statistic(error_A(nn_cut,3))]',

 sum_error1=[ sum_statistic(error_r2(nn_cut,1));sum_statistic(error_r2(nn_cut,2));sum_statistic(error_r2(nn_cut,3))]';

%% plot PVA difference

figure,
ax(1) = subplot(411),plot(t_disp,error_r,'.-'),
ylabel('pos (m)'), grid on;
temp = axis;
hold on ,
% plot(t_disp,fil(:,2:4)+temp(3),'markersize',30),
hold off,
legend('North','East','Height','North STD','East STD','Height STD','Location','EastOutside'), 
if t_shift == 0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t_shift) ' (sec)']);
end
title(gca, 'FR2 vs CPT7 Real time solution')

ax(2) = subplot(412),plot(t_disp,error_v, '.-'),
ylabel('vel (m/s)'), grid on;
temp =axis;
hold on ,
% plot(t_disp,fil(:,5:7)+temp(3),'markersize',30), 
hold off,
legend('North','East','down','North','East','down','Location','EastOutside'),
if t_shift==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t_shift) ' (sec)']);
end
ax(3) = subplot(413),plot(t_disp,error_A, '.-'),
ylabel('att (deg)'), grid on;
temp =axis;
hold on ,
% plot(t_disp,fil(:,8:10)+temp(3),'markersize',30), 
hold off,
legend('Roll','Pitch','Yaw','Roll STD','Pitch STD','Yaw STD' ,'Location','EastOutside'), 
if t_shift==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t_shift) ' (sec)']);
end
% plot ant heading error
ax(4) = subplot(414),plot(t_disp2,error_h, '.-'),
ylabel('ant heading (deg)'), grid on;
temp =axis;
hold on ,
% plot(t_disp,fil(:,8:10)+temp(3),'markersize',30), 
hold off,
legend('AntH','Pitch','Yaw','Roll STD','Pitch STD','Yaw STD' ,'Location','EastOutside'), 
linkaxes(ax,'x');
if t_shift==0
    xlabel(' Time (sec)');
else
    xlabel(['Time -' int2str(t_shift) ' (sec)']);
end

% 
% %%
% %�L�u���Sλ�á��ٶȡ��ˑB�����ٶ�Ӌ��ƫ��������ƫ
% figure,
% ax(1) = subplot(311),plot(t_disp,[r_E_sol r_N_sol sol(:,4)]),
% ylabel('Position (m)'), grid on;
% legend('North','East','Height','Location','EastOutside'), 
% if t_shift == 0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% 
% ax(2) = subplot(312),plot(t_disp,sol(:,5:7)),
% ylabel('V (m/s)'), grid on;
% legend('North','East','down','Location','EastOutside'), 
% if t_shift == 0    xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% ax(3) = subplot(313),plot(t_disp,sol(:,8:10)*180/pi),
% ylabel('Attitude  (deg)'), grid on;
% legend('Roll','Pitch','Yaw','Location','EastOutside'), 
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% 

%%
% figure, %draw the sensor biases estimations
% ax(1) = subplot(211),
% plot(t_disp,sol(:,11:13)/3600),ylabel('Gyro bias (deg/s)'), legend('x','y','z'), grid on; 
% temp =axis;
% hold on ,plot(t_disp,fil(:,11:13)/3600+temp(3),':','markersize',2), hold off,
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% ax(2) = subplot(212),
% plot(t_disp,sol(:,14:16)/100000),ylabel('Acc bias (m/s^2)'), legend('x','y','z'), grid on; 
% temp =axis;
% hold on ,plot(t_disp,fil(:,14:16)/100000+temp(3),':','markersize',2), hold off,
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
% % %---------------�L�uIMUԭʼ�������Q----------
% % N_imu = rawimu_0(:,1)>=t_start&rawimu_0(:,1)<=t_end;
% % rawimu = rawimu_0(N_imu,:);
% % [N,~]=size(rawimu);
% % 
% % figure,
% % ax(1)=subplot(211); plot(rawimu(1:N,1)-t_shift,rawimu(1:N,3:5)), ylabel('Acc (m/s^2)'), 
% %     grid on; legend('Forward','Right','Down'),
% % if t_shift==0
% %     xlabel(' Time (sec)');
% % else
% %     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% % end
% % ax(2)=subplot(212); plot(rawimu(1:N,1)-t_shift,rawimu(1:N,6:8)), ylabel('Gyro (deg/s)'), 
% %     grid on; legend('Forward','Right','Down'),
% % if t_shift==0
% %     xlabel(' Time (sec)');
% % else
% %     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% % end
% % 
% % %Ӌ���D�ӽǶ�
% % angle = mean(rawimu(:,6:8))*(t_end - t_start);
% % 
% % %�L�uIMU?G��
% % figure,
% % plot(rawimu(2:N,1)-t_shift,[diff(rawimu(:,1)) diff(rawimu(:,1))]), ylabel('IMU Time Diff (s)'), 
% %     grid on; legend('Systemtime','GNSS time'),
% % if t_shift==0
% %     xlabel(' Time (sec)');
% % else
% %     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% % end
% clear error_horizontal;
% error_horizontal = sqrt(error_r(:,1).^2 + error_r(:,2).^2);
% hor_err_all = error_horizontal;
%     clear error_horizontal1;
%     clear error_Vertical1;
%     num=length(hor_err_all);
%     hor_err_all=sort(hor_err_all);
%     
%     err_50p =hor_err_all(round(0.50*num));
%     err_68p =hor_err_all(round(0.68*num));
%     err_95p =hor_err_all(round(0.95*num));
%     err_99p =hor_err_all(round(0.99*num));
%     
%     rtk_err=[err_50p err_68p err_95p err_99p]
%     id =1;
%    
%     figure
%     cdfplot(hor_err_all);
%     hold on
%     set(gca,'FontWeight','bold','FontSize',12);
%     xlim(gca, [0 25]);
%     xlabel(gca, 'Horizontal error (m)')
%     ylabel(gca, 'CDF')
% %     title('st teseoV all drives');
%     title('OPEN RTK');
% %   legend('st teseoV-all drives');
% %     title('m8 psb all drives');
% %     legend('st teseoV-all drives');
%     box on;
%     
%     %�L�u�`�����Q
% figure,
% plot(t_disp,error_horizontal),
% ylabel('Horizontal error (m)'), grid on;
% if t_shift==0
%     xlabel(' Time (sec)');
% else
%     xlabel(['Time -' int2str(t_shift) ' (sec)']);
% end
%     
%     figure
%     hold on
%     set(gca,'xticklabel',{'','CEP50','CEP68','CEP95','CEP99',''});
%     ylim(gca, [0 10]);
%     set(gca,'FontWeight','bold','FontSize',12);
%     ylabel(gca, 'Horizontal error (m)')
%     bc=bar(rtk_err(id,:)');
%     set(bc,'facecolor',[0 0.5 0]);
%     title('OPEN RTK');
% %     legend('st teseoV-all drives');
%     grid on;
%    
%     for i=1:4
%     if (rtk_err(id,i)<1)
%     text(i,rtk_err(id,i),num2str(rtk_err(id,i),2),...
%     'FontWeight','bold','HorizontalAlignment','center',...
%      'VerticalAlignment','bottom')
%     else
%     text(i,rtk_err(id,i),num2str(rtk_err(id,i),3),...
%     'FontWeight','bold','HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
%     end
%     end 
%        %�߳����?
%     error_vertical = abs(error_r(:,3));
%     vertical_err_all1 = error_vertical;
%     clear error_horizontal1;
%     clear error_Vertical1;
%     num=length(vertical_err_all1);
%     vertical_err_all=sort(vertical_err_all1);
%     
%     err_50p =vertical_err_all(round(0.50*num));
%     err_68p =vertical_err_all(round(0.68*num));
%     err_95p =vertical_err_all(round(0.95*num));
%     err_99p =vertical_err_all(round(0.99*num));
%     
%     rtk_err_1=[err_50p err_68p err_95p err_99p]
%     id =1;
%    
%     figure
%     cdfplot(vertical_err_all);
%     hold on
%     set(gca,'FontWeight','bold','FontSize',12);
%     xlim(gca, [0 20]);
%     xlabel(gca, 'Horizontal error (m)')
%     ylabel(gca, 'CDF')
% %     title('st teseoV all drives');
%     title('OPEN RTK');
% %   legend('st teseoV-all drives');
% %     title('m8 psb all drives');
% %     legend('st teseoV-all drives');
%     box on;
%     
%     figure
%     hold on
%     set(gca,'xticklabel',{'','CEP50','CEP68','CEP95','CEP99',''});
%     ylim(gca, [0 20]);
%     set(gca,'FontWeight','bold','FontSize',12);
%     ylabel(gca, 'Vertical error (m)')
%     bc=bar(rtk_err_1(id,:)');
%     set(bc,'facecolor',[0 0.5 0]);
%     title('OPEN RTK');
% %     legend('st teseoV-all drives');
%     grid on;
%    
%     for i=1:4
%     if (rtk_err_1(id,i)<1)
%     text(i,rtk_err_1(id,i),num2str(rtk_err_1(id,i),2),...
%     'FontWeight','bold','HorizontalAlignment','center',...
%      'VerticalAlignment','bottom')
%     else
%     text(i,rtk_err_1(id,i),num2str(rtk_err_1(id,i),3),...
%     'FontWeight','bold','HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
%     end
%     end  
%     %��̬���?
%     error_A(:,1) = error_A(:,1) - mean(error_A(:,1));
%     error_A(:,2) = error_A(:,2) - mean(error_A(:,2));
%     clear error_horizontal1;
%     clear error_Vertical1;
%     hor_err_all1 = abs(error_A(:,1));
%     hor_err_all2 = abs(error_A(:,2));
%     hor_err_all3 = abs(error_A(:,3));
%    
%     num=length(hor_err_all1);
%     hor_err_all1=sort(hor_err_all1);
%     
%     err_50p =hor_err_all1(round(0.50*num));
%     err_68p =hor_err_all1(round(0.68*num));
%     err_95p =hor_err_all1(round(0.95*num));
%     err_99p =hor_err_all1(round(0.99*num));
%     
%     roll_error=[err_50p err_68p err_95p err_99p]
%     
%     
%         hor_err_all2=sort(hor_err_all2);
%     
%     err_50p =hor_err_all2(round(0.50*num));
%     err_68p =hor_err_all2(round(0.68*num));
%     err_95p =hor_err_all2(round(0.95*num));
%     err_99p =hor_err_all2(round(0.99*num));
%     
%     pitch_error=[err_50p err_68p err_95p err_99p]
%     
%    hor_err_all3=sort(hor_err_all3);
%     
%     err_50p =hor_err_all3(round(0.50*num));
%     err_68p =hor_err_all3(round(0.68*num));
%     err_95p =hor_err_all3(round(0.95*num));
%     err_99p =hor_err_all3(round(0.99*num));
%     
%      heading_error=[err_50p err_68p err_95p err_99p]
% 
%     id =1;
%    
%     figure
%     cdfplot(hor_err_all1),
%     hold on,cdfplot(hor_err_all2),
%     hold on,cdfplot(hor_err_all3),
% 
%     hold on,
%     set(gca,'FontWeight','bold','FontSize',12);
%     legend('roll','pitch','heading'), 
%     xlim(gca, [0 5]);
%     xlabel(gca, 'Attitude error (deg)')
%     ylabel(gca, 'CDF')
% %     title('st teseoV all drives');
%     title('OPEN RTK');
% %   legend('st teseoV-all drives');
% %     title('m8 psb all drives');
% %     legend('st teseoV-all drives');
%     box on;
%      

% fclose all;


