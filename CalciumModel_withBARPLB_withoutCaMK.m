function [scvrmse,AIC,output] = CalciumModel_withBARPLB_withoutCaMK()

par = [];
load parset_withBARPLB_withoutCaMK

%% Parameter Values
% physical constants & geometry 
Cm = 155.4; %pF - from Walden et al.
Vtot = Cm/6.76*1e-12; %L - capacitance to volume ratio from Walden et al.
F = 9.6485e4*1e-3; %C*mmol^-1
R = 8.3145*1e-3; %J*mmol^-1*K^-1
T = 310.15; %K
lambda_myo = 0.65; %dimensionless - Vmyo/Vtotal - from Page et al.
lambda_sr = 0.035; %dimensionless - Vsr/Vtotal - from Page et al.
vds = par(2); %s^-1 - fitted
vlcc = par(3); %s^-1 - fitted
vryr = par(4); %s^-1 - fitted

% fixed ionic concentrations
Nai = 10; %mM - from Luo and Rudy
Nao = 140; %mM - from Luo and Rudy
Cao = 1.8; %mM - from Luo and Rudy

% buffering
Bfluo = 5e-3; %mM - experiment condition from Dibb et al.
Kfluo = 8.64e-4; %mM - from Ito et al.
Bbuff_myo = par(21); %mM - fitted
Kdmyo = par(29); %mM - fitted
Bbuff_sr = par(1); %mM - fitted
Kdsr = par(30); %mM - fitted

% lcc
Vhalf_lcc = -2e-3; %V - from Hinch et al.
deltaV_lcc = 7e-3; %V - from Hinch et al.
kopen_lcc = 1e4; %s^-1 - value chosen to be arbitrarily fast
kclose_lcc = par(7); %s^-1 - fitted
kinactmin_lcc = par(8); %mM^-1*s^-1 - fitted 
kinactmax_lcc = par(9); %mM^-1*s^-1 - fitted
kdeinactmin_lcc = par(10); %s^-1 - fitted
kdeinactmax_lcc = par(11); %s^-1 - fitted

% ryr
Kryr = par(22); %mM - fitted
kopen_ryr = par(12); %s^-1 - fitted
kclose_ryr = par(13); %s^-1 - fitted
kinactmin_ryr = par(14); %s^-1 - fitted
kinactmax_ryr = par(15); %s^-1 - fitted
kdeinactmin_ryr = par(16); %s^-1 - fitted
kdeinactmax_ryr = par(17); %s^-1 - fitted

% serca
Kserca = par(24); %mM - fitted
Jsercao = par(5); %mmol*s^-1*(liter cell vol)^-1 - fitted

% ncx
Jncxo = par(18); %mmol*s^-1*(liter cell vol)^-1 - fitted
Kncxn = 87.5; %mM - from Luo and Rudy
Kncxc = 1.38; %mM - from Luo and Rudy
kncxsat = 0.1; %dimensionless - from Luo and Rudy
eta_ncx = 0.35; %dimensionless - from Luo and Rudy

% pmca
Jpmca_max = par(19); %mmol*s^-1*(liter cell vol)^-1 - fitted
Kpmca = par(23); %mM - fitted

% leak currents
vsrleak = par(6); %s^-1 - fitted
gpmleak = par(20); %mmol*V^-1*s^-1*(liter cell vol)^-1 - fitted

% beta-AR system regulation
%lcc
%phi_lcc_BAR_max = par(26); %dimensionless - fitted
phi_lcc_BAR_max = 0;
Kiso_lcc = 2.8083e-5; %mM
Hiso_lcc = 0.7984; %dimensionless
tauBAR_lcc = 36/3; %s
%zetalcc = par(27) %dimensionless - fitted
zetalcc = 1;
%plb
phi_plb_BAR_max = par(25); %dimensionless - fitted
Kiso_plb = 1.8299e-4/10; %mM - reduced by a factor of 10 to account for 10 fold higher affinity of isoprenaline compared to norepinephrine
Hiso_plb = 0.7785; %dimensionless
tauBAR_plb = 31.0148/3; %s
zetaserca = par(28); %dimensionless - fitted

% Electrophysiology Model 
thalfo_min = 0.0218; %s
thalfo_max = 0.0337; %s
a_max = 0.0034; % s^2
eta_min = 0.6918; %dimensionless
eta_max = 2.0608; %dimensionless
Eo_min = -115.8435e-3; %V
Eo_max = -86.0027e-3; %V
Emax = 30.1721e-3; %V
Kiso_AP = 2.4040e-006; %mM
Hiso_AP = 1.1337; %dimensionless
tauAPD_H = 0.2; %s
tauAPD_BAR = 6; %s

% CaMK
phi_plb_4Hz_CaMK_con = 0;
phi_plb_6Hz_CaMK_con = 0;
phi_plb_8Hz_CaMK_con = 0;

phi_plb_4Hz_CaMK_iso = 0;
phi_plb_6Hz_CaMK_iso = 0;
phi_plb_8Hz_CaMK_iso = 0;

phi_lcc_4Hz_CaMK_con = 0;
phi_lcc_6Hz_CaMK_con = 0;
phi_lcc_8Hz_CaMK_con = 0;

phi_lcc_4Hz_CaMK_iso = 0;
phi_lcc_6Hz_CaMK_iso = 0;
phi_lcc_8Hz_CaMK_iso = 0;

% Refs.
% [1] Dibb KM, et al. (2007) J Physiol. 585:579-592
% [2] Hartman JM, et al. (2010) Am J Physiol Heart Cir Physiol. 299:H1996-H2008
% [3] Hinch, et al. (2004) Biophys J. 87:3723-3736
% [4] Ito K, et al. (2000) Circ Res. 87:588-595
% [5] Luo CH and Rudy Y. (1994) Circ Res. 74:1071-1096
% [6] Page E, et al. (1971) Proc Nat Acad Sci. 68:1465-1466
% [7] Walden AP, et al. (2009) J of Mol and Cell Cardiology. 46:463-473
% [8] Yatani A and Brown AM. (1989) Science. 245:71-74

%% Simulate current clamp protocol
nstates = 6;
y0 = [1; zeros(nstates-1,1); 1e-6; 0.1; 0; thalfo_min; 0; 0]; %initial condition for odes
t0 = [0; 120; 120*2; 120*3; 120*4; 120*5; 120*6; 120*7; 120*8; 120*9];
tf = [120; 120*2; 120*3; 120*4; 120*5; 120*6; 120*7; 120*8; 120*9; 120*10];
HRstep = [4; 6; 4; 8; 4; 4; 6; 4; 8; 4];
ciso_step = [0, 0, 0, 0, 0, 30e-6, 30e-6, 30e-6, 30e-6, 30e-6];
phi_lcc_CaMK_step = [phi_lcc_4Hz_CaMK_con, phi_lcc_6Hz_CaMK_con, phi_lcc_4Hz_CaMK_con, phi_lcc_8Hz_CaMK_con, phi_lcc_4Hz_CaMK_con, ...
                phi_lcc_4Hz_CaMK_iso, phi_lcc_6Hz_CaMK_iso, phi_lcc_4Hz_CaMK_iso, phi_lcc_8Hz_CaMK_iso, phi_lcc_4Hz_CaMK_iso];
phi_plb_CaMK_step = [phi_plb_4Hz_CaMK_con, phi_plb_6Hz_CaMK_con, phi_plb_4Hz_CaMK_con, phi_plb_8Hz_CaMK_con, phi_plb_4Hz_CaMK_con, ...
                phi_plb_4Hz_CaMK_iso, phi_plb_6Hz_CaMK_iso, phi_plb_4Hz_CaMK_iso, phi_plb_8Hz_CaMK_iso, phi_plb_4Hz_CaMK_iso];

tout = 0;
yout = y0;
options = odeset('InitialStep',1e-6);

out_cmyo = cell(1,10);
out_csr = cell(1,10);
out_theta = cell(1,10);
out_thalf = cell(1,10);
out_Ia = cell(1,10);
out_Ilcc = cell(1,10);
out_Incx = cell(1,10);
out_Ipmca = cell(1,10);
out_Ileak = cell(1,10);
out_Jds = cell(1,10);
out_Jserca = cell(1,10);
out_Jncx = cell(1,10);
out_Jpmca = cell(1,10);
out_Jpmleak = cell(1,10);
out_Jsrleak = cell(1,10);
out_Jryr = cell(1,10);
out_Em = cell(1,10);
tstep = cell(1,10);

for i=1:length(t0)
    HR = HRstep(i);
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    while tout(end) < tf(i) - 1e-2 %note: the term 1e-2 is an arbitrary 
                                   %small number to account for rounding errors
        if tout(end) - 1/HR >= t0(i)
            ix = find(tout > tout(end) - 1/HR,1,'first');
            options = odeset('InitialStep',tout(ix) - tout(ix-1));
        end
        sol = ode15s(@DYDT,[tout(end),tout(end) + 1/HR],y0,options);
        if sol.x(end) ~= tout(end) + 1/HR % solver failed
            return
        end
        tout = [tout, sol.x(2:end)];
        yout = [yout, sol.y(:,2:end)]; 
        y0 = sol.y(:,end);
    end
end

output = yout(:,end);

Em_model = zeros(1,length(tout));
for j=1:length(tout)
    ix = find(t0 <= tout(j), 1, 'last');
    HR = HRstep(ix);
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    [~,out2] = DYDT(tout(j),yout(:,j));
    Em_model(j) = out2(13);
end

%% Load Outputs
for i=1:length(t0);
    HR = HRstep(i);
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    ixi = tout >= tf(i) - 1/HR & tout < tf(i);
    youti = yout(:,ixi);
    tstep{i} = tout(ixi);
    out_cmyo{i} = youti(nstates+1,:);
    out_csr{i} = youti(nstates+2,:);
    out_theta{i} = youti(nstates+3,:);
    out_thalf{i} = youti(nstates+4,:);
    out_phi_lcc_BAR{i} = youti(nstates+5,:);
    out_phi_plb_BAR{i} = youti(nstates+6,:);
    
    out_Ia{i} = zeros(length(tstep{i}),1); 
    out_Ilcc{i} = zeros(length(tstep{i}),1); 
    out_Incx{i} = zeros(length(tstep{i}),1); 
    out_Ipmca{i} = zeros(length(tstep{i}),1); 
    out_Ileak{i} = zeros(length(tstep{i}),1); 
    out_Jds{i} = zeros(length(tstep{i}),1); 
    out_Jserca{i} = zeros(length(tstep{i}),1); 
    out_Jncx{i} = zeros(length(tstep{i}),1); 
    out_Jpmca{i} = zeros(length(tstep{i}),1); 
    out_Jpmleak{i} = zeros(length(tstep{i}),1); 
    out_Jsrleak{i} = zeros(length(tstep{i}),1);
    out_Jryr{i} = zeros(length(tstep{i}),1);
    out_Em{i} = zeros(length(tstep{i}),1);
    
    for j=1:length(tstep{i})
        [~,out2] = DYDT(tstep{i}(j),youti(:,j));
        out_Ia{i}(j) = out2(1); 
        out_Ilcc{i}(j) = out2(2);
        out_Incx{i}(j) = out2(3);
        out_Ipmca{i}(j) = out2(4);
        out_Ileak{i}(j) = out2(5);
        out_Jds{i}(j) = out2(6);
        out_Jserca{i}(j) = out2(7);
        out_Jncx{i}(j) = out2(8);
        out_Jpmca{i}(j) = out2(9);
        out_Jpmleak{i}(j) = out2(10);
        out_Jsrleak{i}(j) = out2(11);
        out_Jryr{i}(j) = out2(12);
        out_Em{i}(j) = out2(13);
    end
    tstep{i} = tstep{i} - tstep{i}(1);
    out_Em{i}(1) = out_Em{i}(end);
end

%% Compare to Dibb Figure 1
Cai_amplitude_model = [(max(out_cmyo{1})-out_cmyo{1}(1)); ...
                       (max(out_cmyo{2})-out_cmyo{2}(1)); ...
                       (max(out_cmyo{4})-out_cmyo{4}(1))]*1e6;
                   
Cai_diastolic_model = [out_cmyo{1}(1); out_cmyo{2}(1); out_cmyo{4}(1)]*1e6;
           
Cai_4Hz = []; Cai_6Hz = []; Cai_8Hz = [];
Cai_amplitude = []; Cai_amplitude_err = [];
Cai_diastolic = []; Cai_diastolic_err = [];
SRC_control = []; SRC_control_err = []; t_switch = [];
load Dibb_Fig1_data

AP_4Hz_ctl = []; AP_6Hz_ctl = []; AP_8Hz_ctl = [];
AP_4Hz_iso = []; AP_6Hz_iso = []; AP_8Hz_iso = [];
load ActionPotentials

Cai_4Hz_norm = (Cai_4Hz(:,2)-Cai_4Hz(1,2))/(max(Cai_4Hz(:,2)) - Cai_4Hz(1,2));
Cai_6Hz_norm = (Cai_6Hz(:,2)-Cai_6Hz(1,2))/(max(Cai_6Hz(:,2)) - Cai_6Hz(1,2));
Cai_8Hz_norm = (Cai_8Hz(:,2)-Cai_8Hz(1,2))/(max(Cai_8Hz(:,2)) - Cai_8Hz(1,2));
Cai_4Hz_norm_model = (out_cmyo{1} - out_cmyo{1}(1))/(max(out_cmyo{1})-out_cmyo{1}(1));
Cai_6Hz_norm_model = (out_cmyo{2} - out_cmyo{2}(1))/(max(out_cmyo{2})-out_cmyo{2}(1));
Cai_8Hz_norm_model = (out_cmyo{4} - out_cmyo{4}(1))/(max(out_cmyo{4})-out_cmyo{4}(1));

cvrmse1 = (sqrt(sum((interp1(tstep{1},Cai_4Hz_norm_model,Cai_4Hz(:,1)*1e-3) - Cai_4Hz_norm).^2)/length(Cai_4Hz_norm))/mean(Cai_4Hz_norm) + ...
           sqrt(sum((interp1(tstep{2},Cai_6Hz_norm_model,Cai_6Hz(:,1)*1e-3) - Cai_6Hz_norm).^2)/length(Cai_6Hz_norm))/mean(Cai_6Hz_norm) + ...
           sqrt(sum((interp1(tstep{4},Cai_8Hz_norm_model,Cai_8Hz(:,1)*1e-3) - Cai_8Hz_norm).^2)/length(Cai_8Hz_norm))/mean(Cai_8Hz_norm))/3;
  
cvrmse2 = sqrt(sum((Cai_amplitude_model - Cai_amplitude(:,2)).^2)/length(Cai_amplitude(:,2)))/mean(Cai_amplitude(:,2));
cvrmse3 = sqrt(sum((Cai_diastolic_model - Cai_diastolic(:,2)).^2)/length(Cai_diastolic(:,2)))/mean(Cai_diastolic(:,2));

chisquare1 = sum(((Cai_amplitude_model - Cai_amplitude(:,2))./(Cai_amplitude_err(:,2)*sqrt(12))).^2);
chisquare2 = sum(((Cai_diastolic_model - Cai_diastolic(:,2))./(Cai_diastolic_err(:,2)*sqrt(12))).^2);

figure(1);
plot(tout,yout(nstates+1,:)*1e6,'k');
title('Figure 5A');
xlabel('time (s)');
ylabel('[Ca]_i (nM)');
ylim([0,2600]);

figure(2);
hold on;
plot(tstep{1}*1e3,Cai_4Hz_norm_model,'k',Cai_4Hz(:,1),Cai_4Hz_norm,'ok');
plot(tstep{2}*1e3,Cai_6Hz_norm_model,'b',Cai_6Hz(:,1),Cai_6Hz_norm,'ob');
plot(tstep{4}*1e3,Cai_8Hz_norm_model,'g',Cai_8Hz(:,1),Cai_8Hz_norm,'og');
hold off;
title('Figure 5B');
xlabel('time (ms)');
ylabel('Normalized [Ca]_i');

figure(3);
hold on;
plot(tstep{1}*1e3,out_Em{1}*1e3,'k',AP_4Hz_ctl(:,1),AP_4Hz_ctl(:,2),'ok');
plot(tstep{2}*1e3,out_Em{2}*1e3,'b',AP_6Hz_ctl(:,1),AP_6Hz_ctl(:,2),'ob');
plot(tstep{4}*1e3,out_Em{4}*1e3,'g',AP_8Hz_ctl(:,1),AP_8Hz_ctl(:,2),'og');
hold off;
title('Figure S2A');
xlabel('time (ms)');
ylabel('V_m (mV)');

figure(4);
hold on;
errorbar(1,Cai_amplitude(1,2),Cai_amplitude_err(1,2),'ok','MarkerFaceColor','auto','LineStyle','none');
errorbar(2,Cai_amplitude(2,2),Cai_amplitude_err(2,2),'ob','MarkerFaceColor','auto','LineStyle','none');
errorbar(3,Cai_amplitude(3,2),Cai_amplitude_err(3,2),'og','MarkerFaceColor','auto','LineStyle','none');
plot([1;2;3],Cai_amplitude_model,'-k');
errorbar(1,Cai_diastolic(1,2),Cai_diastolic_err(1,2),'sk','MarkerFaceColor','none','LineStyle','none');
errorbar(2,Cai_diastolic(2,2),Cai_diastolic_err(2,2),'sb','MarkerFaceColor','none','LineStyle','none');
errorbar(3,Cai_diastolic(3,2),Cai_diastolic_err(3,2),'sg','MarkerFaceColor','none','LineStyle','none');
plot([1;2;3],Cai_diastolic_model,'--k');
set(gca,'XTick',[1,2,3],'XTickLabel',{'4Hz','6Hz','8Hz'});
hold off;
title('Figure 5D');
ylabel('[Ca]_i (nM)');
ylim([0,400]);
xlim([0.8,3.2]);


%% Compare to Dibb Figure 8
Cai_amplitude_model = [(max(out_cmyo{1})-out_cmyo{1}(1)); ...
                       (max(out_cmyo{6})-out_cmyo{6}(1)); ...
                       (max(out_cmyo{7})-out_cmyo{7}(1)); ...
                       (max(out_cmyo{9})-out_cmyo{9}(1))]*1e6;
                   
Cai_diastolic_model = [out_cmyo{1}(1); out_cmyo{6}(1); out_cmyo{7}(1); out_cmyo{9}(1)]*1e6;
           
Cai_4Hz_con = []; Cai_4Hz_iso = []; Cai_6Hz_iso = []; Cai_8Hz_iso = [];
Cai_amplitude = []; Cai_amplitude_err = [];
Cai_diastolic = []; Cai_diastolic_err = [];
SRC_iso = []; SRC_iso_err = []; t_switch = [];
load Dibb_Fig8_data

Cai_4Hz_con_norm = (Cai_4Hz_con(:,2)-Cai_4Hz_con(1,2))/(max(Cai_4Hz_con(:,2)) - Cai_4Hz_con(1,2));
Cai_4Hz_iso_norm = (Cai_4Hz_iso(:,2)-Cai_4Hz_iso(1,2))/(max(Cai_4Hz_iso(:,2)) - Cai_4Hz_iso(1,2));
Cai_6Hz_iso_norm = (Cai_6Hz_iso(:,2)-Cai_6Hz_iso(1,2))/(max(Cai_6Hz_iso(:,2)) - Cai_6Hz_iso(1,2));
Cai_8Hz_iso_norm = (Cai_8Hz_iso(:,2)-Cai_8Hz_iso(1,2))/(max(Cai_8Hz_iso(:,2)) - Cai_8Hz_iso(1,2));
Cai_4Hz_con_norm_model = (out_cmyo{1} - out_cmyo{1}(1))/(max(out_cmyo{1})-out_cmyo{1}(1));
Cai_4Hz_iso_norm_model = (out_cmyo{6} - out_cmyo{6}(1))/(max(out_cmyo{6})-out_cmyo{6}(1));
Cai_6Hz_iso_norm_model = (out_cmyo{7} - out_cmyo{7}(1))/(max(out_cmyo{7})-out_cmyo{7}(1));
Cai_8Hz_iso_norm_model = (out_cmyo{9} - out_cmyo{9}(1))/(max(out_cmyo{9})-out_cmyo{9}(1));

cvrmse4 = (sqrt(sum((interp1(tstep{1},Cai_4Hz_con_norm_model,Cai_4Hz_con(:,1)*1e-3) - Cai_4Hz_con_norm).^2)/length(Cai_4Hz_con_norm))/mean(Cai_4Hz_con_norm) + ...
           sqrt(sum((interp1(tstep{6},Cai_4Hz_iso_norm_model,Cai_4Hz_iso(:,1)*1e-3) - Cai_4Hz_iso_norm).^2)/length(Cai_4Hz_iso_norm))/mean(Cai_4Hz_iso_norm) + ...
           sqrt(sum((interp1(tstep{7},Cai_6Hz_iso_norm_model,Cai_6Hz_iso(:,1)*1e-3) - Cai_6Hz_iso_norm).^2)/length(Cai_6Hz_iso_norm))/mean(Cai_6Hz_iso_norm) + ...
           sqrt(sum((interp1(tstep{9},Cai_8Hz_iso_norm_model,Cai_8Hz_iso(:,1)*1e-3) - Cai_8Hz_iso_norm).^2)/length(Cai_8Hz_iso_norm))/mean(Cai_8Hz_iso_norm))/4;
  
cvrmse5 = sqrt(sum((Cai_amplitude_model - Cai_amplitude(:,2)).^2)/length(Cai_amplitude(:,2)))/mean(Cai_amplitude(:,2));
cvrmse6 = sqrt(sum((Cai_diastolic_model - Cai_diastolic(:,2)).^2)/length(Cai_diastolic(:,2)))/mean(Cai_diastolic(:,2));

chisquare3 = sum(((Cai_amplitude_model - Cai_amplitude(:,2))./(Cai_amplitude_err(:,2)*sqrt(8))).^2);
chisquare4 = sum(((Cai_diastolic_model - Cai_diastolic(:,2))./(Cai_diastolic_err(:,2)*sqrt(8))).^2);

figure(5);
hold on;
plot(tstep{1}*1e3,Cai_4Hz_con_norm_model,'c',Cai_4Hz_con(:,1),Cai_4Hz_con_norm,'oc');
plot(tstep{6}*1e3,Cai_4Hz_iso_norm_model,'k',Cai_4Hz_iso(:,1),Cai_4Hz_iso_norm,'ok');
plot(tstep{7}*1e3,Cai_6Hz_iso_norm_model,'b',Cai_6Hz_iso(:,1),Cai_6Hz_iso_norm,'ob');
plot(tstep{9}*1e3,Cai_8Hz_iso_norm_model,'g',Cai_8Hz_iso(:,1),Cai_8Hz_iso_norm,'og');
hold off;
title('Figure 5C');
xlabel('time (ms)');
ylabel('Normalized [Ca]_i');

figure(6);
hold on;
plot(tstep{6}*1e3,out_Em{6}*1e3,'k',AP_4Hz_iso(:,1),AP_4Hz_iso(:,2),'ok');
plot(tstep{7}*1e3,out_Em{7}*1e3,'b',AP_6Hz_iso(:,1),AP_6Hz_iso(:,2),'ob');
plot(tstep{9}*1e3,out_Em{9}*1e3,'g',AP_8Hz_iso(:,1),AP_8Hz_iso(:,2),'og');
hold off;
title('Figure S2B');
xlabel('time (ms)');
ylabel('V_m (mV)');

figure(7);
hold on;
errorbar(1,Cai_amplitude(1,2),Cai_amplitude_err(1,2),'oc','MarkerFaceColor','auto','LineStyle','none');
errorbar(2,Cai_amplitude(2,2),Cai_amplitude_err(2,2),'ok','MarkerFaceColor','auto','LineStyle','none');
errorbar(3,Cai_amplitude(3,2),Cai_amplitude_err(3,2),'ob','MarkerFaceColor','auto','LineStyle','none');
errorbar(4,Cai_amplitude(4,2),Cai_amplitude_err(4,2),'og','MarkerFaceColor','auto','LineStyle','none');
plot([1;2;3;4],Cai_amplitude_model,'-k');
errorbar(1,Cai_diastolic(1,2),Cai_diastolic_err(1,2),'sc','MarkerFaceColor','none','LineStyle','none');
errorbar(2,Cai_diastolic(2,2),Cai_diastolic_err(2,2),'sk','MarkerFaceColor','none','LineStyle','none');
errorbar(3,Cai_diastolic(3,2),Cai_diastolic_err(3,2),'sb','MarkerFaceColor','none','LineStyle','none');
errorbar(4,Cai_diastolic(4,2),Cai_diastolic_err(4,2),'sg','MarkerFaceColor','none','LineStyle','none');
plot([1;2;3;4],Cai_diastolic_model,'--k');
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'4Hz ctl','4Hz+iso','6Hz+iso','8Hz+iso'});
hold off;
title('Figure 5E');
ylabel('[Ca]_i (nM)');
ylim([0,1800]);
xlim([0.8,4.2]);

%% Simulate caffeine experiments
yout1 = yout;
tout1 = tout;
tf1 = tf;
ix1 = [1,2,4,6,7,9];
t0 = [0; 0; 0; 0; 0; 0];
tf = [2; 2; 2; 2; 2; 2];
HRstep = [4; 6; 8; 4; 6; 8];
ciso_step = [0, 0, 0, 30e-6, 30e-6, 30e-6];
phi_lcc_CaMK_step = [phi_lcc_4Hz_CaMK_con, phi_lcc_6Hz_CaMK_con, phi_lcc_8Hz_CaMK_con, ...
                     phi_lcc_4Hz_CaMK_iso, phi_lcc_6Hz_CaMK_iso, phi_lcc_8Hz_CaMK_iso];
phi_plb_CaMK_step = [phi_plb_4Hz_CaMK_con, phi_plb_6Hz_CaMK_con, phi_plb_8Hz_CaMK_con, ...
                     phi_plb_4Hz_CaMK_iso, phi_plb_6Hz_CaMK_iso, phi_plb_8Hz_CaMK_iso];
options = odeset('InitialStep',1e-6);

youti = cell(6,1);
tstep = cell(6,1);
youti_long = cell(6,1);
tstep_long = cell(6,1);

for i=1:length(t0)
    HR = HRstep(i);
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    ix = find(tout1 <= tf1(ix1(i)), 1, 'last');
    y0 = yout1(:,ix); %initial condition for odes
    tout = 0;
    yout = y0;
    while tout(end) < tf(i) - 1e-2 %note: the term 1e-2 is an arbitrary 
                                   %small number to account for rounding errors
        if tout(end) - 1/HR >= t0(i)
            ix = find(tout > tout(end) - 1/HR,1,'first');
            options = odeset('InitialStep',tout(ix) - tout(ix-1));
        end
        sol = ode15s(@DYDT,[tout(end),tout(end) + 1/HR],y0,options);
        if sol.x(end) ~= tout(end) + 1/HR % solver failed
            return
        end
        tout = [tout, sol.x(2:end)];
        yout = [yout, sol.y(:,2:end)];
        y0 = sol.y(:,end);
    end
    HR = 0;
    vsrleak0 = vsrleak;
    vsrleak = vsrleak0*1e4;
    sol = ode15s(@DYDT,[tout(end),tout(end) + 6],y0,options);
    tout = [tout, sol.x(2:end)];
    yout = [yout, sol.y(:,2:end)];
    youti{i} = sol.y;
    tstep{i} = sol.x;
    youti_long{i} = yout;
    tstep_long{i} = tout;
    vsrleak = vsrleak0;
end

out_Jncx = cell(6,1);
SRC_model = zeros(6,1);
for i=1:length(t0);
    HR = 0;
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    out_Jncx{i} = zeros(length(tstep{i}),1);
    for j=1:length(tstep{i})
        [~,out2] = DYDT(tstep{i}(j),youti{i}(:,j));
        out_Jncx{i}(j) = out2(8);
    end
    fun = @(tprime) interp1(tstep{i},out_Jncx{i},tprime);
    SRC_model(i) = 1.4706*1e3*quad(fun,tstep{i}(1),tstep{i}(end));
    tstep{i} = tstep{i} - tstep{i}(1);
end

cvrmse7 = sqrt(sum((SRC_model(1:3) - SRC_control(:,2)).^2)/length(SRC_control(:,2)))/mean(SRC_control(:,2));
cvrmse8 = sqrt(sum((SRC_model([1,4:6]) - SRC_iso(:,2)).^2)/length(SRC_iso(:,2)))/mean(SRC_iso(:,2));

chisquare5 = sum(((SRC_model(1:3) - SRC_control(:,2))./(SRC_control_err(:,2)*sqrt(12))).^2);
chisquare6 = sum(((SRC_model([1,4:6]) - SRC_iso(:,2))./(SRC_iso_err(:,2)*sqrt(8))).^2);

Em_model = cell(6,1);
for i = 1:length(t0)
    Em_model{i} = zeros(1,length(tstep_long{i}));
    HR = HRstep(i);
    ciso = ciso_step(i);
    phi_lcc_CaMK = phi_lcc_CaMK_step(i);
    phi_plb_CaMK = phi_plb_CaMK_step(i);
    for j=1:length(tstep_long{i})
        if tstep_long{i}(j) >= tf(i)
            HR = 0;
        end
        [~,out2] = DYDT(tstep_long{i}(j),youti_long{i}(:,j));
        Em_model{i}(j) = out2(13);
    end
end

figure(8);
subplot(2,1,2); plot(tstep_long{1},youti_long{1}(nstates+1,:),'k');
xlabel('time (s)');
ylabel('[Ca]_i (mM)');
subplot(2,1,1); plot(tstep_long{1},Em_model{1},'k');
ylabel('V_m (V)');
title('Figure 6A');

figure(9);
hold on;
bar(1,SRC_control(1,2),'w','EdgeColor','k');
bar(2,SRC_control(2,2),'w','EdgeColor','b');
bar(3,SRC_control(3,2),'w','EdgeColor','g');
errorbar(1,SRC_control(1,2),0,SRC_control_err(1,2),'k','Marker','none','LineStyle','none');
errorbar(2,SRC_control(2,2),0,SRC_control_err(2,2),'b','Marker','none','LineStyle','none');
errorbar(3,SRC_control(3,2),0,SRC_control_err(3,2),'g','Marker','none','LineStyle','none');
plot([1;2;3],SRC_model(1:3),'-k');
set(gca,'XTick',[1,2,3],'XTickLabel',{'4Hz','6Hz','8Hz'});
hold off;
title('Figure 6B');
ylabel('SR Ca content (micromol/L)');
ylim([0,140]);

figure(10);
hold on;
bar(1,SRC_iso(1,2),'w','EdgeColor','c');
bar(2,SRC_iso(2,2),'w','EdgeColor','k');
bar(3,SRC_iso(3,2),'w','EdgeColor','b');
bar(4,SRC_iso(4,2),'w','EdgeColor','g');
errorbar(1,SRC_iso(1,2),0,SRC_iso_err(1,2),'c','Marker','none','LineStyle','none');
errorbar(2,SRC_iso(2,2),0,SRC_iso_err(2,2),'k','Marker','none','LineStyle','none');
errorbar(3,SRC_iso(3,2),0,SRC_iso_err(3,2),'b','Marker','none','LineStyle','none');
errorbar(4,SRC_iso(4,2),0,SRC_iso_err(4,2),'g','Marker','none','LineStyle','none');
plot([1;2;3;4],SRC_model([1,4:6]),'-k');
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'4Hz ctl','4Hz+iso','6Hz+iso','8Hz+iso'});
hold off;
title('Figure 6C');
ylabel('SR Ca content (micromol/L)');
ylim([0,200]);

scvrmse = cvrmse1 + cvrmse2 + cvrmse3 + cvrmse4 + cvrmse5 + cvrmse6 + cvrmse7 + cvrmse8
AIC = 2*28 + chisquare1 + chisquare2 + chisquare3 + chisquare4 + chisquare5 + chisquare6

%% DYDT function
function [ydot,out2] = DYDT(t,y)
z1_mode1 = y(1);
z2_mode1 = y(2);
z3_mode1 = y(3);
z4_mode1 = 1 - y(1) - y(2) - y(3);
z1_mode2 = y(4);
z2_mode2 = y(5);
z3_mode2 = y(6);
z4_mode2 = 1 - y(4) - y(5) - y(6);

cmyo = y(7);
csr = y(8);
theta = y(9);
thalf = y(10);
phi_lcc_BAR = y(11);
phi_plb_BAR = y(12);

ttilde = mod(t,1/HR);

% electrophysiology
thalfo = (thalfo_max - thalfo_min)*theta + thalfo_min;
a = a_max*theta;
eta = (eta_max - eta_min)*theta + eta_min;
Eo = (Eo_max - Eo_min)*theta + Eo_min;
theta_ss = 1/(1 + (Kiso_AP/ciso)^Hiso_AP);
thalf_ss = thalfo - a*HR;
if isnan(ttilde)
    Em = -80e-3;
else
    Em = Eo + (Emax - Eo)*(0.5 - 1/pi*atan(log((ttilde/thalf)^eta)));
end

% calcium handling
cds00 = cmyo;
cds01 = (vds*cmyo + vryr*csr)/(vds + vryr);
if abs(Em) < 1e-10
    cds10 = (vds*cmyo + vlcc*Cao)/(vds + vlcc);
    cds11 = (vds*cmyo + vlcc*Cao + vryr*csr)/(vds + vlcc + vryr); 
else                  
    cds10 = (vds*cmyo + vlcc*Em*2*F/R/T*Cao/(exp(Em*2*F/R/T) - 1))/ ...
            (vds + vlcc*Em*2*F/R/T/(1 - exp(-Em*2*F/R/T)));
    cds11 = (vds*cmyo + vlcc*Em*2*F/R/T*Cao/(exp(Em*2*F/R/T) - 1) + vryr*csr)/ ...
            (vds + vlcc*Em*2*F/R/T/(1 - exp(-Em*2*F/R/T)) + vryr);
end
    
% diad space rate laws
alphaf = kopen_lcc/(1+exp(-(Em-Vhalf_lcc)/deltaV_lcc));
alphar_mode1 = kclose_lcc;
alphar_mode2 = zetalcc*kclose_lcc;
epsf00 = cds00*(kinactmin_lcc + (kinactmax_lcc - kinactmin_lcc)/(1+exp(-(Em - Vhalf_lcc)/deltaV_lcc)));
epsf01 = cds01*(kinactmin_lcc + (kinactmax_lcc - kinactmin_lcc)/(1+exp(-(Em - Vhalf_lcc)/deltaV_lcc)));
epsr = kdeinactmin_lcc + (kdeinactmax_lcc - kdeinactmin_lcc)/(1+exp((Em - Vhalf_lcc)/deltaV_lcc));
betaf00 = kopen_ryr/(1 + (Kryr/cds00)^2);
betaf10 = kopen_ryr/(1 + (Kryr/cds10)^2);
betar = kclose_ryr;
muf00 = kinactmin_ryr + (kinactmax_ryr - kinactmin_ryr)/(1 + (Kryr/cds00)^2);
muf10 = kinactmin_ryr + (kinactmax_ryr - kinactmin_ryr)/(1 + (Kryr/cds10)^2);
mur00 = kdeinactmax_ryr - (kdeinactmax_ryr - kdeinactmin_ryr)/(1 + (Kryr/cds00)^2);
mur10 = kdeinactmax_ryr - (kdeinactmax_ryr - kdeinactmin_ryr)/(1 + (Kryr/cds10)^2);

% note: the below ugly expressions were obtained from MAPLE
Pyocz1_mode1 = (alphaf+alphar_mode1+betaf00+betar)*alphaf*betar/...
               ((alphaf+alphar_mode1)*(alphaf*betaf10+alphaf*betar+betaf00*betar+betar^2+betar*alphar_mode1+betaf00*betaf10+alphar_mode1*betaf00+betar*betaf10));
Pycoz1_mode1 = alphar_mode1*(alphaf*betaf10+betaf00*betaf10+alphar_mode1*betaf00+betaf00*betar)/...
               (alphaf*betaf00*betar+alphar_mode1*alphaf*betaf00+alphaf^2*betar+alphaf*betar^2+2*alphaf*betar*alphar_mode1+alphar_mode1*betaf00*betar+...
               betaf00*alphar_mode1^2+alphar_mode1*betar^2+betar*alphar_mode1^2+betaf00*betaf10*alphaf+betaf00*betaf10*alphar_mode1+alphaf^2*betaf10+...
               alphaf*betaf10*alphar_mode1+betar*betaf10*alphaf+betar*betaf10*alphar_mode1);
Pyccz1_mode1 = alphar_mode1*betar*(alphaf+betar+betaf10+alphar_mode1)/...
               (alphaf*betaf00*betar+alphar_mode1*alphaf*betaf00+alphaf^2*betar+alphaf*betar^2+2*alphaf*betar*alphar_mode1+alphar_mode1*betaf00*betar+...
               betaf00*alphar_mode1^2+alphar_mode1*betar^2+betar*alphar_mode1^2+betaf00*betaf10*alphaf+betaf00*betaf10*alphar_mode1+alphaf^2*betaf10+...
               alphaf*betaf10*alphar_mode1+betar*betaf10*alphaf+betar*betaf10*alphar_mode1);
Pyooz1_mode1 = 1 - Pyocz1_mode1 - Pycoz1_mode1 - Pyccz1_mode1;
Pyoiz2_mode1 = alphaf/(alphaf+alphar_mode1);
Pyciz2_mode1 = alphar_mode1/(alphaf+alphar_mode1);
Pyocz1_mode2 = (alphaf+alphar_mode2+betaf00+betar)*alphaf*betar/...
               ((alphaf+alphar_mode2)*(alphaf*betaf10+alphaf*betar+betaf00*betar+betar^2+betar*alphar_mode2+betaf00*betaf10+alphar_mode2*betaf00+betar*betaf10));
Pycoz1_mode2 = alphar_mode2*(alphaf*betaf10+betaf00*betaf10+alphar_mode2*betaf00+betaf00*betar)/...
               (alphaf*betaf00*betar+alphar_mode2*alphaf*betaf00+alphaf^2*betar+alphaf*betar^2+2*alphaf*betar*alphar_mode2+alphar_mode2*betaf00*betar+...
               betaf00*alphar_mode2^2+alphar_mode2*betar^2+betar*alphar_mode2^2+betaf00*betaf10*alphaf+betaf00*betaf10*alphar_mode2+alphaf^2*betaf10+...
               alphaf*betaf10*alphar_mode2+betar*betaf10*alphaf+betar*betaf10*alphar_mode2);
Pyccz1_mode2 = alphar_mode2*betar*(alphaf+betar+betaf10+alphar_mode2)/...
               (alphaf*betaf00*betar+alphar_mode2*alphaf*betaf00+alphaf^2*betar+alphaf*betar^2+2*alphaf*betar*alphar_mode2+alphar_mode2*betaf00*betar+...
               betaf00*alphar_mode2^2+alphar_mode2*betar^2+betar*alphar_mode2^2+betaf00*betaf10*alphaf+betaf00*betaf10*alphar_mode2+alphaf^2*betaf10+...
               alphaf*betaf10*alphar_mode2+betar*betaf10*alphaf+betar*betaf10*alphar_mode2);
Pyooz1_mode2 = 1 - Pyocz1_mode2 - Pycoz1_mode2 - Pyccz1_mode2;
Pyoiz2_mode2 = alphaf/(alphaf+alphar_mode2);
Pyciz2_mode2 = alphar_mode2/(alphaf+alphar_mode2);
Pyioz3 = betaf00/(betaf00 + betar);
Pyicz3 = betar/(betaf00 + betar);

r1_mode1 = Pyocz1_mode1*muf10 + Pyccz1_mode1*muf00;
r1_mode2 = Pyocz1_mode2*muf10 + Pyccz1_mode2*muf00;
r2_mode1 = Pyoiz2_mode1*mur10 + Pyciz2_mode1*mur00;
r2_mode2 = Pyoiz2_mode2*mur10 + Pyciz2_mode2*mur00;
r3 = Pyicz3*muf00;
r4 = mur00;
r5_mode1 = Pyccz1_mode1*epsf00 + Pycoz1_mode1*epsf01;
r5_mode2 = Pyccz1_mode2*epsf00 + Pycoz1_mode2*epsf01;
r6 = epsr;
r7_mode1 = Pyciz2_mode1*epsf00;
r7_mode2 = Pyciz2_mode2*epsf00;
r8 = epsr;

phi_lcc = 1 - (1 - phi_lcc_BAR)*(1 - phi_lcc_CaMK);
phi_plb = 1 - (1 - phi_plb_BAR)*(1 - phi_plb_CaMK);

P01 = (1 - phi_lcc)*(z1_mode1*Pycoz1_mode1 + z3_mode1*Pyioz3) + ...
      phi_lcc*(z1_mode2*Pycoz1_mode2 + z3_mode2*Pyioz3);
P10 = (1 - phi_lcc)*(z1_mode1*Pyocz1_mode1 + z2_mode1*Pyoiz2_mode1) + ...
      phi_lcc*(z1_mode2*Pyocz1_mode2 + z2_mode2*Pyoiz2_mode2);
P11 = (1 - phi_lcc)*z1_mode1*Pyooz1_mode1 + ...
      phi_lcc*z1_mode2*Pyooz1_mode2;

% diad space flux
Jds = P01*vds*(cds01 - cmyo) + P10*vds*(cds10 - cmyo) + P11*vds*(cds11 - cmyo);
Jryr = P01*vryr*(csr - cds01) + P11*vryr*(csr - cds11);
Jlcc = P10*vlcc*Em*2*F/R/T*(Cao*exp(-Em*2*F/R/T) - cds10)/(1 - exp(-Em*2*F/R/T)) + ...
       P11*vlcc*Em*2*F/R/T*(Cao*exp(-Em*2*F/R/T) - cds11)/(1 - exp(-Em*2*F/R/T));

% serca flux
Jserca = (1 - phi_plb)*Jsercao/(1 + (Kserca/cmyo)^2) + ...
          phi_plb*Jsercao/(1 + (zetaserca*Kserca/cmyo)^2);
     
% other fluxes
Jncx = -Jncxo*(Nai^3*exp(eta_ncx*F*Em/R/T)*Cao - Nao^3*exp((eta_ncx-1)*F*Em/R/T)*cmyo)/ ...
               ((Kncxn^3 + Nao^3)*(Kncxc + Cao)*(1 + kncxsat*exp((eta_ncx-1)*F*Em/R/T)));
Jpmca = Jpmca_max/(1 + (Kpmca/cmyo));
ECa = R*T/2/F*log(Cao/cmyo);
Jpmleak = gpmleak*(ECa - Em);
Jsrleak = vsrleak*(csr - cmyo);

% currents
Ilcc = 2*Jlcc*Vtot*F*1e9; %nA
Incx = Jncx*Vtot*F*1e9; %nA
Ipmca = 0; %pmca is electroneutral
Ileak = 2*Jpmleak*Vtot*F*1e9; %nA

% buffering
beta_myo = 1 + Bfluo*Kfluo/(cmyo + Kfluo)^2 + Bbuff_myo*Kdmyo/(cmyo + Kdmyo)^2;
beta_sr = 1 + Bbuff_sr*Kdsr/(csr + Kdsr)^2;

% beta-AR regulation
phi_lcc_BAR_ss = phi_lcc_BAR_max/(1 + (Kiso_lcc/ciso)^Hiso_lcc);
phi_plb_BAR_ss = phi_plb_BAR_max/(1 + (Kiso_plb/ciso)^Hiso_plb);

% derivatives
ydot = zeros(12,1);
ydot(1) = z2_mode1*r2_mode1 + z3_mode1*r6 - z1_mode1*(r5_mode1 + r1_mode1);
ydot(2) = z1_mode1*r1_mode1 + z4_mode1*r8 - z2_mode1*(r2_mode1 + r7_mode1);
ydot(3) = z1_mode1*r5_mode1 + z4_mode1*r4 - z3_mode1*(r6 + r3);
ydot(4) = z2_mode2*r2_mode2 + z3_mode2*r6 - z1_mode2*(r5_mode2 + r1_mode2);
ydot(5) = z1_mode2*r1_mode2 + z4_mode2*r8 - z2_mode2*(r2_mode2 + r7_mode2);
ydot(6) = z1_mode2*r5_mode2 + z4_mode2*r4 - z3_mode2*(r6 + r3);

ydot(7) = (Jds + Jpmleak + Jsrleak - Jncx - Jpmca - Jserca)/lambda_myo/beta_myo;
ydot(8) = (Jserca - Jryr - Jsrleak)/lambda_sr/beta_sr;
ydot(9) = (-theta + theta_ss)/tauAPD_BAR;
ydot(10) = (-thalf + thalf_ss)/tauAPD_H;
ydot(11) = (-phi_lcc_BAR + phi_lcc_BAR_ss)/tauBAR_lcc;
ydot(end) = (-phi_plb_BAR + phi_plb_BAR_ss)/tauBAR_plb;

if ttilde > 5e-3
    Ia = 0;
else
    Ia = 1;
end

% other output
out2 = [Ia; Ilcc; Incx; Ipmca; Ileak; Jds; Jserca; Jncx; Jpmca; Jpmleak; Jsrleak; Jryr; Em];
end
end