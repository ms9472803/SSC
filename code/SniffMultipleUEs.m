clear variables;
close all;
load('SatCH.mat');
load('NoiseParam.mat');

TransmitterPower_dB = 60; % dBm
%LNAGain_dB = 35; % dBm
PowerVar = 10^(TransmitterPower_dB/10);
%PowerVar = 10^((TransmitterPower_dB+LNAGain_dB)/10);
NumTransmitAntennas = 3;
NumReceiveAntennas = 3;
%H = complex(randn(NumTransmitAntennas,NumReceiveAntennas), randn(NumTransmitAntennas,NumReceiveAntennas)); % random genereate downlink CSI
H = SatCH(randi(length(SatCH),NumReceiveAntennas,NumTransmitAntennas)); % random from Sat Channel

%imperfectH = H + err;


%% SNR calculation

Eb_N0_dB = cell2mat(Eb_N0_dB);
Eb_N0_dB_MAX = Eb_N0_dB;
Eb_N0_dB = Eb_N0_dB_MAX-30:2:Eb_N0_dB_MAX-10; % Es/N0 in dB

RcvrPower_dB = cell2mat(RcvrPower_dB);
Eb_N0 = 10.^(Eb_N0_dB./10);
RcvrPower = 10.^(RcvrPower_dB./10); 
NoiseVar = RcvrPower./Eb_N0; 
snr = 20;


%% mask CSI 
theta_max = [45 22.5 15 11.25];
theta = [50 100 200];
alpha = [complex(cosd(theta(1)),sind(theta(1))) complex(cosd(theta(2)),sind(theta(2))) complex(cosd(theta(3)),sind(theta(3)))];
mask = alpha;

mask_CSI = zeros(1,3); % initialization ; BS receive from UE's feedback

for i = 1:3
    mask_CSI(i,:) = H(i,:)*mask(i);
end
precoding_mask_ByCh = H*pinv(mask_CSI); % H*H^+ = diagonal matrix: [bar(alpha1); bar(alpha2); bar(alpha3)]

%% mask pilot symbol (refer to "Securing Channel State Information in Multiuser MIMO With Limited Feedback")
pilot = complex(1/sqrt(2),1/sqrt(2)); % I=1/sqrt(2), Q=1/sqrt(2) DVB-S2
mask_pilot = mask*pilot;
Y_mask_pilot = zeros(3,3);
for i = 1:3
    Y_mask_pilot(i,:) = H(i,:)*mask_pilot(i); % Y=HX
end
CSI_mask_pilot = Y_mask_pilot/pilot; % UE estimate CSI_mask_pilot and feedback to BS

BS_recover_CSI = zeros(3,3);
for i = 1:3
    BS_recover_CSI(i,:) = CSI_mask_pilot(i,:)/mask(i); % BS recover CSI for AMC and precoding (using: CSI_mask_pilot=mask*H)
end
%ChByPrecoding = H*pinv(BS_recover_CSI); % H*H^+ = identity matrix

%% BS transmit data using zero forcing precoding
symbolNum = 32000; %一共有幾個符號
M_QPSK = 4; % QPSK order
M_8PSK = 8;
M_16APSK = 16; 
M_32APSK = 32;
MOD = M_QPSK;

x = [];
txSig = [];

spacing = symbolNum/MOD;
temp = [1:spacing:symbolNum];
for i = 1:NumReceiveAntennas
   %x(i,:) = randi([0 MOD-1],1,symbolNum); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
   
   for k = 0:length(temp)-1
       x(i,1+spacing*k:spacing*(k+1)) = k;
   end

   
   if MOD== 16 || MOD == 32
       txSig(i,:) = dvbsapskmod(x(i,:),MOD,'s2'); 
   else
       txSig(i,:) = pskmod(x(i,:),MOD,pi/MOD); 
   end

end

precoding = pinv(H);
ChByPrecoding = H*precoding;

m = txSig;
%[txSig1;txSig2;txSig3];
rxSig = ChByPrecoding*m; % Y = HX + n = CPm + n

%% Sniffing attack (sniff multiple UEs)

M = 3;
Var_x = [1 1 1]; 
F = H;
F(3,:) = 0;

%竊聽者估計自己的CSI 還有竊聽UE_1, UE_2的CSI
eve_CSI(1,:) = H(1,:);
eve_CSI(2,:) =  H(2,:);
eve_CSI(3,:) = H(3,:); 

% 用公式計算F
for i=1:3
    temp = 0;
    for k=1:M-1
        temp = temp + Var_x(k)*eve_CSI(k,i);
    end
   F(3,i) = (eve_CSI(M,i)-temp);
end

pinvF = pinv(F);
pinvH = pinv(H);

totalPower = RcvrPower*norm(pinvH,'fro')^2; % P

%power_i = totalPower/ (norm(pinvF,'fro')^2 + Var_x(2)^2*norm(pinvF(:,M),'fro')^2);
power_i = totalPower/ norm(pinvF,'fro')^2;
power_i

SINR_UT_sniffing_attack = (10*log10(power_i ./ NoiseVar)).' ;

Eb_N0_Eve = 10.^( (Eb_N0_dB-5)./10);
NoiseVar_Eve = RcvrPower./Eb_N0_Eve; 


% 先預設T, 用T計算

eve_CSI_ = H;
T_ = eye(3);
Var_x_ = [1 1 1];
T_(3,:) = Var_x;

%Var_x = T(3,:); % selected by attacker

F_Eve = pinv(T_) * eve_CSI_;
F_ = H;
F_(3,:) = F_Eve(3,:);
pinvF_multi = pinv(F_);

rxSig_SniffAttack2 = H*pinvF_multi*m;
rxSig_SniffAttack2_AWGN = awgn(rxSig_SniffAttack2,snr,'measured');
sniffMessage2 = rxSig_SniffAttack2_AWGN(M,:) - m(M,:);



if MOD== 16 || MOD == 32
    sniffMessage_Label2_strong= dvbsapskdemod(sniffMessage2,MOD,'s2');
    [errSym_strong, SER_strong] = symerr(sniffMessage_Label2_strong, x(2,:));
    
    StrongSig = dvbsapskmod(sniffMessage_Label2_strong,MOD,'s2');
    sniffMessage2_remain = sniffMessage2 - StrongSig*Var_x(2);
    sniffMessage_Label2_weak = pskdemod(sniffMessage2_remain,MOD,pi/MOD);
    [errSym_weak, SER_weak] = symerr(sniffMessage_Label2_strong, x(1,:));
    
else
    sniffMessage_Label2_strong = pskdemod(sniffMessage2,MOD,pi/MOD);
    [errSym_strong, SER_strong] = symerr(sniffMessage_Label2_strong, x(2,:));
    
    StrongSig = pskmod(sniffMessage_Label2_strong,MOD,pi/MOD);
    sniffMessage2_remain = sniffMessage2 - StrongSig;
    sniffMessage_Label2_weak = pskdemod(sniffMessage2_remain,MOD,pi/MOD);
    [errSym_weak, SER_weak] = symerr(sniffMessage_Label2_strong, x(1,:));
end

power_i_multi = totalPower/ norm(pinvF_multi,'fro')^2;

SINR_UT_sniffing_attack_multi_strong = (10*log10(power_i_multi ./ NoiseVar(:))).' ;

Eb_N0_weak = 10.^( (Eb_N0_dB-3)./10);
NoiseVar_weak = RcvrPower./Eb_N0_weak; 
SINR_UT_sniffing_attack_multi_weak = (10*log10(power_i_multi ./ NoiseVar_weak(:))).' ;

SINR_multi_strong = [];
for i = 1:length(NoiseVar_Eve)
    SINR_multi_strong(i) = power_i_multi*(Var_x_(2)^2)/(NoiseVar_Eve(i)/2 + power_i_multi + power_i_multi*Var_x_(1)^2);
end
SINR_Eve_sniffing_attack_multi_strong = 10*log10(SINR_multi_strong);
%PlotSINR(Eb_N0_dB, Eb_N0_dB, SINR_UT_sniffing_attack_multi_larger, Eb_N0_dB-5, SINR_Eve_sniffing_attack_multi_larger);

SINR_multi_weak = [];
for i = 1:length(NoiseVar_Eve)
    SINR_multi_weak(i) = power_i_multi*(Var_x_(1)^2)/(NoiseVar_Eve(i)/2 + power_i_multi + power_i_multi*Var_x_(2)^2);
end
SINR_Eve_sniffing_attack_multi_weak = 10*log10(SINR_multi_weak);

% function PlotSINR_Multi(Eb_N0_dB, SINR_smaller, SINR_larger, SINR_Eve, SINR_smaller_SA, SINR_larger_SA, SINR_smaller_Eve_SA, SINR_larger_Eve_SA)
PlotSINR_Multi(Eb_N0_dB, Eb_N0_dB-3,Eb_N0_dB,Eb_N0_dB-5,SINR_UT_sniffing_attack_multi_weak,SINR_UT_sniffing_attack_multi_strong, SINR_Eve_sniffing_attack_multi_weak,SINR_Eve_sniffing_attack_multi_strong);


%% Countermeasure for sniffing attack, mask pilot symbol (SCF)

% satellite sends the same pilot symbol with same mask for each UT to estimate CSI
Y_mask_pilot = zeros(3,3); % received pilot
mask_SCF = complex(cosd(30),sind(30));
for i = 1:3
    Y_mask_pilot(i,:) = H(i,:)*(pilot*mask_SCF); % Y=HX
end
CSI_mask_pilot = Y_mask_pilot/pilot; % UE estimate CSI_mask_pilot and feedback to Sat

%竊聽者估計自己的CSI 還有竊聽UE_1, UE_2的CSI
eve_CSI_mask_pilot(1,:) = CSI_mask_pilot(1,:);
eve_CSI_mask_pilot(2,:) =  CSI_mask_pilot(2,:);
eve_CSI_mask_pilot(3,:) = CSI_mask_pilot(3,:); 

F_SCF = CSI_mask_pilot; % UE feedback estimated CSI
F_SCF(3,:) = 0;
% 用公式計算F
for i = 1:3
    temp = 0;
    for k = 1:M-1
        temp = temp + Var_x(k)*eve_CSI_mask_pilot(k,i);
    end
   F_SCF(3,i) = (eve_CSI_mask_pilot(M,i)-temp);
end

BS_recover_F_SCF = zeros(3,3);
for i = 1:3
    BS_recover_F_SCF(i,:) = F_SCF(i,:)/mask_SCF; % BS recover CSI for AMC and precoding (using: CSI_mask_pilot=mask*H)
end
T_SCF = H*pinv(BS_recover_F_SCF); % H*H^+ 

rxSig_SCF = T_SCF*m;
%sniffMessage_Countermeasure = (rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/Var_x(2);
%sniffMessage_Countermeasure = awgn(sniffMessage_Countermeasure,snr,'measured');
%sniffMessage_Label_Countermeasure = pskdemod(sniffMessage_Countermeasure,MOD,pi/MOD);
%[errSym_QPSK_Countermeasure, SER_QPSK_Countermeasure] = symerr(sniffMessage_Label_Countermeasure, x(2,:));   
% [a,b] = symerr(pskdemod(awgn((rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/((Var_x(2) * alpha3)*-1),snr,'measured'),M_QPSK,pi/M_QPSK),x2)

for snr_i = 1:length(Eb_N0_dB)
    UE2_message = awgn(rxSig_SCF(2,:),Eb_N0_dB(snr_i),'measured');
    
    sniffMessage_SCF_AWGN = awgn(rxSig_SCF(M,:),Eb_N0_dB(snr_i),'measured');
    sniffMessage_SCF = sniffMessage_SCF_AWGN - m(M,:);
    
    %sniffMessage_Countermeasure_ = awgn(sniffMessage_Countermeasure,Eb_N0_dB(snr_i),'measured');
    
    if MOD== 16 || MOD == 32
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2');
        
        sniffLabel_SCF_strong = dvbsapskdemod(sniffMessage_SCF,MOD,'s2');
        StrongSig = dvbsapskmod(sniffLabel_SCF_strong,MOD,'s2');
        remain_sniffMessage_SCF = sniffMessage_SCF - StrongSig*Var_x(2);
        sniffLabel_SCF_weak = dvbsapskdemod(remain_sniffMessage_SCF,MOD,'s2');
    else
        UE2_message_Label = pskdemod(UE2_message,MOD,pi/MOD);
        
        sniffLabel_SCF_strong = pskdemod(sniffMessage_SCF,MOD,pi/MOD);
        StrongSig = pskmod(sniffLabel_SCF_strong,MOD,pi/MOD);
        remain_sniffMessage_SCF = sniffMessage_SCF - StrongSig;
        sniffLabel_SCF_weak = pskdemod(remain_sniffMessage_SCF,MOD,pi/MOD);
    end
    
    [SCF_errSym(snr_i), SCF_SER(snr_i)] = symerr(UE2_message_Label, x(2,:));  
    
    [SCF_errSym_Eve_strong(snr_i), SCF_SER_Eve_strong(snr_i)] = symerr(sniffLabel_SCF_strong, x(2,:));  
    [SCF_errSym_Eve_weak(snr_i), SCF_SER_Eve_weak(snr_i)] = symerr(sniffLabel_SCF_weak, x(1,:));  
end

%PlotConstellation(UE2_message, txSig(2,:), 'UE_2 receive message with SCF', MOD, spacing);
%PlotConstellation(sniffMessage_SCF, txSig(2,:), 'Eavesdropping UE_2 message with SCF', MOD, spacing);


%% Countermeasure for sniffing attack, mask CSI feedback (SSC)
%竊聽者估計自己的CSI 還有竊聽UE_1, UE_2的CSI
eve_CSI_mask_CSI(1,:) = mask_CSI(1,:);
eve_CSI_mask_CSI(2,:) =  mask_CSI(2,:);
eve_CSI_mask_CSI(3,:) = mask_CSI(3,:); 

F_OUR = mask_CSI; % UE feedback estimated CSI 再做mask

% 攻擊者用公式計算F
for i = 1:3
    temp = 0;
    for k = 1:M-1
        temp = temp + Var_x(k)*eve_CSI_mask_CSI(k,i);
    end
   F_OUR(3,i) = (eve_CSI_mask_CSI(M,i)-temp);
end

BS_recover_F_OUR = F_OUR;
for i = 1:3
    BS_recover_F_OUR(i,:) = F_OUR(i,:)/mask(i);
end

T_OUR = H*pinv(BS_recover_F_OUR)
%(Var_x(2)*alpha(2))/alpha(3)
%Var_x(2) * alpha(3)'

rxSig_OUR = T_OUR*m; % rxSig_SniffAttack_Countermeasure2(M): 攻擊者收到的訊號

%sniffMessage_Countermeasure2 = rxSig_SniffAttack_Countermeasure2(M,:) - m(M,:) ; %攻擊者消掉自己的訊息
%sniffMessage_Countermeasure2 = rxSig_SniffAttack_Countermeasure2(M,:) - (m(M,:)*mask_bar(M)) ; %攻擊者消掉自己的訊息
%sniffMessage_Countermeasure2 = sniffMessage_Countermeasure2 / Var_x(2); 
%sniffMessage_Countermeasure2 = sniffMessage_Countermeasure2 / (mask_bar(M)/Var_x(2)); 
%sniffMessage_Countermeasure2 = awgn(sniffMessage_Countermeasure2,snr,'measured');
%sniffMessage_Label_Countermeasure2 = pskdemod(sniffMessage_Countermeasure2,M_QPSK,pi/M_QPSK);
%[errSym_QPSK_Countermeasure2, SER_QPSK_Countermeasure2] = symerr(sniffMessage_Label_Countermeasure2, x2);  




for snr_i = 1:length(Eb_N0_dB)
    UE2_message = awgn(rxSig_OUR(2,:),Eb_N0_dB(snr_i),'measured');
    
    sniffMessage_OUR_AWGN = awgn(rxSig_OUR(M,:),Eb_N0_dB(snr_i),'measured');
    sniffMessage_OUR = sniffMessage_OUR_AWGN - m(M,:);
    
    if MOD== 16 || MOD == 32
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2');
        
        sniffLabel_OUR_strong = dvbsapskdemod(sniffMessage_OUR,MOD,'s2');
        StrongSig = dvbsapskmod(sniffLabel_OUR_strong,MOD,'s2');
        remain_sniffMessage_OUR = sniffMessage_OUR - StrongSig;
        sniffLabel_OUR_weak = dvbsapskdemod(remain_sniffMessage_OUR,MOD,'s2');
    else
        UE2_message_Label = pskdemod(UE2_message,MOD,pi/MOD);
        
        sniffLabel_OUR_strong = pskdemod(sniffMessage_OUR,MOD,pi/MOD);
        StrongSig = pskmod(sniffLabel_OUR_strong,MOD,pi/MOD);
        remain_sniffMessage_OUR = sniffMessage_OUR - StrongSig;
        sniffLabel_OUR_weak = pskdemod(remain_sniffMessage_OUR,MOD,pi/MOD);
        
    end
    [errSym(snr_i), OUR_SER(snr_i)] = symerr(UE2_message_Label, x(2,:));  
    
    [OUR_errSym_Eve_strong(snr_i), OUR_SER_Eve_strong(snr_i)] = symerr(sniffLabel_OUR_strong, x(2,:));  
    [OUR_errSym_Eve_weak(snr_i), OUR_SER_Eve_weak(snr_i)] = symerr(sniffLabel_OUR_weak, x(1,:));  
    
end

PlotConstellation(UE2_message, txSig(2,:), 'UE_2 receive message with our method', MOD, spacing);
PlotConstellation(sniffMessage_OUR, txSig(2,:), 'Eavesdropping UE_2 message with our method', MOD, spacing);


%% Plot SER

MODName = '';
if MOD == 4
    MODName = 'QPSK';
elseif MOD == 8
    MODName = '8PSK';
elseif MOD == 16
    MODName = '16APSK';
elseif MOD == 32
    MODName = '32APSK';
end

%PlotSER(Eb_N0_dB, SCF_SER_Eve, OUR_SER_Eve, MODName);

PlotSER_multi(Eb_N0_dB, SCF_SER_Eve_weak, SCF_SER_Eve_strong, OUR_SER_Eve_weak, OUR_SER_Eve_strong, MODName)



function PlotConstellation(receivedSym, referenceSym, title_name, MOD, spacing)

    %figure();
    %plot(receivedSym,'b*'); hold on;
    %plot(referenceSym,'r+','LineWidth',2,'MarkerSize',10); hold on;
    
    figure();
    Color = ["#FF0000", "#0000FF", "#FFFF00", "#00FF00", "#FFBB00", "#7700FF", "#FF00FF", "#00FFFF",...
             "#FF8888", "#9999FF", "#BBBB00", "#00AA00", "#AA7700", "#B088FF", "#990099", "#00AAAA",...
             "#FF0088", "#CCBBFF", "#88AA00", "#008800", "#BB5500", "#FFDD55", "#009FCC", "#00FF99",...
             "#FFB7DD", "#9955FF", "#FFFFBB", "#FF8888", "#A42D00", "#AAFFEE", "#8C0044", "#CCFF33"];
    
    i = 1;
    plot(receivedSym(1+spacing*(i-1):spacing*i),'*','color',Color(i)); hold on;
    if MOD == 16 || MOD == 32
        plot(dvbsapskmod(i-1,MOD,'s2'),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    else
        plot(pskmod(i-1,MOD,pi/MOD),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    end
    
    temp = [1:spacing:spacing*MOD];
    for i = 1:length(temp)
        plot(receivedSym(1+spacing*(i-1):spacing*i),'*','color',Color(i)); hold on;
    end
    
    for i = 1:length(temp)
        if MOD== 16 || MOD == 32
            temp2 = complex(real(dvbsapskmod(i-1,MOD,'s2')),imag(dvbsapskmod(i-1,MOD,'s2')));
        else
            temp2 = complex(real(pskmod(i-1,MOD,pi/MOD)),imag(pskmod(i-1,MOD,pi/MOD)));
        end
        plot(temp2,'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    end

    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    xlabel('In-Phase');
    ylabel('Quadrature');
    xlim([-2,2]);
    ylim([-2,2]);
    legend('Received constellation','Reference constellation');
    title(title_name);
    saveas(gcf,strcat(title_name,'.png'));
end

function PlotSER_multi(Eb_N0_dB, Eve_SER_weak_SCF, Eve_SER_strong_SCF, Eve_SER_weak_OUR, Eve_SER_strong_OUR, title_name)
    figure();
    %figure('position',[300, 300, 800, 600]);
    semilogy(Eb_N0_dB,Eve_SER_weak_OUR,'r-o','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_strong_OUR,'r-+','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_weak_SCF,'b-s','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_strong_SCF,'b-d','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_weak_SCF,'g-^','LineWidth',2.5,'MarkerSize',12,'Color',[0 170/255 0]);hold on;
    semilogy(Eb_N0_dB,Eve_SER_strong_SCF,'g-v','LineWidth',2.5,'MarkerSize',12,'Color',[0 170/255 0]);hold on;
    
    
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('E_b/N_0 (dB)', 'FontSize', 18, 'fontweight','bold','Interpreter','tex');
    ylabel('Symbol Error Rate', 'FontSize', 18, 'fontweight','bold','Interpreter','tex');
    legend('sniffing UT1 w/ SSC','sniffing UT2 w/ SSC','sniffing UT1 w/ SCF','sniffing UT2 w/ SCF','sniffing UT1','sniffing UT2','Location','southwest');
    title(title_name, 'FontSize', 18, 'fontweight','bold','Interpreter','tex');
    title_name = strcat('SER_',title_name);
    saveas(gcf,strcat(title_name,'.png'));
    
    
end

function PlotSINR_Multi(Eb_N0_dB, SINR_weak, SINR_strong, SINR_Eve, SINR_weak_SA, SINR_strong_SA, SINR_weak_Eve_SA, SINR_strong_Eve_SA)
    figure();
    
    MarkerSize = 12;
    
    plot(Eb_N0_dB,SINR_weak,'b-^','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    plot(Eb_N0_dB,SINR_strong,'b-v','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    plot(Eb_N0_dB,SINR_Eve,'b-*','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    
    plot(Eb_N0_dB,SINR_weak_SA,'r-s','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    plot(Eb_N0_dB,SINR_strong_SA,'r-d','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    plot(Eb_N0_dB,SINR_weak_Eve_SA,'r-o','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    plot(Eb_N0_dB,SINR_strong_Eve_SA,'r-+','LineWidth',2.5,'MarkerSize',MarkerSize);hold on;
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('E_b/N_0 (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    ylabel('Received SINR (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    xlim([-3,20]);
    ylim([-10,30]);
    legend('UT1 w/o SA','UT2 w/o SA','Eve w/o SA','UT1 w/ SA','UT2 w/ SA','Eve sniffing UT1','Eve sniffing UT2','Location','northwest','FontSize',13);
    %legend('Eve with sniffing attack','Eve without sniffing attack','UT with sniffing attack','UT without sniffing attack','Location','northwest');
    saveas(gcf,'SINR_multi.png');
end