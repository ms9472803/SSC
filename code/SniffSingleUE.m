clear variables;
close all;
load('SatCH.mat');
load('NoiseParam.mat');

load('tmp.mat');
Random_SER_Eve = OUR_SER_Eve_M;

TransmitterPower_dB = 60; % dBm
%LNAGain_dB = 35; % dBm
PowerVar = 10^(TransmitterPower_dB/10);
%PowerVar = 10^((TransmitterPower_dB+LNAGain_dB)/10);
NumTransmitAntennas = 3;
NumReceiveAntennas = 3;


%% SNR calculation

Eb_N0_dB = cell2mat(Eb_N0_dB);
Eb_N0_dB_MAX = Eb_N0_dB;
Eb_N0_dB = Eb_N0_dB_MAX-30:2:Eb_N0_dB_MAX-10; % Es/N0 in dB

RcvrPower_dB = cell2mat(RcvrPower_dB);
Eb_N0 = 10.^(Eb_N0_dB./10);
RcvrPower = 10.^(RcvrPower_dB./10); 
NoiseVar = RcvrPower./Eb_N0; 
snr = 20;

%% H

%H = complex(randn(NumTransmitAntennas,NumReceiveAntennas), randn(NumTransmitAntennas,NumReceiveAntennas)); % random genereate downlink CSI
H = SatCH(randi(length(SatCH),NumReceiveAntennas,NumTransmitAntennas)); % random from Sat Channel
%imperfectH = H + err;

%% Simulation
% Number of Monte-Carlo iterations
NumIter = 1;

for it = 1:NumIter
disp(it)
%% mask CSI 
theta_max = [45 22.5 15 11.25];
candidate = [45, 90, 135, 180, 225, 270, 315, 360];
%theta = [50 randi([90, 360]) randi([90, 360])];
%theta = [50 candidate(randi([1,8])) candidate(randi([1,8]))];
theta = [50 100 200];
alpha = [complex(cosd(theta(1)),sind(theta(1))) complex(cosd(theta(2)),sind(theta(2))) complex(cosd(theta(3)),sind(theta(3)))];
Amplitude = [1 1 1];
mask = Amplitude .* alpha;

mask_CSI = zeros(1,3); % initialization ; BS receive from UE's feedback

for i = 1:3
    mask_CSI(i,:) = H(i,:)*mask(i);
end

precoding_mask_ByCh = H*pinv(mask_CSI); % H*H^+ = diagonal matrix: [bar(alpha1); bar(alpha2); bar(alpha3)]

%% mask pilot symbol (SCF)
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
symbolNum = 30000; %一共有幾個符號
M_QPSK = 4; % QPSK order
M_8PSK = 8;
M_16APSK = 16; 
M_32APSK = 32;
M_64APSK = 64;
MOD = M_32APSK;

x = [];
txSig = [];

spacing = symbolNum/MOD;
temp = [1:spacing:symbolNum];
for i = 1:NumReceiveAntennas
   x(i,:) = randi([0 MOD-1],1,symbolNum); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
   
   for k = 0:length(temp)-1
       x(i,1+spacing*k:spacing*(k+1)) = k;
   end

   
   if MOD== 16 || MOD == 32
       txSig(i,:) = dvbsapskmod(x(i,:),MOD,'s2x'); 
   elseif MOD == 64
       txSig(i,:) = dvbsapskmod(x(i,:),MOD,'s2x','132/180'); 
   else
       txSig(i,:) = pskmod(x(i,:),MOD,pi/MOD); 
   end

end

precoding = pinv(H);
ChByPrecoding = H*precoding;

m = txSig;
%[txSig1;txSig2;txSig3];
rxSig = ChByPrecoding*m; % Y = HX + n = CPm + n

%% Sniffing attack (sniff single UE)

M = 3;
Var_x = [0 1 1]; 
F = H;
F(3,:) = 0;

%竊聽者估計自己的CSI 還有竊聽UE_2的CSI
eve_CSI = zeros(3,3);
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

syms f(xk);
f = xk-1+totalPower/(NoiseVar_Eve(11)*(norm(pinvH,'fro')^2)+xk*norm(pinvH(:,M),'fro')^2);
xk = solve(f);

SINR = [];
for i = 1:length(NoiseVar_Eve)
    SINR(i) = power_i*(Var_x(2)^2)/(NoiseVar_Eve(i)/2 + power_i);
end
SINR_Eve_sniffing_attack = 10*log10(SINR);

PlotSINR(Eb_N0_dB, Eb_N0_dB, SINR_UT_sniffing_attack, Eb_N0_dB-5, SINR_Eve_sniffing_attack);

 
ChByForgedCSIPrecoding = H*pinvF;
T = ChByForgedCSIPrecoding;
rxSig_SniffAttack = ChByForgedCSIPrecoding*m;
rxSig_SniffAttack_AWGN = awgn(rxSig_SniffAttack,snr,'measured');

sniffMessage = (rxSig_SniffAttack_AWGN(M,:)-m(M,:))/Var_x(2);
if MOD== 16 || MOD == 32 || MOD == 64
    sniffMessage_Label = dvbsapskdemod(sniffMessage,MOD,'s2x');
elseif MOD == 64
    sniffMessage_Label = dvbsapskdemod(sniffMessage,MOD,'s2x','132/180');
else
    sniffMessage_Label = pskdemod(sniffMessage,MOD,pi/MOD);
end
[errSym_QPSK, SER_QPSK] = symerr(sniffMessage_Label, x(2,:));
fprintf("SER QPSK Sniff attack: %d/%d(%f%%)\n", errSym_QPSK, length(sniffMessage_Label), SER_QPSK);

%PlotConstellation(sniffMessage, txSig(2,:), 'Sniffing attack Eavesdropping UE_2 message', MOD, spacing);

%% PROOF OF THEOREM 1
%{
% T*F*F'(3,1) = Sum1 + Sum2 = Sum1 + Sum3
l = 1;
M = 3;
Sum1 = 0;
for j=1:M-1
    Sum1 = Sum1 + Var_x(j) * sum(H(j,:).'.*H(l,:)');
end

Sum2 = 0;
for i=1:3
    Sum2 = Sum2 + F(3,i).'.*H(l,i)';
    %Sum2 = Sum2 + F(3,i);
    %F(3,i)
end

F__ = H;
F__(3,:) = 0;
Sum3 = 0;
for i=1:3
    temp = 0;
    for k=1:M-1
        temp = temp + Var_x(k)*H(k,i);
    end
   Sum3 = Sum3 + (H(M,i)-temp).'.*H(1,i)'; 
   %Sum3 = Sum3 + (H(M,i)-temp);
   F__(3,i) = (H(M,i)-temp);
end
%}

%% Countermeasure for sniffing attack, mask pilot symbol (SCF)

% satellite sends the same pilot symbol with same mask for each UT to estimate CSI
Y_mask_pilot = zeros(3,3); % received pilot
mask_SCF = complex(cosd(50),sind(50));
for i = 1:3
    Y_mask_pilot(i,:) = H(i,:)*(pilot*mask_SCF); % Y=HX
end
CSI_mask_pilot = Y_mask_pilot/pilot; % UE estimate CSI_mask_pilot and feedback to Sat

%竊聽者估計自己的CSI 還有竊聽UE_2的CSI
eve_CSI_mask_pilot = zeros(3,3);
eve_CSI_mask_pilot(2,:) =  CSI_mask_pilot(2,:);
eve_CSI_mask_pilot(3,:) = CSI_mask_pilot(3,:); 

F_countermeasure = CSI_mask_pilot; % UE feedback estimated CSI
F_countermeasure(3,:) = 0;
% 用公式計算F
for i = 1:3
    temp = 0;
    for k = 1:M-1
        temp = temp + Var_x(k)*eve_CSI_mask_pilot(k,i);
    end
   F_countermeasure(3,i) = (eve_CSI_mask_pilot(M,i)-temp);
end

BS_recover_F_countermeasure = zeros(3,3);
for i = 1:3
    BS_recover_F_countermeasure(i,:) = F_countermeasure(i,:)/mask_SCF; % BS recover CSI for AMC and precoding (using: CSI_mask_pilot=mask*H)
end
ChByCountermeasurePrecoding = H*pinv(BS_recover_F_countermeasure); % H*H^+ 

rxSig_SniffAttack_Countermeasure = ChByCountermeasurePrecoding*m;
%sniffMessage_Countermeasure = (rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/Var_x(2);
%sniffMessage_Countermeasure = awgn(sniffMessage_Countermeasure,snr,'measured');
%sniffMessage_Label_Countermeasure = pskdemod(sniffMessage_Countermeasure,MOD,pi/MOD);
%[errSym_QPSK_Countermeasure, SER_QPSK_Countermeasure] = symerr(sniffMessage_Label_Countermeasure, x(2,:));   
% [a,b] = symerr(pskdemod(awgn((rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/((Var_x(2) * alpha3)*-1),snr,'measured'),M_QPSK,pi/M_QPSK),x2)

for snr_i = 1:length(Eb_N0_dB)
    UE2_message = awgn(rxSig_SniffAttack_Countermeasure(2,:),Eb_N0_dB(snr_i),'measured');
    
    sniffMessage_Countermeasure_ = awgn(rxSig_SniffAttack_Countermeasure(M,:),Eb_N0_dB(snr_i),'measured');
    sniffMessage_Countermeasure_ = (sniffMessage_Countermeasure_ - m(M,:)) / Var_x(2);
    
    %sniffMessage_Countermeasure_ = awgn(sniffMessage_Countermeasure,Eb_N0_dB(snr_i),'measured');
    
    if MOD== 16 || MOD == 32 
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x');
        sniffMessage_Label_Countermeasure = dvbsapskdemod(sniffMessage_Countermeasure_,MOD,'s2x');
    elseif MOD == 64
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x','132/180');
        sniffMessage_Label_Countermeasure = dvbsapskdemod(sniffMessage_Countermeasure_,MOD,'s2x','132/180');
    else
        UE2_message_Label = pskdemod(UE2_message,MOD,pi/MOD);
        sniffMessage_Label_Countermeasure = pskdemod(sniffMessage_Countermeasure_,MOD,pi/MOD);
    end
    
    [SCF_errSym(snr_i,it), SCF_SER(snr_i,it)] = symerr(UE2_message_Label, x(2,:));  
    [SCF_errSym_Eve(snr_i,it), SCF_SER_Eve(snr_i,it)] = symerr(sniffMessage_Label_Countermeasure, x(2,:));  
    
end

%PlotConstellation(UE2_message, txSig(2,:), 'UE_2 receive message with SCF', MOD, spacing);
%PlotConstellation(sniffMessage_Countermeasure_, txSig(2,:), 'Eavesdropping UE_2 message with SCF', MOD, spacing);


%% Countermeasure for sniffing attack, mask CSI feedback (SSC)
%竊聽者估計自己的CSI 還有竊聽UE_2的CSI
eve_CSI_mask_CSI = zeros(3,3);
eve_CSI_mask_CSI(2,:) =  mask_CSI(2,:);
eve_CSI_mask_CSI(3,:) = mask_CSI(3,:); 

F_countermeasure2 = mask_CSI; % UE feedback estimated CSI 再做mask

% 攻擊者用公式計算F
for i = 1:3
    temp = 0;
    for k = 1:M-1
        temp = temp + Var_x(k)*eve_CSI_mask_CSI(k,i);
    end
   F_countermeasure2(3,i) = (eve_CSI_mask_CSI(M,i)-temp);
end

BS_recover_F_countermeasure2 = F_countermeasure2;
for i = 1:3
    BS_recover_F_countermeasure2(i,:) = F_countermeasure2(i,:)/mask(i);
end

ChByCountermeasure2Precoding = H*pinv(BS_recover_F_countermeasure2)
%(Var_x(2)*alpha(2))/alpha(3)
%Var_x(2) * alpha(3)'

rxSig_SniffAttack_Countermeasure2 = ChByCountermeasure2Precoding*m; % rxSig_SniffAttack_Countermeasure2(M): 攻擊者收到的訊號

sniffMessage_Countermeasure2 = rxSig_SniffAttack_Countermeasure2(M,:) - m(M,:) ; %攻擊者消掉自己的訊息
%sniffMessage_Countermeasure2 = rxSig_SniffAttack_Countermeasure2(M,:) - (m(M,:)*mask_bar(M)) ; %攻擊者消掉自己的訊息
sniffMessage_Countermeasure2 = sniffMessage_Countermeasure2 / Var_x(2); 
%sniffMessage_Countermeasure2 = sniffMessage_Countermeasure2 / (mask_bar(M)/Var_x(2)); 
%sniffMessage_Countermeasure2 = awgn(sniffMessage_Countermeasure2,snr,'measured');
%sniffMessage_Label_Countermeasure2 = pskdemod(sniffMessage_Countermeasure2,M_QPSK,pi/M_QPSK);
%[errSym_QPSK_Countermeasure2, SER_QPSK_Countermeasure2] = symerr(sniffMessage_Label_Countermeasure2, x2);  




for snr_i = 1:length(Eb_N0_dB)
    
    UE2_message = awgn(rxSig_SniffAttack_Countermeasure2(2,:),Eb_N0_dB(snr_i),'measured');
    
    sniffMessage_Countermeasure2_ = awgn(rxSig_SniffAttack_Countermeasure2(M,:),Eb_N0_dB(snr_i),'measured');
    %sniffMessage_Countermeasure2_ = (sniffMessage_Countermeasure2_ - m(M,:)) / Var_x(2);
    
    sniffMessage_Countermeasure2_ = awgn(sniffMessage_Countermeasure2,Eb_N0_dB(snr_i),'measured');
    
    if MOD== 16 || MOD == 32 
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x');
        sniffMessage_Label_Countermeasure2 = dvbsapskdemod(sniffMessage_Countermeasure2_,MOD,'s2x');
    elseif MOD == 64
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x','132/180');
        sniffMessage_Label_Countermeasure2 = dvbsapskdemod(sniffMessage_Countermeasure2_,MOD,'s2x','132/180');
    else
        UE2_message_Label = pskdemod(UE2_message,MOD,pi/MOD);
        sniffMessage_Label_Countermeasure2 = pskdemod(sniffMessage_Countermeasure2_,MOD,pi/MOD);
    end
    
    [errSym(snr_i,it), OUR_SER(snr_i,it)] = symerr(UE2_message_Label, x(2,:));  
    [errSym_Eve(snr_i,it), OUR_SER_Eve(snr_i,it)] = symerr(sniffMessage_Label_Countermeasure2, x(2,:));  
    
end

PlotConstellation(UE2_message, txSig(2,:), 'UE_2 receive message with our method', MOD, spacing);
PlotConstellation(sniffMessage_Countermeasure2_, txSig(2,:), 'Eavesdropping UE_2 message with our method', MOD, spacing);



%{
%% Countermeasure for sniffing attack, mask pilot symbol and mask CSI (ver 2.0)

% satellite sends the same pilot symbol with same mask for each UT to estimate CSI
Y_mask_pilot = zeros(3,3); % received pilot
mask_SAT = complex(cosd(50),sind(50));
for i = 1:3
    Y_mask_pilot(i,:) = H(i,:)*(pilot*mask_SAT); % Y=HX
end
CSI_mask_pilot = Y_mask_pilot/pilot; % UE estimate CSI_mask_pilot

for i = 1:3
    CSI_mask_pilot_mask_CSI(i,:) = CSI_mask_pilot(i,:)*mask(i);
end

%竊聽者估計自己的CSI 還有竊聽UE_2的CSI
eve_CSI_mask_pilot_mask_CSI = zeros(3,3);
eve_CSI_mask_pilot_mask_CSI(2,:) =  CSI_mask_pilot_mask_CSI(2,:);
eve_CSI_mask_pilot_mask_CSI(3,:) = CSI_mask_pilot_mask_CSI(3,:); 

F_countermeasure = CSI_mask_pilot_mask_CSI; % UE feedback estimated CSI
F_countermeasure(3,:) = 0;
% 用公式計算F
for i = 1:3
    temp = 0;
    for k = 1:M-1
        temp = temp + Var_x(k)*eve_CSI_mask_pilot_mask_CSI(k,i);
    end
   F_countermeasure(3,i) = (eve_CSI_mask_pilot_mask_CSI(M,i)-temp);
end

BS_recover_F_countermeasure = zeros(3,3);
BS_recover_F_countermeasure2 = zeros(3,3);
for i = 1:3
    BS_recover_F_countermeasure(i,:) = F_countermeasure(i,:) /mask_SAT; % BS recover CSI for AMC and precoding (using: CSI_mask_pilot=mask*H)
    BS_recover_F_countermeasure2(i,:) = (F_countermeasure(i,:).*F_countermeasure(i,:)) /mask_SAT; % BS recover CSI for AMC and precoding (using: CSI_mask_pilot=mask*H)
end
ChByCountermeasurePrecoding = H*pinv(BS_recover_F_countermeasure); % H*H^+ 

rxSig_SniffAttack_Countermeasure = ChByCountermeasurePrecoding*m;
%sniffMessage_Countermeasure = (rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/Var_x(2);
%sniffMessage_Countermeasure = awgn(sniffMessage_Countermeasure,snr,'measured');
%sniffMessage_Label_Countermeasure = pskdemod(sniffMessage_Countermeasure,MOD,pi/MOD);
%[errSym_QPSK_Countermeasure, SER_QPSK_Countermeasure] = symerr(sniffMessage_Label_Countermeasure, x(2,:));   
% [a,b] = symerr(pskdemod(awgn((rxSig_SniffAttack_Countermeasure(M,:)-m(M,:))/((Var_x(2) * alpha3)*-1),snr,'measured'),M_QPSK,pi/M_QPSK),x2)

for snr_i = 1:length(Eb_N0_dB)
    UE2_message = awgn(rxSig_SniffAttack_Countermeasure(2,:),Eb_N0_dB(snr_i),'measured');
    
    sniffMessage_Countermeasure_ = awgn(rxSig_SniffAttack_Countermeasure(M,:),Eb_N0_dB(snr_i),'measured');
    sniffMessage_Countermeasure_ = (sniffMessage_Countermeasure_ - m(M,:)) / Var_x(2);
    
    %sniffMessage_Countermeasure_ = awgn(sniffMessage_Countermeasure,Eb_N0_dB(snr_i),'measured');
    
    if MOD== 16 || MOD == 32 
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x');
        sniffMessage_Label_Countermeasure = dvbsapskdemod(sniffMessage_Countermeasure_,MOD,'s2x');
    elseif MOD == 64
        UE2_message_Label = dvbsapskdemod(UE2_message,MOD,'s2x','132/180');
        sniffMessage_Label_Countermeasure = dvbsapskdemod(sniffMessage_Countermeasure_,MOD,'s2x','132/180');
    else
        UE2_message_Label = pskdemod(UE2_message,MOD,pi/MOD);
        sniffMessage_Label_Countermeasure = pskdemod(sniffMessage_Countermeasure_,MOD,pi/MOD);
    end
    
    [SCF_errSym(snr_i,it), SCF_SER(snr_i,it)] = symerr(UE2_message_Label, x(2,:));  
    [SCF_errSym_Eve(snr_i,it), SCF_SER_Eve(snr_i,it)] = symerr(sniffMessage_Label_Countermeasure, x(2,:));  
    
end
%}

%PlotConstellation(UE2_message, txSig(2,:), 'UE_2 receive message with SCF', MOD, spacing);
%PlotConstellation(sniffMessage_Countermeasure_, txSig(2,:), 'Eavesdropping UE_2 message with SCF', MOD, spacing);

end


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
elseif MOD == 64
    MODName = '64APSK';
end


SCF_SER_Eve_M = mean(SCF_SER_Eve,2).';
OUR_SER_Eve_M = mean(OUR_SER_Eve,2).';
%Random_SER_Eve = mean(Random_SER_Eve,2).';

PlotSER(Eb_N0_dB, SCF_SER_Eve_M, OUR_SER_Eve_M, SCF_SER_Eve_M, MODName);
%PlotSER2(Eb_N0_dB, Random_SER_Eve, OUR_SER_Eve_M, MODName);

function PlotConstellation(receivedSym, referenceSym, title_name, MOD, spacing)

    %figure();
    %plot(receivedSym,'b*'); hold on;
    %plot(referenceSym,'r+','LineWidth',2,'MarkerSize',10); hold on;
    
    figure();
    Color = ["#000000"];
    for i = 2:MOD
       Color = [Color, "#000000"]; 
    end
    Color(1:32) = ["#FF0000", "#0000FF", "#FFFF00", "#00FF00", "#FFBB00", "#7700FF", "#FF00FF", "#00FFFF",...
             "#FF8888", "#9999FF", "#BBBB00", "#00AA00", "#AA7700", "#B088FF", "#990099", "#00AAAA",...
             "#FF0088", "#CCBBFF", "#88AA00", "#008800", "#BB5500", "#FFDD55", "#009FCC", "#00FF99",...
             "#FFB7DD", "#9955FF", "#FFFFBB", "#FF8888", "#A42D00", "#AAFFEE", "#8C0044", "#CCFF33"];
    
    i = 1;
    plot(receivedSym(1+spacing*(i-1):spacing*i),'*','color',Color(i),'MarkerSize',10); hold on;
    if MOD == 16 || MOD == 32
        plot(dvbsapskmod(i-1,MOD,'s2x'),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    elseif MOD == 64
        plot(dvbsapskmod(i-1,MOD,'s2x','132/180'),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    else
        plot(pskmod(i-1,MOD,pi/MOD),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    end
    
    temp = [1:spacing:spacing*MOD];
    
    for i = 1:length(temp)
        plot(receivedSym(1+spacing*(i-1):spacing*i),'*','color',Color(i)); hold on;
    end
    
    %{
    i=10;
    if MOD == 16 || MOD == 32
        plot(dvbsapskmod(i-1,MOD,'s2'),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
        disp(dvbsapskmod(i-1,MOD,'s2'));
    else
        plot(pskmod(i-1,MOD,pi/MOD),'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    end
    %}
    
    
    for i = 1:length(temp)
        if MOD== 16 || MOD == 32 
            temp2 = complex(real(dvbsapskmod(i-1,MOD,'s2x')),imag(dvbsapskmod(i-1,MOD,'s2x')));
        elseif MOD == 64
            temp2 = complex(real(dvbsapskmod(i-1,MOD,'s2x','132/180')),imag(dvbsapskmod(i-1,MOD,'s2x','132/180')));
        else
            temp2 = complex(real(pskmod(i-1,MOD,pi/MOD)),imag(pskmod(i-1,MOD,pi/MOD)));
        end
        plot(temp2,'d','MarkerFaceColor',Color(i),'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10); hold on;
    end
    
    
    
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('In-Phase', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    ylabel('Quadrature', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    xlim([-2.3,2.3]);
    ylim([-2.3,2.3]);
    %legend('Received signal','Reference constellation');
    legend('Received signal','32APSK constellation');
    %title(title_name);
    saveas(gcf,strcat(title_name,'.png'));
    
    %RGB = imread(strcat(title_name,'.png'));
    %I = rgb2gray(RGB);
    %figure();
    %imshow(I);
    
end

function PlotSER(Eb_N0_dB, Eve_SER_SCF, Eve_SER_OUR, Eve_SER_NO, title_name)
    figure();
    
    semilogy(Eb_N0_dB,Eve_SER_OUR,'r-*','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_SCF,'b-o','LineWidth',2.5,'MarkerSize',12);hold on;
    %semilogy(Eb_N0_dB,Eve_SER_NO,'y-+','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_NO,'-+','LineWidth',2.5,'MarkerSize',12,'Color',[0 170/255 0]);hold on;
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('E_b/N_0 (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    ylabel('Symbol Error Rate', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    legend('Eve w/ SSC','Eve w/ SCF','Eve w/o countermeasure','Location','southwest');
    title(title_name, 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    title_name = strcat('SER_',title_name);
    saveas(gcf,strcat(title_name,'.png'));
end

function PlotSER2(Eb_N0_dB, Eve_SER_Random, Eve_SER_OUR, title_name)
    figure();
    semilogy(Eb_N0_dB,Eve_SER_OUR,'r-o','LineWidth',2.5,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,Eve_SER_Random,'b-*','LineWidth',2.5,'MarkerSize',12);hold on;
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('E_b/N_0 (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    ylabel('Symbol Error Rate', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    legend('Eve with SSC','Eve with random','Location','southwest');
    title(title_name, 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    title_name = strcat('SER_',title_name);
    saveas(gcf,strcat(title_name,'.png'));
end

function PlotSINR(Eb_N0_dB, SINR_UT2, SINR_UT2_SA, SINR_Eve, SINR_Eve_SA)
    figure();
    plot(Eb_N0_dB,SINR_UT2,'b-o','LineWidth',2.5,'MarkerSize',12);hold on;
    plot(Eb_N0_dB,SINR_Eve,'b-h','LineWidth',2.5,'MarkerSize',12);hold on;
    plot(Eb_N0_dB,SINR_UT2_SA,'r-s','LineWidth',2.5,'MarkerSize',12);hold on;
    plot(Eb_N0_dB,SINR_Eve_SA,'r-*','LineWidth',2.5,'MarkerSize',12);hold on;
    
    set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
    xlabel('E_b/N_0 (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    ylabel('Received SINR (dB)', 'FontSize', 22, 'fontweight','bold','Interpreter','tex');
    xlim([-5,20]);
    ylim([-10,30]);
    legend('UT w/o SA','Eve w/o SA','UT w/ SA', 'Eve w/ SA','Location','northwest');
    saveas(gcf,'SINR.png');
end

