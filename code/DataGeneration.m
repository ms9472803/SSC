%% Training Data Generation 
%
% This script is created to generate training and validation data for the deep
% learning model in a single-user system. 

% The training data and the validation data is collected for a single
% subcarrier selected based on a pre-defined metric. The transmitter sends
% packets to the receiver, where each packet contains one pilot symbol and 
% one data symbol. Data symbols can be interleaved in
% the pilot sequence. 

%% Clear workspace

clear variables;
close all;
load('EAngle.mat');
load('SatChannelParam.mat');


%% OFDM system parameters

NumSC = 64; % Number of subcarriers
NumPilot = 64; % Number of pilot subcarriers
PilotSpacing = NumSC/NumPilot;
NumPilotSym = 1;
NumDataSym = 1;
NumSym = NumPilotSym + NumDataSym;

%% Channel generation

NumPath = 20;
LengthCP = 16; % The length of the cyclic prefix
% The channel matrix generated using the 3GPP TR38.901 channel model of the
% writer's own implementation, which is saved and loaded:
load('SavedChan.mat');  
load('NoiseParam.mat');

% One can replace the 3GPP channel with the narrowband Rayleigh fading channel:  
% h = 1/sqrt(2)*complex(randn(NumPath,1),randn(NumPath,1)); 

H = fft(h,NumSC,1); 
AbsH = abs(H); % get different sub-carriers channel gain

%% Transmitter Model Spec

TransmitterPower_dB = 60; % dBm
LNAGain_dB = 35; % dBm
PowerVar = 10^(TransmitterPower_dB/10);
%PowerVar = 10^((TransmitterPower_dB+LNAGain_dB)/10);

%% SNR calculation

Eb_N0_dB = cell2mat(Eb_N0_dB);
Eb_N0_dB_MAX = Eb_N0_dB;
Eb_N0_dB = Eb_N0_dB_MAX-20:2:Eb_N0_dB_MAX; % Es/N0 in dB

RcvrPower_dB = cell2mat(RcvrPower_dB);
Eb_N0 = 10.^(Eb_N0_dB./10);
RcvrPower = 10.^(RcvrPower_dB./10);
NoiseVar = RcvrPower./Eb_N0;
% variance = 10^(-snr/10);

%% Subcarrier selection

% This is the subcarrier selected for the saved channel h, which can be
% replaced by the following process.
%idxSC = 7;
idxSC = 26;

% Select the subcarrier randomly whose gain is above the median value
%MedianGain = median(abs(H).^2);
%[PossibleSC,~] = find(logical(abs(H).^2 >= MedianGain) == 1);
%idxSC = PossibleSC(randi(length(PossibleSC)));
%disp(idxSC);

%% Simulation

% Size of dataset to be defined
NumPacket = 45000; % Number of packets per LEO track
%NumPacket = length(EAngle{2});

% Pilot symbols - Fixed during the whole transmission
FixedPilot = sqrt(PowerVar/2)*complex(sign(rand(1,NumPilot)-0.5),sign(rand(1,NumPilot)-0.5)); 

CSIRS = [0.7071+0.7071i, -0.7071+0.7071i, -0.7071-0.7071i, 0.7071-0.7071i];
theta_rm = [45, 45+90,  45+180, 45+270];
knownSeq1 = zeros(1,NumPacket);

for j = 1:length(knownSeq1)
    knownSeq1(j) = CSIRS(1);
    %knownSeq1(j) = CSIRS(randi([1,4],1,1));
end

theta = 50;
knownSeq2 = knownSeq1.*complex(cosd(theta),sind(theta)); 
%{
knownSeq2 = zeros(1,NumPacket);
for j = 1:length(knownSeq2)
    for i = 1:length(CSIRS)
       if knownSeq1(j) == CSIRS(i)
           theta1 = theta_rm(i);
           theta2 = theta1 + theta;
           knownSeq2(j) = complex(cosd(theta2),sind(theta2));
           break;
       end
    end
end
%}


knownSeq1 = sqrt(PowerVar/2) * knownSeq1;
knownSeq2 = sqrt(PowerVar/2) * knownSeq2;
% Same pilot sequences used in all packets
%FixedPilotAll = repmat(FixedPilot,1,1,NumPacket); 
FixedPilotAll_nonMask = repmat(FixedPilot,1,1,NumPacket); 
FixedPilotAll_Mask = repmat(FixedPilot,1,1,NumPacket); 

for i = 1:length(knownSeq2)
    FixedPilotAll_nonMask(:,:,i) = knownSeq1(i);
    FixedPilotAll_Mask(:,:,i) = knownSeq2(i);
end

% Pilot symbol (can be interleaved with random data symbols)

PilotSym = FixedPilotAll_nonMask;

% Data symbol QPSK

M_QPSK = 4;
DataLabel_QPSK = randi([0 M_QPSK-1],NumSC,NumPacket);
DataLabel_idxSc_QPSK = DataLabel_QPSK(idxSC,:).';
DataSym_QPSK(1,:,:) = sqrt(PowerVar/2)*pskmod(DataLabel_QPSK,M_QPSK,pi/M_QPSK);

% Data symbol 8PSK
M_8PSK = 8;
DataLabel_8PSK = randi([0 M_8PSK-1],NumSC,NumPacket);
DataLabel_idxSc_8PSK = DataLabel_8PSK(idxSC,:).';
DataSym_8PSK(1,:,:) = sqrt(PowerVar/2)*pskmod(DataLabel_8PSK,M_8PSK,pi/M_8PSK);

% Data symbol 16APSK    
M_16APSK = 16; 
DataLabel_16APSK = randi([0 M_16APSK-1],NumSC,NumPacket); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
DataLabel_idxSc_16APSK = DataLabel_16APSK(idxSC,:).';
DataSym_16APSK(1,:,:) = sqrt(PowerVar/2)*dvbsapskmod(DataLabel_16APSK,M_16APSK,'s2'); % Apply APSK modulation to the data.

% Data symbol 32APSK    
M_32APSK = 32;
DataLabel_32APSK = randi([0 M_32APSK-1],NumSC,NumPacket); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
DataLabel_idxSc_32APSK = DataLabel_32APSK(idxSC,:).';
DataSym_32APSK(1,:,:) = sqrt(PowerVar/2)*dvbsapskmod(DataLabel_32APSK,M_32APSK,'s2'); % Apply APSK modulation to the data.

% Transmitted frame
TransmittedPacket_QPSK = [PilotSym;DataSym_QPSK];
TransmittedPacket_8PSK = [PilotSym;DataSym_8PSK];
TransmittedPacket_16APSK = [PilotSym;DataSym_16APSK];
TransmittedPacket_32APSK = [PilotSym;DataSym_32APSK];

%% Get training feature 

% Mode = S, only CSI signal
% Mode = SE, training model with CSI and Elevation angle
% Mode = SL, CSI signal, local max and local min for elevation angle
% replacement
Mode = 'S';

% Scenario = 1, Suburban and rural scenario
% Scenario = 2, Urban scenario
% Scenario = 3, Dense urban scenario
Scenario = 1;

CSV = [1];
NumCSV = length(CSV);

SatCH = [];
for snr = length(NoiseVar)
    for n = CSV
        disp(snr);
        disp(Name{n});
        
        NoiseVar(n) = NoiseVar(snr);
        
        % Received frame
        ReceivedPacket_QPSK = getLEOChannel(Scenario,TransmittedPacket_QPSK,LengthCP,h,NoiseVar(n),n);
        ReceivedPacket_8PSK = getLEOChannel(Scenario,TransmittedPacket_8PSK,LengthCP,h,NoiseVar(n),n);
        ReceivedPacket_16APSK = getLEOChannel(Scenario,TransmittedPacket_16APSK,LengthCP,h,NoiseVar(n),n);
        ReceivedPacket_32APSK = getLEOChannel(Scenario,TransmittedPacket_32APSK,LengthCP,h,NoiseVar(n),n);
        
        % LS Channel Estimation
        wrapper = @(y,x) lsChanEstimation(y,x,NumPilot,NumSC,idxSC);
        
        % QPSK
        ReceivedPilot = mat2cell(ReceivedPacket_QPSK(1,:,:),1,NumSC,ones(1,NumPacket)); %  get y_pilot
        PilotSeq_nonMask = mat2cell(FixedPilotAll_nonMask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_QPSK_nonMask = cellfun(wrapper,ReceivedPilot,PilotSeq_nonMask,'UniformOutput',false);
        EstChanLS_QPSK_nonMask = cell2mat(squeeze(EstChanLS_QPSK_nonMask));
        SatCH = [SatCH EstChanLS_QPSK_nonMask];
        ReceivedDataSymbol_QPSK = ReceivedPacket_QPSK(2,idxSC,1:end); %  get y_data
        EstSym_QPSK_nonMask = squeeze(ReceivedDataSymbol_QPSK)./EstChanLS_QPSK_nonMask; % y_data/H
        DecLabel_QPSK_nonMask = pskdemod(EstSym_QPSK_nonMask/sqrt(PowerVar/2),M_QPSK,pi/M_QPSK);
        [errSym_QPSK(snr), SER_QPSK(snr)] = symerr(DataLabel_idxSc_QPSK, DecLabel_QPSK_nonMask);      
        fprintf("SER QPSK: %d/%d(%f%%)\n", errSym_QPSK(snr), length(EstSym_QPSK_nonMask), SER_QPSK(snr));
        
        PlotConstellation(DataSym_QPSK, EstSym_QPSK_nonMask, PowerVar, idxSC);
        
        PilotSeq_Mask = mat2cell(FixedPilotAll_Mask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_QPSK_Mask = cellfun(wrapper,ReceivedPilot,PilotSeq_Mask,'UniformOutput',false);
        EstChanLS_QPSK_Mask = cell2mat(squeeze(EstChanLS_QPSK_Mask));
        EstSym_QPSK_Mask = squeeze(ReceivedDataSymbol_QPSK)./EstChanLS_QPSK_Mask; % y_data/H_mask
        DecLabel_QPSK = pskdemod(EstSym_QPSK_Mask/sqrt(PowerVar/2),M_QPSK,pi/M_QPSK);
        [errSym_QPSK_Eve(snr), SER_QPSK_Eve(snr)] = symerr(DataLabel_idxSc_QPSK, DecLabel_QPSK);
        fprintf("SER QPSK Eve: %d/%d(%f%%)\n", errSym_QPSK_Eve(snr), length(EstSym_QPSK_Mask), SER_QPSK_Eve(snr));
        
        PlotConstellation(DataSym_QPSK, EstSym_QPSK_Mask, PowerVar, idxSC);

        
        
        % 8PSK
        ReceivedPilot = mat2cell(ReceivedPacket_8PSK(1,:,:),1,NumSC,ones(1,NumPacket));
        PilotSeq_nonMask = mat2cell(FixedPilotAll_nonMask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_8PSK_nonMask = cellfun(wrapper,ReceivedPilot,PilotSeq_nonMask,'UniformOutput',false);
        EstChanLS_8PSK_nonMask = cell2mat(squeeze(EstChanLS_8PSK_nonMask));
        ReceivedDataSymbol_8PSK = ReceivedPacket_8PSK(2,idxSC,1:end);
        EstSym_8PSK_nonMask = squeeze(ReceivedDataSymbol_8PSK)./EstChanLS_8PSK_nonMask; % y_data/H
        %EstSym_8PSK_nonMask = squeeze(ReceivedDataSymbol_8PSK)./EstChanLS_nonMask; % y_data/H
        DecLabel_8PSK = pskdemod(EstSym_8PSK_nonMask/sqrt(PowerVar/2),M_8PSK,pi/M_8PSK);
        [errSym_8PSK(snr), SER_8PSK(snr)] = symerr(DataLabel_idxSc_8PSK, DecLabel_8PSK);
        fprintf("SER 8PSK: %d/%d(%f%%)\n", errSym_8PSK(snr), length(EstSym_8PSK_nonMask), SER_8PSK(snr));
        
        PlotConstellation(DataSym_8PSK, EstSym_8PSK_nonMask, PowerVar, idxSC);
        
        PilotSeq_Mask = mat2cell(FixedPilotAll_Mask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_8PSK_Mask = cellfun(wrapper,ReceivedPilot,PilotSeq_Mask,'UniformOutput',false);
        EstChanLS_8PSK_Mask = cell2mat(squeeze(EstChanLS_8PSK_Mask));
        EstSym_8PSK_Mask = squeeze(ReceivedDataSymbol_8PSK)./EstChanLS_8PSK_Mask; % y_data/H_mask
        %EstSym_8PSK_Mask = squeeze(ReceivedDataSymbol_8PSK)./EstChanLS_QPSK_Mask; % y_data/H_mask
        DecLabel_8PSK = pskdemod(EstSym_8PSK_Mask/sqrt(PowerVar/2),M_8PSK,pi/M_8PSK);
        [errSym_8PSK_Eve(snr), SER_8PSK_Eve(snr)] = symerr(DataLabel_idxSc_8PSK, DecLabel_8PSK);
        fprintf("SER 8PSK Eve: %d/%d(%f%%)\n", errSym_8PSK_Eve(snr), length(EstSym_8PSK_Mask), SER_8PSK_Eve(snr));
        
        PlotConstellation(DataSym_8PSK, EstSym_8PSK_Mask, PowerVar, idxSC);

        % 16APSK
        ReceivedPilot = mat2cell(ReceivedPacket_16APSK(1,:,:),1,NumSC,ones(1,NumPacket));
        PilotSeq_nonMask = mat2cell(FixedPilotAll_nonMask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_16APSK_nonMask = cellfun(wrapper,ReceivedPilot,PilotSeq_nonMask,'UniformOutput',false);
        EstChanLS_16APSK_nonMask = cell2mat(squeeze(EstChanLS_16APSK_nonMask));
        ReceivedDataSymbol_16APSK = ReceivedPacket_16APSK(2,idxSC,1:end);
        EstSym_16APSK_nonMask = squeeze(ReceivedDataSymbol_16APSK)./EstChanLS_16APSK_nonMask; % y_data/H
        DecLabel_16APSK = dvbsapskdemod(EstSym_16APSK_nonMask/sqrt(PowerVar/2),16,'s2');
        [errSym_16APSK(snr), SER_16APSK(snr)] = symerr(DataLabel_idxSc_16APSK, DecLabel_16APSK);
        fprintf("SER 16APSK: %d/%d(%f%%)\n", errSym_16APSK(snr), length(EstSym_16APSK_nonMask), SER_16APSK(snr));
        %plotCSIEAngle(EstChanLS_16APSK,'CSI',['r','b'],Eb_N0_dB(n),n); 
        
        PlotConstellation(DataSym_16APSK, EstSym_16APSK_nonMask, PowerVar, idxSC);
        
        PilotSeq_16APSK_Mask = mat2cell(FixedPilotAll_Mask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_16APSK_Mask = cellfun(wrapper,ReceivedPilot,PilotSeq_16APSK_Mask,'UniformOutput',false);
        EstChanLS_16APSK_Mask = cell2mat(squeeze(EstChanLS_16APSK_Mask));
        EstSym_16APSK_Mask = squeeze(ReceivedDataSymbol_16APSK)./EstChanLS_16APSK_Mask; % y_data/H_mask
        DecLabel_16APSK = dvbsapskdemod(EstSym_16APSK_Mask/sqrt(PowerVar/2),16,'s2');
        [errSym_16APSK_Eve(snr), SER_16APSK_Eve(snr)] = symerr(DataLabel_idxSc_16APSK, DecLabel_16APSK);
        fprintf("SER 16APSK Eve: %d/%d(%f%%)\n", errSym_16APSK_Eve(snr), length(EstSym_16APSK_Mask), SER_16APSK_Eve(snr));
        
        PlotConstellation(DataSym_16APSK, EstSym_16APSK_Mask, PowerVar, idxSC);
        
        % 32APSK
        ReceivedPilot = mat2cell(ReceivedPacket_32APSK(1,:,:),1,NumSC,ones(1,NumPacket));
        PilotSeq_nonMask = mat2cell(FixedPilotAll_nonMask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_32APSK_nonMask = cellfun(wrapper,ReceivedPilot,PilotSeq_nonMask,'UniformOutput',false);
        EstChanLS_32APSK_nonMask = cell2mat(squeeze(EstChanLS_32APSK_nonMask));
        ReceivedDataSymbol_32APSK = ReceivedPacket_32APSK(2,idxSC,1:end);
        EstSym_32APSK_nonMask = squeeze(ReceivedDataSymbol_32APSK)./EstChanLS_32APSK_nonMask; % y_data/H
        DecLabel_32APSK = dvbsapskdemod(EstSym_32APSK_nonMask/sqrt(PowerVar/2),32,'s2');
        [errSym_32APSK(snr), SER_32APSK(snr)] = symerr(DataLabel_idxSc_32APSK, DecLabel_32APSK);
        fprintf("SER 32APSK: %d/%d(%f%%)\n", errSym_32APSK(snr), length(EstSym_32APSK_nonMask), SER_32APSK(snr));
        
        PlotConstellation(DataSym_32APSK, EstSym_32APSK_nonMask, PowerVar, idxSC);
        
        PilotSeq_Mask = mat2cell(FixedPilotAll_Mask,1,NumPilot,ones(1,NumPacket));
        EstChanLS_32APSK_Mask = cellfun(wrapper,ReceivedPilot,PilotSeq_Mask,'UniformOutput',false);
        EstChanLS_32APSK_Mask = cell2mat(squeeze(EstChanLS_32APSK_Mask));
        EstSym_32APSK_Mask = squeeze(ReceivedDataSymbol_32APSK)./EstChanLS_32APSK_Mask; % y_data/H_mask
        DecLabel_32APSK = dvbsapskdemod(EstSym_32APSK_Mask/sqrt(PowerVar/2),32,'s2');
        [errSym_32APSK_Eve(snr), SER_32APSK_Eve(snr)] = symerr(DataLabel_idxSc_32APSK, DecLabel_32APSK);
        fprintf("SER 32APSK Eve: %d/%d(%f%%)\n", errSym_32APSK_Eve(snr), length(EstSym_32APSK_Mask), SER_32APSK_Eve(snr));
       
        PlotConstellation(DataSym_32APSK, EstSym_32APSK_Mask, PowerVar, idxSC);
        
        % Plot CSI ground truth
        %plotCSIEAngle(EstChanLS,'CSI',['r','b'],Eb_N0_dB(n),n);
        %plotCSIEAngle(EstChanLSMask,'CSI',['r','b'],Eb_N0_dB(n),n);
        
        %plotCSIEAngle(EstChanLS_QPSK_nonMask.*complex(cosd(theta),sind(theta)),'CSI',['r','b'],Eb_N0_dB(n),n);
        %plotCSIEAngle(EstChanLS_QPSK_nonMask,'CSI',['r','b'],Eb_N0_dB(n),n);
        %plotCSIEAngle(EstChanLS_QPSK_Mask,'CSI',['r','b'],Eb_N0_dB(n),n);


        % plotCSI(EstChanLS,'CSI',n,['m','c'],Eb_N0_dB(n));
        %plotCSIFluctuation(EstChanLS,'CSI',['r','b'],Eb_N0_dB(n),n,1);
        %plotCSIFluctuation(EstChanLS,'CSI',['r','b'],Eb_N0_dB(n),n,2);

        
    end
end

%save('SatCH.mat','SatCH');

%% Plot SER

PlotSER(Eb_N0_dB, SER_QPSK, SER_QPSK_Eve, 'QPSK');
PlotSER(Eb_N0_dB, SER_8PSK, SER_8PSK_Eve, '8PSK');
PlotSER(Eb_N0_dB, SER_16APSK, SER_16APSK_Eve, '16APSK');
PlotSER(Eb_N0_dB, SER_32APSK, SER_32APSK_Eve, '32APSK');

function PlotConstellation(referenceSym, receivedSym, PowerVar, idxSC)
    figure();
    plot(receivedSym/sqrt(PowerVar/2),'b*'); hold on;
    referenceSym = squeeze(referenceSym);
    plot(referenceSym(idxSC,:)/sqrt(PowerVar/2),'r+','LineWidth',3,'MarkerSize',10);
    xlabel('In-Phase');
    ylabel('Quadrature');
    xlim([-2,2]);
    ylim([-2,2]);
    legend('Received constellation','Reference constellation');
end

function PlotSER(Eb_N0_dB, SER, SER_Eve, title_name)
    figure();
    semilogy(Eb_N0_dB,SER_Eve,'r-o','LineWidth',3,'MarkerSize',12);hold on;
    semilogy(Eb_N0_dB,SER,'b-o','LineWidth',3,'MarkerSize',12);hold on;
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    xlabel('Eb/N0 (dB)', 'FontSize', 18, 'fontweight','bold','Interpreter','tex');
    ylabel('Symbol Error Rate', 'FontSize', 18, 'fontweight','bold','Interpreter','tex');
    title(title_name);
    saveas(gcf,strcat(title_name,'.png'));
end