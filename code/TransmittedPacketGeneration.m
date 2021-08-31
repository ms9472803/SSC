function TransmittedPacket = TransmittedPacketGeneration(MOD,NumSC,NumPacket,idxSC)
% This function is to Generate Transmitted frame for different modulation
% include pilot symbol and data symbol


NumPilot = 64;
% Pilot symbols - Fixed during the whole transmission
FixedPilot = sqrt(PowerVar/2)*complex(sign(rand(1,NumPilot)-0.5),sign(rand(1,NumPilot)-0.5)); 

CSIRS = [0.7071+0.7071i, -0.7071+0.7071i, -0.7071-0.7071i, 0.7071-0.7071i];
theta_rm = [45, 45+90,  45+180, 45+270];
knownSeq1 = zeros(1,NumPacket);

for j = 1:length(knownSeq1)
    knownSeq1(j) = CSIRS(randi([1,4],1,1));
end

knownSeq2 = zeros(1,NumPacket);
theta = 60;
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

knownSeq1 = sqrt(PowerVar/2) * knownSeq1;
knownSeq2 = sqrt(PowerVar/2) * knownSeq2;

FixedPilotAll_nonMask = repmat(FixedPilot,1,1,NumPacket); 
FixedPilotAll_Mask = repmat(FixedPilot,1,1,NumPacket); 

for i = 1:length(knownSeq2)
    FixedPilotAll_nonMask(:,:,i) = knownSeq1(i);
    FixedPilotAll_Mask(:,:,i) = knownSeq2(i);
end

PilotSym = FixedPilotAll_nonMask;

if MOD == 'QPSK'
    disp('QPSK');
    
    TransmittedPacket = [PilotSym;DataSym];
elseif MOD == '8PSK'
    disp('8PSK');
    
    TransmittedPacket = [PilotSym;DataSym];
elseif MOD == '16APSK'
    disp('16APSK');
    M = [4 12];
    radii = [0.5 1];
    modOrder = sum(M);
    DataLabel = randi([0 modOrder-1],NumSC,NumPacket); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
    DataLabel_idxSc = DataLabel(idxSC,:).';
    DataSym(1,:,:) = sqrt(PowerVar/2)*apskmod(DataLabel,M,radii); % Apply APSK modulation to the data.
    
    TransmittedPacket = [PilotSym;DataSym];
    
elseif MOD == '32APSK'
    disp('32APSK');
    M = [4 8 20];
    radii = [0.25 0.5 1];
    modOrder = sum(M);
    DataLabel = randi([0 modOrder-1],NumSC,NumPacket); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
    DataLabel_idxSc = DataLabel(idxSC,:).';
    DataSym(1,:,:) = sqrt(PowerVar/2)*apskmod(DataLabel,M,radii); % Apply APSK modulation to the data.
    
    TransmittedPacket = [PilotSym;DataSym];
else
    disp('err');
end