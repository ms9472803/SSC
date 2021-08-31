
clear all;
close all;


ofdmMod = comm.OFDMModulator('FFTLength',256, ...
    'NumGuardBandCarriers',[12; 11], ...
    'NumSymbols', 5, ...
    'NumTransmitAntennas', 3, ...
    'PilotInputPort',true, ...
    'Windowing', true, ...
    'WindowLength', 6);

pilotIndOdd = [20; 58; 96; 145; 182; 210];
pilotIndEven = [35; 73; 111; 159; 197; 225];

pilotIndicesAnt1 = cat(2, pilotIndOdd, pilotIndEven, pilotIndOdd, pilotIndEven, pilotIndOdd);
%pilotIndicesAnt1 = cat(2, pilotIndOdd, pilotIndEven, pilotIndOdd, pilotIndEven);

pilotIndicesAnt2 = pilotIndicesAnt1 + 5;
pilotIndicesAnt3 = pilotIndicesAnt1 - 5;

%ofdmMod.PilotCarrierIndices = cat(3, pilotIndicesAnt1, pilotIndicesAnt2, pilotIndicesAnt3);
ofdmMod.PilotCarrierIndices = cat(3, pilotIndicesAnt1, pilotIndicesAnt2);
ofdmMod.PilotCarrierIndices = pilotIndicesAnt1;

ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = 3;

dims = info(ofdmMod)

dataIn = complex(randn(dims.DataInputSize), randn(dims.DataInputSize));
pilotIn = complex(randn(dims.PilotInputSize), randn(dims.PilotInputSize));

modOut = ofdmMod(dataIn,pilotIn);
% 1360 = (227+6+12+11+16)*5 (data+pilot+GB+CP)

%chanGain = complex(randn(3,2), randn(3,2));
chanGain = complex(randn(ofdmMod.NumTransmitAntennas,ofdmDemod.NumReceiveAntennas), randn(ofdmMod.NumTransmitAntennas,ofdmDemod.NumReceiveAntennas));

chanOut = modOut * chanGain;
chanOutAWGN = awgn(chanOut,10);

[dataOut,pilotOut] = ofdmDemod(chanOut);

showResourceMapping(ofdmMod)

pilotCompare = abs(pilotIn(:,:,1)*chanGain(1,1)) - abs(pilotOut(:,:,1,1));
max(pilotCompare(:) < 1e-10)


