
%% Clear workspace
clear all;
close all;


r_m = [0.7071+0.7071i, 0.7071-0.7071i, -0.7071+0.7071i, -0.7071-0.7071i];
theta_rm = [45, 45+90,  45+180, 45+270];
r_m2 = complex(cosd(theta_rm),sind(theta_rm));

carrier = nrCarrierConfig;
carrier.NSizeGrid = 25;
carrier.SubcarrierSpacing = 15;
carrier.NSlot = 1;
carrier.NFrame = 0

%% Create a CSI-RS configuration object

csirs = nrCSIRSConfig;
csirs.CSIRSType = {'nzp','zp'};
csirs.CSIRSPeriod = {[5 1],[10 1]};
csirs.Density = {'dot5even','one'};
csirs.RowNumber = [3 5];
csirs.SymbolLocations = {1,6};
csirs.SubcarrierLocations = {6,4};
csirs.NumRB = 25;

powerCSIRS = 0;
disp(['CSI-RS power scaling: ' num2str(powerCSIRS) ' dB']);

sym = nrCSIRS(carrier,csirs);
csirsSym = sym*db2mag(powerCSIRS);

csirsInd = nrCSIRSIndices(carrier,csirs);

ports = max(csirs.NumCSIRSPorts);   % Number of antenna ports
txGrid = nrResourceGrid(carrier,ports);

txGrid(csirsInd) = csirsSym;

plotGrid(size(txGrid),csirsInd,csirsSym);


[txWaveform,ofdmInfo] = nrOFDMModulate(carrier,txGrid);

R = 4;

channel = nrTDLChannel;
channel.NumTransmitAntennas = ports;
channel.NumReceiveAntennas = R;
channel.DelayProfile = 'TDL-C';
channel.MaximumDopplerShift = 10;
channel.DelaySpread = 1e-8

chInfo = info(channel);
maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
%txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];

%{
prbs = [];
PR.x2_init = 64;
x1_init = zeros(1,31);
x1_init(1) = 1;
N_slot_symb = 14;
l = 1;
n_ID = 1;
c_init = mod(2.^10*(N_slot_symb*+l+1)*(2*n_ID+1)+n_ID, 2.^31);
for n = 0:119

    % This is supposed to be used as c_init as described in 36.211 7.2
    % But for simplicity, I just se the arbitrary value for now

    x2_init = de2bi(c_init,31,'left-msb'); %轉二進位從右到左 長度31


    % Mpn is the length of the final sequence c()

    % '1' is added because of Matlab Array Index starts with '1'
    Mpn = n + 1;

    % Nc as defined in 36.211 7.2
    Nc = 1600;

    % Create a vector(array) for x1() and x2() all initialized with 0
    x1 = zeros(1,Nc + Mpn + 31);
    x2 = zeros(1,Nc + Mpn + 31);

    % Create a vector(array) for c() all initialized with 0
    c = zeros(1,Mpn);

    % Initialize x1() and x2()
    x1(1:31) = x1_init;
    x2(1:31) = x2_init;
 
    
    for n = 1 : (Mpn+Nc)    
       x1(n+31) = mod(x1(n+3) + x1(n),2);
    end

    
    for n = 1 : (Mpn+Nc)
       x2(n+31) = mod(x2(n+3) + x2(n+2) + x2(n+1) + x2(n),2);
    end

    n = Mpn;
    c(n)= mod(x1(n+Nc) + x2(n+Nc),2);

    prbs = [prbs c(n)];
end
%}
%stem(prbs);xlim([1 length(prbs)]); ylabel('c(n)');

% Mpn is the length of the final sequence c()
Mpn = 150;
c = zeros(1,Mpn);
Nc = 1600;
N_slot_symb = 14;
l = 1;
n_ID = 1;

% Initialize x1() and x2()
x1 = zeros(1,Nc + Mpn + 31);
x2 = zeros(1,Nc + Mpn + 31);
x1_init = zeros(1,31);
x1_init(1) = 1;
x1(1:31) = x1_init;
dec_c_init = mod(2.^10*(N_slot_symb*+l+1)*(2*n_ID+1)+n_ID, 2.^31);
c_init = de2bi(dec_c_init,31,'left-msb'); %轉二進位從右到左 長度31
x2(1:31) = c_init;

% generate the m-sequence : x1()
for n = 1:(Mpn+Nc)
    x1(n+31) = mod(x1(n+3) + x1(n),2);
end

% generate the m-sequence : x2()
for n = 1:(Mpn+Nc)
    x2(n+31) = mod(x2(n+3) + x2(n+2) + x2(n+1) + x2(n),2);
end

c = zeros(1,Mpn);
for n = 1:Mpn
    c(n)= mod(x1(n+Nc) + x2(n+Nc),2);
end

for m = 0:floor((Mpn-1)/2)
    r(m+1) = 1/sqrt(2)*(1-2*c(2*m+1))+j*(1/sqrt(2))*(1-2*c(2*m+1+1));
end

function plotGrid(gridSize,csirsInd,csirsSym)
%    plotGrid(GRIDSIZE,CSIRSIND,CSIRSSYM) plots the carrier grid of size GRIDSIZE
%    by populating the grid with CSI-RS symbols using CSIRSIND and CSIRSSYM.

    figure()
    cmap = colormap(gcf);
    chpval = {20,2};
    chpscale = 0.25*length(cmap); % Scaling factor

    tempSym = csirsSym;
    tempSym(tempSym ~= 0) = chpval{1}; % Replacing non-zero-power symbols
    tempSym(tempSym == 0) = chpval{2}; % Replacing zero-power symbols
    tempGrid = complex(zeros(gridSize));
    tempGrid(csirsInd) = tempSym;

    image(chpscale*tempGrid(:,:,1)); % Multiplied with scaling factor for better visualization
    axis xy;
    names = {'NZP CSI-RS','ZP CSI-RS'};
    clevels = chpscale*[chpval{:}];
    N = length(clevels);
    L = line(ones(N),ones(N),'LineWidth',8); % Generate lines
    % Index the color map and associate the selected colors with the lines
    set(L,{'color'},mat2cell(cmap( min(1+clevels,length(cmap) ),:),ones(1,N),3)); % Set the colors according to cmap
    % Create legend 
    legend(names{:});

    title('Carrier Grid Containing CSI-RS')
    xlabel('OFDM Symbols');
    ylabel('Subcarriers');
end
