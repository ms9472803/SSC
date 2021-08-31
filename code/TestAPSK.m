clear all;
close all;


M = [4 8 20]; %  number of constellation points per PSK ring
radii = [0.25 0.5 1]; %  radius of each PSK ring
modOrder = sum(M);

symbolNum = 45000; %一共有幾個符號
x = randi([0 modOrder-1],symbolNum,1); %隨機產生symbolNum個 [0~ (modOrder-1)]的symbol
txSig = apskmod(x,M,radii); % Apply APSK modulation to the data.

snr = 30; % dB
tmp = txSig;
%txSig = txSig * complex(cosd(30),sind(30));
rxSig = awgn(txSig,snr,'measured'); % Pass the modulated signal through a noisy channel.

z = apskdemod(rxSig,M,radii); % Demodulate the received signal
z2 = apskdemod(rxSig,M,radii,'OutputType','approxllr');
[errSym, SER] = symerr(x, z);

% Plot the transmitted (reference) signal points and the noisy received signal points.
plot(rxSig,'b*'); hold on;
plot(tmp,'r+');
grid
xlim([-3 3]);
ylim([-3 3]);
xlabel('In-Phase');
ylabel('Quadrature');
legend('Received constellation','Reference constellation');

EsN0 = 1:2:25; %符號能量噪聲比, db
txSig = txSig.';
x = x.';
for i=1:length(EsN0)

    rxSig = awgn(txSig,EsN0(i),'measured');
    z = apskdemod(rxSig,M,radii);
    [errSym(i), SER(i)] = symerr(x,z);
end

figure();
semilogy(EsN0,SER,'-o','LineWidth',3,'MarkerSize',12);
set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
xlabel('Es/N0','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
ylabel('SER','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
legend('SER','FontSize', 12);
saveas(gcf,"16APSK.png");