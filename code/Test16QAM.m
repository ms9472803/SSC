clear all;clc;

symbolNum = 45000; %一共有幾個符號
M = 16; %Order 表示16QAM
graycode = [0 1 3 2 4 5 7 6 12 13 15 14 8 9 11 10]; %定義格雷映射, Decimal
EsN0 = 5:20; %符號能量噪聲比, db
snr = 10.^(EsN0/10); %db轉為非線性
msg = randi([0,15],1,symbolNum); %隨機產生symbolNum個0~15的符號
graycode_msg = graycode(msg+1); %格雷映射
msg_mod = qammod(graycode_msg,M); %使用qamod函數, 得到調變後的符號
scatterplot(msg_mod);%画出星座点图
%{
[-3+3i, -3+1i, -3-1i, -3-3i,
-1+3i, -1+1i, -1-1i, -1-3i,
 1+3i,  1+1i,  1-1i,  1-3i,
 3+3i,  3+1i,  3-1i,  3-3i] 
%}
spow = norm(msg_mod).^2/symbolNum; %a+bj取norm的平方,即功率;功率除以符號數得到平均功率
for i=1:length(EsN0)
    sigma = sqrt(spow/(2*snr(i)));
    % 星座點圖乘以隨機長度高斯白噪聲 
    rx = msg_mod+sigma*(randn(1,length(msg_mod))+1i*randn(1,length(msg_mod))); 
    
    y = qamdemod(rx,M); %使用qamdemod函數, 解調rx, 轉為對應的點
    graycode_decmsg = graycode(y+1); %格雷逆映射
    [errBit(i), BER(i)] = biterr(msg, graycode_decmsg, log2(M)); %比較bit誤差
    [errSym(i), SER(i)] = symerr(msg, graycode_decmsg); %比較符號誤差
end

scatterplot(rx);%画出星座点图

figure()
%ber仿真值，ser比特仿真值，ber1理论误比特率，ser1理论误码率
%semilogy(EsN0,BER,'o' ,EsN0,SER,'*',EsN0,ser1,'-',EsN0,ber1,'-.');
semilogy(EsN0,SER,'-o' ,'LineWidth',3,'MarkerSize',12);
set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
%title('16QAM-AWGN')
xlabel('Es/N0','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
ylabel('SER','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
legend('SER','FontSize', 12);
%legend('ber simulation','ser simulation','ser theory' ,'ber theory');
saveas(gcf,"16QAM.png");