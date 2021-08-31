function plotCSIEAngle(Y1,Y2,Title,Trace,color,snr)
    load('EAngle.mat');
    EAngle = cell2mat(EAngle(Trace));
    figure('Name',strcat('CSI Diagram LEO_Track_',num2str(Trace),'_SNR_',num2str(snr)));
    plot(1:size(real(Y1)),real(Y1),'color',color(1),'LineWidth',1); hold on;
    plot(1:size(real(Y2)),real(Y2),'color','r','LineWidth',1); hold on;
    set(gca,'xtick',1:floor(size(Y1)/10):size(Y1),'xticklabel',round(EAngle(1:floor(size(Y1)/10):size(Y1))));
%    plot(1:size(real(Y)),real(Y),'color',[0 0.4470 0.7410],'LineWidth',1); hold on;
%    plot(1:size(imag(Y)),imag(Y),'color',[0.8500 0.3250 0.0980],'LineWidth',1); hold off;

    %title(strcat(Title,{' '},'of LEO Track',{' '},num2str(Trace)),'Interpreter','latex');    
%    title(strcat(Title,{' '},'of LEO Track from SW to NE'),'Interpreter','latex');
%    title(strcat(Title,{' '},'of Cloud and Fog Condition'),'Interpreter','latex');
    
    %{
    [M, I] = max(real(Y)); % find largest point
    rPoint = [I M];
    [M, I] = min(real(Y)); % find lowest point
    rPoint = [rPoint; I M];
    [M, I] = min(abs(real(Y))); % find the point nearest to zero
    rPoint = [rPoint; I M];
    plot(rPoint(:,1),rPoint(:,2),'s','MarkerSize',10,'MarkerFaceColor',color(1));
    text(rPoint(:,1),rPoint(:,2),strcat({'  '},num2str(rPoint(:,1))));
    
    [M, I] = max(imag(Y)); 
    iPoint = [I M];
    [M, I] = min(imag(Y));
    iPoint = [iPoint; I M];
    [M, I] = min(abs(imag(Y))); % find the point nearest to zero
    iPoint = [iPoint; I M];
    plot(iPoint(:,1),iPoint(:,2),'s','MarkerSize',10,'MarkerFaceColor',color(2));
    text(iPoint(:,1),iPoint(:,2),strcat({'  '},num2str(iPoint(:,1))));
    %}

%    legend('Max Real $(h)$','Max Imag $(h)$','Interpreter','latex');
    xlabel('Elevation Angle (deg)','Interpreter','latex');
    ylabel('CSI Value $(h)$','Interpreter','latex');
    legend('Ground Truth Real $(h)$','GRU Real $(h)$','Interpreter','latex','Location','northwest');
%     legend('Real $(h)$','Interpreter','latex','Location','southwest');
    legend('boxoff');
    
%     figure('Name','Peak and Trough');
%     loc = rPoint(3,1)-50:rPoint(3,1)+50;
%     loci = rPoint(3,1)-50:20:rPoint(3,1)+50;
%     plot(loc,real(Y(loc)),'color',color(1),'LineWidth',1); hold on;
%     %plot(loc,imag(Y(loc)),'color',color(2),'LineWidth',1); hold off;
%     set(gca,'xtick',loci,'xticklabel',EAngle(loci));
%     title('Peak and Trough','Interpreter','latex');    
%     xlabel('Elevation Angle (deg)','Interpreter','latex');
%     ylabel('CSI Value $(h)$','Interpreter','latex');
%     legend('Real $(h)$','Interpreter','latex','Location','southwest');
%     legend('boxoff');

    figure('Name',strcat('CSI Diagram LEO_Track_',num2str(Trace),'_SNR_',num2str(snr)));
    plot(1:size(imag(Y1)),imag(Y1),'--','color',color(1),'LineWidth',1); hold on;
    plot(1:size(imag(Y2)),imag(Y2),'--','color','b','LineWidth',1); hold on;
    set(gca,'xtick',1:floor(size(Y1)/10):size(Y1),'xticklabel',round(EAngle(1:floor(size(Y1)/10):size(Y1))));
    xlabel('Elevation Angle (deg)','Interpreter','latex');
    ylabel('CSI Value $(h)$','Interpreter','latex');
    legend('Ground Truth Imag $(h)$','GRU Imag $(h)$','Interpreter','latex','Location','northwest');
    legend('boxoff');
    saveas(gcf,strcat(string(snr),".png"));
end