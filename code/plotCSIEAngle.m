function plotCSIEAngle(Y,Title,color,snr,n)
    load('EAngle.mat');
    EAngle = cell2mat(EAngle(n));
    
    [M1, I1] = min(abs(real(Y))); % find the point nearest to zero
    [M2, I2] = min(abs(imag(Y))); % find the point nearest to zero
    [M3, I3] = max(EAngle);
    [M4, I4] = min(imag(Y)); % find lowest point
    [M5, I5] = min(real(Y)); % find lowest point
    [M6, I6] = max(imag(Y)); % find lowest point
    [M7, I7] = max(real(Y)); % find lowest point
    
    shift = 50;
    I2 = I5;
    I1 = I7;
    
    I1 = 24200;
    I2 = I1;
    
    figure();
    %figure('Name',strcat('CSI Diagram LEO_Track_',num2str(n),'_SNR_',num2str(snr)));
    %plot(I1-shift:I1+shift,real(Y(I1-shift:I1+shift)),'color',color(1),'LineWidth',3); hold on;
    %plot(I2-shift:I2+shift,imag(Y(I2-shift:I2+shift)),'color',color(2),'LineWidth',3); hold on;
    %plot(I3-shift:I3+shift,imag(Y(I3-shift:I3+shift)),'color',color(2),'LineWidth',3); hold on;
    plot(1:size(real(Y)),real(Y),'color',color(1),'LineWidth',3); hold on;
    plot(1:size(imag(Y)),imag(Y),'color',color(2),'LineWidth',3); hold on;
%     plot(1:size(real(Y)),real(Y),'color',[0.8500 0.3250 0.0980],'LineWidth',3); hold on;
%     plot(1:size(imag(Y)),imag(Y),'color',[0 0.4470 0.7410],'LineWidth',3); hold on;
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    set(gca,'xtick',1:floor(size(Y)/10):size(Y),'xticklabel',round(EAngle(1:floor(size(Y)/10):size(Y))));
   %set(gca,'xtick',I1-shift:floor(shift*2/7):I1+shift,'xticklabel',roundn(EAngle(I1-shift:floor(shift*2/7):I1+shift),-2));
    %set(gca,'xtick',I2-shift:floor(shift*2/7):I2+shift,'xticklabel',roundn(EAngle(I2-shift:floor(shift*2/7):I2+shift),-2));
    
    %{
    aaa = EAngle(min(I1,I2)-shift);
    bbb = EAngle(max(I1,I2)+shift);
    ccc = min(I1,I2)-shift;
    ddd = max(I1,I2)+shift;
    set(gca,'xtick',ccc:floor((ddd-ccc)/10):ddd,'xticklabel',round(aaa):floor((bbb-aaa)/10):round(bbb) );
    %}
%    title(strcat('Subcarrier',{' '},num2str(Title),{' '},'of LEO Track',{' '},num2str(n)),'Interpreter','latex');    
%    title(strcat(Title,{' '},'of LEO Track from E to W'),'Interpreter','latex');
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
%}
   %{
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
    xlabel('Elevation Angle (deg)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    ylabel('CSI Value (h)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
%    legend('Real (h)','Imag (h)','Interpreter','tex','Location','southwest','FontSize', 12);
    legend('Real $(h)$','Imag $(h)$','Interpreter','latex','Location','southwest');
%     legend('Real (h)','Interpreter','tex','Location','southwest', 'FontSize', 12);
%     legend('Imag (h)','Interpreter','tex','Location','southwest', 'FontSize', 12);
    legend('boxoff');
    saveas(gcf,strcat(string(n),".png"));
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
%{
    left = 1;
    right = 1;
    index = 1;
    index = ceil(length(EAngle)/2);
    for i = EAngle
        %disp(i);
        if i > 10 && left == 1
           left = index;
        end
        
        if i > 13 && right == 1
           right = index; 
        end
        index = index + 1;
    end
    disp(left);
    disp(right);
    left = 1;
    right =3000;
    figure();
    plot(EAngle(left:right),real(Y(left:right)),'color',color(1),'LineWidth',3); hold on;
    plot(EAngle(left:right),imag(Y(left:right)),'color',color(2),'LineWidth',3); hold on;
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    xlabel('Elevation Angle (deg)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    ylabel('CSI Value (h)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    legend('Real (h)','Imag (h)','Interpreter','tex','Location','southwest','FontSize', 12);
    legend('boxoff');
    saveas(gcf,"4.png");
%}
%{
    figure();
    rev_real = real(Y) * -1;
    rev_imag = imag(Y) * -1;
    plot(1:size(real(Y)),rev_real,'color',color(1),'LineWidth',3); hold on;
    plot(1:size(imag(Y)),rev_imag,'color',color(2),'LineWidth',3); hold on;
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    set(gca,'xtick',1:floor(size(Y)/10):size(Y),'xticklabel',round(EAngle(1:floor(size(Y)/10):size(Y))));
    xlabel('Elevation Angle (deg)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    ylabel('CSI Value (h)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    legend('Real $(h)$','Imag $(h)$','Interpreter','latex','Location','southwest');
    legend('boxoff');
%}

end