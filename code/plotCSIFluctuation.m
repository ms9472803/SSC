function plotCSIFluctuation(Y,Title,color,snr,n,ch)
    load('EAngle.mat');
    EAngle = cell2mat(EAngle(n));
    %{
    left = 1;
    right = 1;
    index = 1;
    index = ceil(length(EAngle)/2);
    for i = EAngle(ceil(length(EAngle)/2):end)
        %disp(i);
        if i < 30.5 && left == 1
           left = index;
        end
        
        if i < 30.1 && right == 1
           right = index; 
        end
        index = index + 1;
    end
    right = right;
    left = left;
    disp(left);
    disp(right);
    newY = Y(left:right);
    
    figure('Name',strcat('CSI Diagram LEO_Track_',num2str(n),'_SNR_',num2str(snr)));
    if ch == 1
        plot(1:size(real(newY)),real(newY),'color',color(1),'LineWidth',1); hold on;
    end
    if ch == 2
        plot(1:size(imag(newY)),imag(newY),'color',color(2),'LineWidth',1); hold on;
    end
    %}
    
     %[M, I] = max(imag(Y)); % find largest point
     %rPoint = [I M];
     %[M, I] = min(real(Y)); % find lowest point
     %rPoint = [rPoint; I M];
     %[M, I] = min(abs(real(Y))); % find the point nearest to zero
     %rPoint = [rPoint; I M];
     [M, I] = min(abs(imag(Y))); % find the point nearest to zero
     %rPoint = [rPoint; I M];
     
    figure();
    tmp = I;
    shift = 500;
    plot(EAngle(tmp-shift:tmp+shift),imag(Y(tmp-shift:tmp+shift)),'color',color(2),'LineWidth',3); hold on;
    %set(gca,'xticklabel',EAngle(tmp-shift:tmp+shift));
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    %set(gca,'ytick',min(newY):0.002:1.718,'yticklabel',min(newY):0.002:1.718);
    
    xlabel('Elevation Angle (deg)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    ylabel('CSI Value (h)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    legend('Imag(h)','Interpreter','tex','Location','southwest', 'FontSize', 12);
    %legend('Real $(h)$','Imag $(h)$','Interpreter','latex','Location','southwest');
    legend('boxoff');
    saveas(gcf,strcat(string(ch),".png"));
    
    
    [M, I] = min(abs(real(Y))); % find the point nearest to zero
     %rPoint = [rPoint; I M];
    
    figure();
    tmp = I;
    disp(tmp);
    plot(EAngle(tmp-shift:tmp+shift),real(Y(tmp-shift:tmp+shift)),'color','r','LineWidth',3); hold on;
    %set(gca,'xticklabel',EAngle(tmp-shift:tmp+shift));
    set(gca,'FontSize',12,'fontweight','bold','linewidth',1.5);
    %set(gca,'ytick',min(newY):0.002:1.718,'yticklabel',min(newY):0.002:1.718);
    
    %xlabel('Elevation Angle (deg)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    ylabel('CSI Value (h)','Interpreter','tex', 'FontSize', 18, 'fontweight','bold');
    legend('Real (h)','Interpreter','tex','Location','southwest', 'FontSize', 12);
    %legend('Real $(h)$','Imag $(h)$','Interpreter','latex','Location','southwest');
    legend('boxoff');
    saveas(gcf,strcat(string(ch),".png"));
    
end