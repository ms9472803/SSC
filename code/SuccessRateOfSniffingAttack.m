clear all;
close all;
a = complex(0.2245 ,0.2245);
b = complex(-0.2245, 0.2245);

c = complex(0.2588 ,0.9659);
d = complex(-0.2588 ,0.9659);

theta1 = 90;
m1 = complex(cosd(theta1),sind(theta1)); 

theta2 = 30;
m2 = complex(cosd(theta2),sind(theta2)); 

m3 = complex(cosd(theta1-theta2),sind(theta1-theta2)); 

disp(a*m1);
disp(c*m2);


title_name = '8PSK';
NumIter = 50000;
cnt1 = 0;
cnt2 = 0;
w = 8;
theta_max = 360/w/2;
for i = 1:NumIter
    theta_i = randi([1, 360]);
    theta_j = randi([1, 360]);
    if abs(theta_i - theta_j) < theta_max
        cnt1 = cnt1+1;
    end
end
P1 = cnt1/NumIter;

candidate = [theta_max:theta_max:theta_max*2*w];
for i = 1:NumIter
    theta_i = candidate(randi([1,2*w]));
    theta_j = candidate(randi([1,2*w]));
    %if abs(theta_i - theta_j) < theta_max
    if theta_i == theta_j
        cnt2 = cnt2+1;
    end
end
P2 = cnt2/NumIter;

x = [P1, P2];
figure();
b = bar(x);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(roundn(b(1).YData,-3));
text(xtips1,ytips1,labels1,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',16);
ylabel('Success rate of sniffing attack');
ylim([0,0.4]);
set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
set(gca, 'xticklabel', {'SCF','SSC'});
title(title_name, 'FontSize', 20, 'fontweight','bold','Interpreter','tex');
title_name = strcat('Bar_',title_name);
saveas(gcf,strcat(title_name,'.png'));


theta = [50 100 145];
alpha = [complex(cosd(theta(1)),sind(theta(1))) complex(cosd(theta(2)),sind(theta(2))) complex(cosd(theta(3)),sind(theta(3)))];
alpha_SAT = complex(cosd(45),sind(45));
x_2 = 1.5;
H = complex(randn(3,3), randn(3,3));
F = zeros(3,3);
F(1,:) = H(1,:)*alpha(1)*alpha_SAT;  % 1-row of F
F(2,:) = H(2,:)*alpha(2)*alpha_SAT;  % 2-row of F
F(3,:) = -1*x_2*F(2,:) + H(3,:)*alpha(3)*alpha_SAT; % 3-row of F
F = F + (F/alpha_SAT);
T = H * pinv(F) % pinv: pseudo inverse
1/ (alpha(1)*alpha_SAT^2)
x_2/(alpha(3)*alpha_SAT^2)



figure();
x = [0.232 0.126; 0.124 0.063 ; 0.229 0.124 ; 0.236 0.124];
b = bar(x)
b(1).FaceColor = [0 0 1];
b(2).FaceColor = [155/255 187/255 89/255];
%xtips1 = b(1).XEndPoints;
%ytips1 = b(1).YEndPoints;
%labels1 = string(roundn(b(1).YData,-3));
%text(xtips1,ytips1,labels1,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',16);

legend('SCF','SSC');
ylabel('Success rate of sniffing attack');
ylim([0,0.35]);
set(gca,'FontSize',16,'fontweight','bold','linewidth',2);
set(gca, 'xticklabel', {'QPSK','8PSK','16APSK','32APSK'});

saveas(gcf,'Bar.png');