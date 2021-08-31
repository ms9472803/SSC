
clear all;
close all;


c = zeros(1, 360);
for theta = 1:360
    c(theta) = complex(cosd(theta),sind(theta));
end

a = 0.7071+0.7071i;
c = c*a;
abs_c = abs(c);

mask = complex(cosd(225),sind(225)) / complex(cosd(45),sind(45));
%mask = complex(cosd(55),sind(55)) * complex(cosd(45),sind(45))
r_mask = real(mask);
i_mask = imag(mask);
if r_mask >= 0 && i_mask >= 0
    cos_mask = mod(acosd(r_mask),360);
    sin_mask = mod(asind(i_mask),360);
    theta = cos_mask;
elseif r_mask >= 0 && i_mask <= 0
    cos_mask = mod(acosd(r_mask),360);
    sin_mask = mod(asind(i_mask),360);
    theta = sin_mask;
elseif r_mask <= 0 && i_mask >= 0
    cos_mask = mod(acosd(r_mask),360);
    sin_mask = mod(asind(i_mask),360);
    theta = cos_mask;
elseif r_mask <= 0 && i_mask <= 0
    cos_mask = mod(acosd(r_mask),360);
    sin_mask = mod(asind(i_mask),360);
    theta = cos_mask + 2*(180-cos_mask);
end
disp(theta);

angle = 10;
angle_3dB = 0.8;
u = 2.07123 * (sind(angle)/sind(angle_3dB));
b_max = 45.6;
BeamGain_angle = b_max*( (besselj(1,u)/2*u) + (besselj(3,u)/u^3) )^2
ChannelFadingCoefficient = 1;
H = ChannelFadingCoefficient*sqrt(BeamGain_angle);

 x=[1:100];
 for i=1:5
 y(:,i)=i*log(x);
 end
colorstring = 'kbgyr';
figure(1); cla;
hold on
for i = 1:5
  plot(x,y(:, i), 'Color', colorstring(i))
end
