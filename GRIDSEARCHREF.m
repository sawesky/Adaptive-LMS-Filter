clear all;
close all;
clc

Fs = 1000;
mi1 = 0.001;
mi2 = [0.001 0.00125 0.0025 0.00375 0.005 0.01];
mi3 = mi2;
mi4 = mi2;

w = zeros(4,216);
J = zeros(1,216);

mi = 0.01;

load('Z.mat');
load('D.mat');

z = Z(1,:);
d = D(1,:);

t = 0:1/Fs:(6000 - 1)/Fs;
rng(1) 
x_f50_sum = cos(2*pi*50*t);
% x_f100 = 2*x_f50.^2 - 1;
% x_f150 = 4*x_f50.^3 - 3*x_f50;
% x_f200 = 8*x_f50.^4 - 4*x_f100 - 3;

for b = 1:6
    for c = 1:6
        for f = 1:6
w1 = zeros(1,length(t)); w1(6) = 0;
w2 = zeros(1,length(t)); w2(6) = 0;
w3 = zeros(1,length(t)); w3(6) = 0;
w4 = zeros(1,length(t)); w4(6) = 0;
x_f50 = zeros(1, length(t));
x_f100 = zeros(1, length(t));
x_f150 = zeros(1, length(t));
x_f200 = zeros(1, length(t));

for i = 6:length(t)
    x_f50(i) = x_f50_sum(i) - w1(i);
    w1(i + 1) = w1(i) + mi1*x_f50(i);
%     xac = x_f50_sum(i:-1:i-4);
%     x_f(i) = 2*20/(pi*(xac*xac' + 0.000001))/xac(1);
    x_f100(i) = 2*x_f50(i)^2 - w2(i);
    w2(i + 1) = w2(i) + mi2(b)*x_f100(i);
    x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(i);
    w3(i + 1) = w3(i) + mi3(c)*x_f150(i);
    x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(i);
    w4(i + 1) = w4(i) + mi4(f)*x_f200(i);
end
w(:, (b-1)*36 + (c-1)*6 + f) = [w1(end); w2(end); w3(end); w4(end)];
J((b-1)*36 + (c-1)*6 + f) = (w1(end)^2 + (1-w2(end))^2 + w3(end)^2 + (3-w4(end))^2)/4;

        end
    end
end

figure, stem(J)

figure
hold on;
plot(w(1,:));
plot(w(2,:));
plot(w(3,:));
plot(w(4,:));
hold off;

% w1 = zeros(1,length(t)); w1(6) = 0;
% w2 = zeros(1,length(t)); w2(6) = 0;
% w3 = zeros(1,length(t)); w3(6) = 0;
% w4 = zeros(1,length(t)); w4(6) = 0;
% x_f50 = zeros(1, length(t));
% x_f100 = zeros(1, length(t));
% x_f150 = zeros(1, length(t));
% x_f200 = zeros(1, length(t));
% 
% for i = 6:length(t)-1
%     x_f50(i) = x_f50_sum(i) - w1(i);
%     w1(i + 1) = w1(i) + mi1*x_f50(i);
% %     xac = x_f50_sum(i:-1:i-4);
% %     x_f(i) = 2*20/(pi*(xac*xac' + 0.000001))/xac(1);
%     x_f100(i) = 2*x_f50(i)^2 - w2(i);
%     w2(i + 1) = w2(i) + mi2*x_f100(i);
%     x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(i);
%     w3(i + 1) = w3(i) + mi3*x_f150(i);
%     x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(i);
%     w4(i + 1) = w4(i) + mi4*x_f200(i);
% end
% 
% figure(4)
% hold on;
% plot(w1);
% plot(w2);
% plot(w3);
% plot(w4);
% hold off;
% 
% disp(w1)
% disp(w2)
% disp(w3)
% disp(w4)