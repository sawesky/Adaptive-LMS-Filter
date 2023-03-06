clear all;
close all;
clc

Fs = 1000;
mi1 = 0.001;
mi2 = 0.00125;
mi3 = 0.001;
mi4 = 0.0025;


mi = 0.05;

load('Z.mat');
load('D.mat');

z = Z(1,:);
d = D(1,:);

t = 0:1/Fs:(6000 - 1)/Fs;
rng(1);
x_f50_sum = cos(2*pi*50*t);
% x_f100 = 2*x_f50.^2 - 1;
% x_f150 = 4*x_f50.^3 - 3*x_f50;
% x_f200 = 8*x_f50.^4 - 4*x_f100 - 3;
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
    w2(i + 1) = w2(i) + mi2*x_f100(i);
    x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(i);
    w3(i + 1) = w3(i) + mi3*x_f150(i);
    x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(i);
    w4(i + 1) = w4(i) + mi4*x_f200(i);
end

figure(4)
hold on;
plot(w1);
plot(w2);
plot(w3);
plot(w4);
hold off;

%dodaj
% for i = 6:length(t)
%     x_f50(i) = x_f50_sum(i) - w1(end);
%     x_f100(i) = 2*x_f50(i)^2 - w2(end);
%     x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(end);
%     x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(end);
% end



%%
e = zeros(1, 6000);
W = zeros(8, 6000);

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,i)'*xn;
    W(:,i + 1) = W(:, i) + mi*e(i)*xn;
    
end

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,end-1)'*xn;
    
end

figure(5)
plot(z);
figure(6)
plot(d);
figure(7)
plot(e);
figure(8)
hold all;
plot(W(1,:));
plot(W(2,:));
plot(W(3,:));
plot(W(4,:));
plot(W(5,:));
plot(W(6,:));
plot(W(7,:));
plot(W(8,:));




%%
% idealni slucaj
% x_f50 = cos(2*pi*50*t);
% x_f100 = 2*x_f50.^2 - 1;
% x_f150 = 4*x_f50.^3 - 3*x_f50;
% x_f200 = 8*x_f50.^4 - 4*x_f100 - 3;
% 
% e = zeros(1, 6000);
% W = zeros(8, 6000);
% 
% for i = 6:6000
%   
%     x_q50 =  x_f50(i - round(Fs./(4*50)));
%     x_q100 =  x_f100(i - round(Fs./(4*100)));
%     x_q150 =  x_f150(i - round(Fs./(4*150)));
%     x_q200 =  x_f200(i - round(Fs./(4*200)));
%     xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
%     
%     e(i) = d(i) - W(:,i)'*xn;
%     W(:,i + 1) = W(:, i) + mi*e(i)*xn;
%     
% end
% 
% for i = 6:6000
%   
%     x_q50 =  x_f50(i - round(Fs./(4*50)));
%     x_q100 =  x_f100(i - round(Fs./(4*100)));
%     x_q150 =  x_f150(i - round(Fs./(4*150)));
%     x_q200 =  x_f200(i - round(Fs./(4*200)));
%     xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
%     
%     e(i) = d(i) - W(:,end)'*xn;
%     
% end
% 
% figure(9)
% plot(z);
% figure(10)
% plot(d);
% figure(11)
% plot(e(7:end));
% figure(12)
% hold all;
% plot(W(1,:));
% plot(W(2,:));
% plot(W(3,:));
% plot(W(4,:));
% plot(W(5,:));
% plot(W(6,:));
% plot(W(7,:));
% plot(W(8,:));

