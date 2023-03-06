clear all;
close all;
clc

Fs = 1000;
mi1 = 0.01;
mi2 = mi1;
mi3 = mi1;
mi4 = mi1;


mi = 0.01;

load('Z.mat');
load('D.mat');

z = Z(1,:);
d = D(1,:);

t = 0:1/Fs:(6000 - 1)/Fs;

x_f50 = cos(2*pi*50*t);
x_f100 = 2*x_f50.^2 - 1;
x_f150 = 4*x_f50.^3 - 3*x_f50;
x_f200 = 8*x_f50.^4 - 4*x_f100 - 3;

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

figure(1)
plot(z);
figure(2)
plot(d);
figure(3)
plot(e);
figure(4)
hold all;
plot(W(1,:));
plot(W(2,:));
plot(W(3,:));
plot(W(4,:));
plot(W(5,:));
plot(W(6,:));
plot(W(7,:));
plot(W(8,:));
