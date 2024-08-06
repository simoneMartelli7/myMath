clear 
close all 
clc
format short e 

heun = readmatrix("heun.csv");
euler = readmatrix("euler.csv");
n = 101;

hu1 = heun(1:2:end);
hu2 = heun(2:2:end);

eu1 = euler(1:2:end);
eu2 = euler(2:2:end);

t = linspace(0, 10, n);

figure()
plot(t, eu1, t, hu1);
grid on

figure()
plot(t, eu2, t, hu2);
grid on 