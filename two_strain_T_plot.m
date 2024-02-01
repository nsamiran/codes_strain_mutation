clc 
clear all
close all
format long
set(0,'DefaultAxesFontSize',13);

beta=0.82; gamma=0.78;
mu=0.0001:0.00001:0.004;
a=3:0.01:20;
R0=beta/gamma;
for i=1:length(a)
    for j=1:length(mu)
T(i,j)=(1.53*a(i))/(2*sqrt(mu(j)*gamma*(R0-1)));
    end
end
surf(mu,a,T); shading interp;
hold on
xlabel('$\mu$','interpreter','latex');
ylabel('$a$','interpreter','latex');
zlabel('$\mathcal{T}$','interpreter','latex');
