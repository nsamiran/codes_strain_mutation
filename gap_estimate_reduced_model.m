clc
clear all
close all
set(0,'DefaultAxesFontSize',15);






L=100; %domain size
dx=0.1;
P=1+ L/dx;
eps=0.01;
Mid=(P+1)/2;
% Parameter values
D=0.0002; h1=10;HH1=h1/dx;h2=9;HH2=h2/dx;
tmax=20000;
dt=0.5;

% for j=Mid-20:Mid+20
%    D(j)=DD;
% end
% for j=1:Mid-21
%    D(j)=(1+0.2*rand*(-1)^j)*DD+0.0002*10^-4;%(1+0.001*rand)*DD-0.0002*10^-5;
% end
% for j=Mid+21:P
%    D(j)=(1+0.2*rand*(-1)^j)*DD-0.0002*10^-4;
% end
% for j=Mid-40:Mid+40
%    D(j)=DD;
% end


% for i=1:length(D)
%     D(i)=(1+0.1*rand*(-1)^i)*DD;
% end



% for j=1:Mid-21
%    D(j)=10*DD;%(1+0.001*rand)*DD-0.0002*10^-5;
% end
% for j=Mid-20:Mid
%    D(j)=DD;%(1+0.001*rand)*DD-0.0002*10^-5;
% end
% for j=Mid+1:Mid+20
%    D(j)=0.1*DD;
% end
% for j=Mid+21:P
%    D(j)=10*DD;
% end
% for j=Mid+1:P
%    D(j)=0.1*DD;
% end




for j=1:P
    x(j)=-L/2+(j-1)*dx;
    r(j)=D*dt/(dx)^2;
end
% PPP=floor(P/10);
% for m=1:PPP
% for j=(m-1)*PPP+1:m*PPP
%     r(j)=D(j)*dt/(dx)^2;
% end
% end

% PPP=floor(P/10);
% for m=1:PPP/2
% for j=(m-1)*PPP+1:m*PPP
%     r(j)=(1+.001)*D*dt/(dx)^2;
% end
% end
% for m=1+PPP/2:PPP
% for j=(m-1)*PPP+1:m*PPP
%     r(j)=(1-.001)*D*dt/(dx)^2;
% end
% end

a=0.001;
%  sigma2=0.*x;
%  sigma2(Mid-xx0:Mid-xx0)=(a*(x(1:Mid-xx0)+x0).^2)./(1+(a*(x(1:Mid-xx0)+x0).^2));
% sigma2(Mid+xx0:length(x))=(a*(x(Mid+xx0:length(x))-x0).^2)./(1+(a*(x(Mid+xx0:length(x))-x0).^2));

% sigma1=@(zz) 1-exp(-0.0000001*zz.^2);
% sigma2=@(zz) 1-exp(-0.0000001*zz.^2);
% syms qq
% sigma11=piecewise(qq <= -x0,0,(-x0 < qq) & (qq < x0),1/x0,qq >=x0,0);
% sigma1=@(zz) subs(sigma11,qq,[zz]);
% sigma22=piecewise(qq <= -x00,0,(-x00 < qq) & (qq < x00),1/x00,qq >=x00,0);
% sigma2=@(zz) subs(sigma22,qq,[zz]);

tmax=200000;
dt=0.05;
% r=D*dt/(dx)^2;
MM=tmax/dt;
N1=10^8;N2=5*10^8;
beta11=0.87/N1; beta12=0.000001/N1; beta2=0.82/N2; gamma1=0.85; nu1=0.0; gamma2=0.78; nu2=0.0;
%beta11=1.1/N1; beta12=0.000001/N1; beta2=0.98/N2; gamma1=0.85; nu1=0.01; gamma2=0.78; nu2=0.03;
 %beta11=0.01; beta12=0.001; beta2=0.01; gamma1=0.004; nu1=0.001; gamma2=0.003; nu2=0.001;
%beta11=0.4/N1; beta12=0.001/N1; beta2=0.5/N2; gamma1=0.4; nu1=0.004; gamma2=0.1; nu2=0.02;

%beta11=0.32/N1; beta12=0.001/N1; beta2=0.37/N2; gamma1=0.5; nu1=0.005; gamma2=0.3; nu2=0.04;
%beta11=0.32/N1; beta12=0.001/N1; beta2=0.37/N2; gamma1=0.5; nu1=0.005; gamma2=0.3; nu2=0.04;

%Y(1:2000)=0.1;Z(1:2000)=0.0000002;Y(2001:4000)=0.1+0.001;Z(2001:4000)=0.0000002+0.001;Y(4001:6000)=0.1;Z(4001:6000)=0.0000002;Y(6001:P)=0.1+.001;Z(6001:P)=0.0000002+.001;



% Y=0.1+0.1*cos(8*pi*x/L);Z=0.000011+0.000001*cos(8*pi*x/L);Q= 0.000011 +0.000001*cos(8*pi*x/L);
% Y=1.0000001*Y; Z=1.00000001*Z; Q=1.00000001*Q;
S1=N1-0*(1-x+x); I1=0*(1-sin(x));R1=0*(1-x+x);D1=0*(1-x+x);I1(Mid-1:Mid+1)=0;
S11=S1; I11=I1;R11=R1;D11=D1;
S2=N2-0*rand*(1-x+x); I2=0*(1-sin(x));R2=0*(1-x+x);D2=0*(1-x+x);
for j=Mid-1:Mid+1
I2(j)=20;
end

S22=S2; I22=I2;R22=R2;D22=D2;
W11=[];W12=[];W13=[];W14=[];W21=[];W22=[];W23=[];W24=[];

CC=0;

for n=1:MM

    

      PP1=P+2*HH1;
%    U(1:HH+1)=0;U(HH+2:P+HH+1)=Z(1:P);U(P+HH+2:P+2*HH+1)=0;
%    
%    J(1)=(dx/2)*(U(HH+1-HH)+U(HH+1+HH)+2*sum(U(HH+1-HH+1:HH+1+HH-1)))/(2*h);
%    J(P)=(dx/2)*(U(HH+P+1-HH)+U(HH+P+1+HH)+2*sum(U(HH+P+1-HH+1:HH+P+1+HH-1)))/(2*h);
%    for j=HH+2:HH+P
%      J(j-HH)=(dx/2)*(U(j-HH)+U(j+HH)+2*sum(U(j-HH+1:j+HH-1)))/(2*h);
%    end

  U1(1:HH1)=I1(1); U1(HH1+1:P+HH1)=I1(1:P);U1(P+HH1+1:P+2*HH1)=I1(P);
   
%    J(1)=(dx/2)*(U(HH+1-HH)+U(HH+1+HH)+2*sum(U(HH+1-HH+1:HH+1+HH-1)))/(2*h);
%    J(P)=(dx/2)*(U(HH+P-HH)+U(HH+P+HH)+2*sum(U(HH+P-HH+1:HH+P+HH-1)))/(2*h);
   for j=HH1+1:HH1+P
%      J1(j-HH1)=(dx/2)*(U1(j-HH1)+U1(j+HH1)+2*sum(U1(j-HH1+1:j+HH1-1)))/(2*h1);
J1(j-HH1)=(dx/2)*(U1(j-HH1)+U1(j+HH1)+2*sum(U1(j-HH1+1:j+HH1-1)));
   end

   U2(1:HH2)=I2(1); U2(HH2+1:P+HH2)=I2(1:P);U2(P+HH2+1:P+2*HH2)=I2(P);
   
%    J(1)=(dx/2)*(U(HH+1-HH)+U(HH+1+HH)+2*sum(U(HH+1-HH+1:HH+1+HH-1)))/(2*h);
%    J(P)=(dx/2)*(U(HH+P-HH)+U(HH+P+HH)+2*sum(U(HH+P-HH+1:HH+P+HH-1)))/(2*h);
   for j=HH2+1:HH2+P
%      J2(j-HH2)=(dx/2)*(U2(j-HH2)+U2(j+HH2)+2*sum(U2(j-HH2+1:j+HH2-1)))/(2*h2);
 J2(j-HH2)=(dx/2)*(U2(j-HH2)+U2(j+HH2)+2*sum(U2(j-HH2+1:j+HH2-1)));
   end


    S10=N1-(dx/2)*(I1(1) +I1(P) + 2*sum(I1(2:P-1)));

    S20=N2-(dx/2)*(I2(1) +I2(P) + 2*sum(I2(2:P-1)));
 S101=(dx/2)*(I1(1) +I1(P) + 2*sum(I1(2:P-1)));
 S201=(dx/2)*(I2(1) +I2(P) + 2*sum(I2(2:P-1)));
    for i=1:P
    S1(i)=S10+S101-J1(i);

      S2(i)=S20+S201-J2(i);
    end


   
  

 I22(1)=I2(1)+2*r(1)*(I2(2)-2*I2(1)+I2(1))+dt*(beta2*S2(1)*I2(1)-(gamma2+nu2)*I2(1)); %Noflux boundary condition
 R22(1)=R2(1)+dt*(gamma2*I2(1));
 D22(1)=D2(1)+dt*(nu2*I2(1));
  
   
    for i=2:P-1
       I22(i)=I2(i)+r(i)*(I2(i+1)-2*I2(i)+I2(i-1))+dt*(beta2*S2(i)*I2(i)-(gamma2+nu2)*I2(i));
        R22(i)=R2(i)+dt*(gamma2*I2(i));
 D22(i)=D2(i)+dt*(nu2*I2(i));
    end
%       YY(P)=Y(P)+2*r*(Y(P-1)-Y(P))+dt*(mu*(1-p0-min(c*J(P),1-p0))-mu*Y(P)-beta*Z(P)*Y(P));
%        ZZ(P)=Z(P)+2*r*(Z(P-1)-Z(P))+dt*(beta*Z(P)*Y(P)-(mu+nu)*Z(P));
       I22(P)=I2(P)+2*r(P)*(I2(P-1)-2*I2(P)+I2(P))+dt*(beta2*S2(P)*I2(P)-(gamma2+nu2)*I2(P));
        R22(P)=R2(P)+dt*(gamma2*I2(P));
       D22(P)=D2(P)+dt*(nu2*I2(P));


        for i=1:P
      I11(i)=I1(i)+dt*(beta11*S1(i)*I1(i)+beta12*S1(i)*I2(i)-(gamma1+nu1)*I1(i));
      R11(i)=R1(i)+dt*(gamma1*I1(i));
      D11(i)=D1(i)+dt*(nu1*I1(i));
    end
%       Y=YY;
%       Z=ZZ;
%       Q=QQ; 

%     for j=1:P
%        if S1(j) < 0
%                 S1(j)=0;
%             end
%             if I11(j) < 0
%                 I11(j)=0;
%             end
%             if R11(j) < 0
%                 R11(j)=0;
%             end
%              if D11(j) < 0
%                 D11(j)=0;
%             end
%         
%     end
% 
%      for j=1:P
%          if S2(j) < 0
%                 S2(j)=0;
%             end
%         
%             if I22(j) < 0
%                 I22(j)=0;
%             end
%             if R22(j) < 0
%                 R22(j)=0;
%             end
%              if D22(j) < 0
%                 D22(j)=0;
%             end
%    
%     end


      I1=I11;
      R1=R11;
      D1=D11;
      I2=I22;
      R2=R22;
      D2=D22;
      

    
    if mod(n,10)==0
       CC=CC+1
        n=n
         subplot(2,1,1)
        plot(x,I1,'b-','LineWidth',2);
         subplot(2,1,2)
%         plot(x,S1,'b-','LineWidth',2); 
%         subplot(4,1,3)
         plot(x,I2,'r-','LineWidth',2);
%         subplot(4,1,4)
%         plot(x,S2,'r-','LineWidth',2); 
        pause(0.001)
        W11=[W11; S1];
        W12=[W12; I1];
        W13=[W13; R1];
        W14=[W14; D1];
        W21=[W21; S2];
        W22=[W22; I2];
        W23=[W23; R2];
        W24=[W24; D2];
%         WWW1(n)=sum(I1);
%         WWW2(n)=sum(I2);
    end
    
   
end


for j=1:length(W12(:,1))
WWW1(j)=sum(W12(j,:));
WWW2(j)=sum(W22(j,:));
end

pcolor(W12);
 shading flat



ttt(1)=0;
for j=2:length(WWW1)
ttt(j)=ttt(j-1)+10*dt;
end
surf(x,ttt,W12);shading flat



 surf(x,ttt,W12);shading flat
 xlabel('genotype $x$','interpreter','latex');
 ylabel('time $t$','interpreter','latex');
zlabel('$I_1(x,t)$','interpreter','latex');




 plot(ttt,WWW1,'b','LineWidth',2);hold on
xlabel('time $t$','interpreter','latex');
ylabel('$\bar{I}_1(t)$','interpreter','latex');

pcolor(x(1:end),ttt(1:6000),W12(1:6000,1:length(x)));shading flat
xlabel('genotype $x$','interpreter','latex');
 ylabel('time $t$','interpreter','latex');
