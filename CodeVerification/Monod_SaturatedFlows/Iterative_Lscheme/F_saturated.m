function [fs,fg,fh,VxE,VyE]=F_saturated(t,x,y)

[X,Y] = meshgrid(x,y) ; 
    % % with Matlab:
%     syms p(t,x,y)
%     p = t*x*(2-x)*y*(3-y);
%     syms c1(t,x,y)
%     c1 = x*(2-x)*y^3*exp(-0.1*t)/27;
%     syms c2(t,x,y)
%     c2 = (x-1)^2*y^2*exp(-0.1*t)/9;
%     theta = 1/(3.3333-p);
%     K= 0.05;
%     theta_sat=0.3;
%     D1=0.025; D2=D1;
%     Vx=-K*(diff(p, x))
%     Vx= - (t*x*y*(y - 3))/20 - (t*y*(x - 2)*(y - 3))/20
      solVx = @(t,x,y) -(t*x.*y.*(y - 3))/20 - (t*y.*(x - 2).*(y - 3))/20;
%     Vy =-K*(diff(p+y, y))
%     Vy=(t*x*y*(x - 2))/20 + (t*x*(x - 2)*(y - 3))/20 - 1/20
      solVy = @(t,x,y) -(t*x.*y.*(x - 2))/20 - (t*x.*(x - 2).*(y - 3))/20 - 1/20;
      VxE = solVx(t,X,Y);
      VyE = solVy(t,X,Y);    
%     F2 = diff(theta_sat,t)+diff(Vx,x)+diff(Vy,y); %% Vx=-K*diff(p,x); Vy=-K*diff(p+y,y) 
      F2 = @(t,x,y) - (t*x.*(x - 2))/10 - (t*y.*(y - 3))/10;

% % a1=1; a2=3; M1=2; M2=0.2;
% % m=1e-3*(c1/(2+c1))*(c2/(0.2+c2));
% % c1(2:J-1,2:I-1)=c01(2:J-1,2:I-1)-dt*a1*m(2:J-1,2:I-1);
% % c2(2:J-1,2:I-1)=c02(2:J-1,2:I-1)-dt*a2*m(2:J-1,2:I-1);
%     m = (x*y^5*exp(-t/5)*(x - 1)^2*(x - 2))/(243000*((y^2*exp(-t/10)*(x - 1)^2)/9 + 1/5)*((x*exp(-t/10)*(x - 2)*y^3)/27 - 2))
%     R1=-m*theta_sat;
%     R2=-3*m*theta_sat;
%     G2 = diff(theta_sat*c1,t)+diff(Vx*c1,x)+diff(Vy*c1,y)-D1*diff(diff(c1,x),x)-D2*diff(diff(c1,y),y)-R1;
      G2 = @(t,x,y) (y.^3*exp(-t/10))/540 + (y.^3*exp(-t/10).*(x - 2).*((t*x.*y.*(y - 3))/20 + (t*y.*(x - 2).*(y - 3))/20))/27 + (x.*y*exp(...
          -t/10).*(x - 2))/180 + (x.*y.^3*exp(-t/10).*(x - 2))/900 + (x.*y.^3*exp(-t/10).*((t*x.*y.*(y - 3))/20 + (t*y.*(x - 2).*(y ...
          - 3))/20))/27 + (x.*y.^2*exp(-t/10).*(x - 2).*((t*x.*y.*(x - 2))/20 + (t*x.*(x - 2).*(y - 3))/20 + 1/20))/9 + (t*x.^2.*y.^3*exp(...
          -t/10).*(x - 2).^2)/270 + (t*x.*y.^4*exp(-t/10).*(x - 2).*(y - 3))/270 + (x.*y.^5*exp(-t/5).*(x - 1).^2.*(x - 2))./(810000*((y.^2*exp(...
          -t/10).*(x - 1).^2)/9 + 1/5).*((x*exp(-t/10).*(x - 2).*y.^3)/27 - 2));
 
%     H2 = diff(theta_sat*c2,t)+diff(Vx*c2,x)+diff(Vy*c2,y)-D1*diff(diff(c2,x),x)-D2*diff(diff(c2,y),y)-R2;
      H2 = @(t,x,y) (x.*y.^5*exp(-t/5).*(x - 1).^2.*(x - 2))./(270000*((y.^2*exp(-t/10).*(x - 1).^2)/9 + 1/5).*((x*exp(-t/10).*(x ...
          - 2).*y.^3)/27 - 2)) - (exp(-t/10).*(x - 1).^2)/180 - (y.^2*exp(-t/10).*(x - 1).^2)/300 - (2*y*exp(-t/10).*(x - 1).^2.*((t*x.*y.*(x ...
          - 2))/20 + (t*x.*(x - 2).*(y - 3))/20 + 1/20))/9 - (y.^2*exp(-t/10).*(2*x - 2).*((t*x.*y.*(y - 3))/20 + (t*y.*(x - 2).*(y ...
          - 3))/20))/9 - (t*y.^3*exp(-t/10).*(x - 1).^2.*(y - 3))/90 - (t*x.*y.^2*exp(-t/10).*(x - 1).^2.*(x - 2))/90 - (y.^2*exp(-t/10))/180;
     
    fs=F2(t,X,Y); 
    fg=G2(t,X,Y); 
    fh=H2(t,X,Y); 
