%% 2D BGRW code for constant velocity and nolinear bimolecular reactions

clear all; 
close all
tic

Itest = 1;
L2_c1 = zeros(4,1); eoc_c1=zeros(3,1);
L2_c2 = zeros(4,1); eoc_c2=zeros(3,1);
for it = [10 20 40 80]% 160 320]
    %%   Grid Initialization
    I = it+1; J = 3*it/2+1;
    a=0; b=2;
    c=0; d=3;
    dx = (b-a)/(I-1);
    x = a:dx:b;
    dy=(d-c)/(J-1);
    y=c:dy:d;
    [X,Y] = meshgrid(x,y) ;    
    %% Parameters
    D1=0.1; D2=D1;
    Vy=-1*ones(J,I); Vx=zeros(J,I); 
    T = 1;
    Dfactor=2;
    dt=Dfactor*(2*D1/dx^2+2*D2/dy^2); dt=1./dt; 
    %% Exact solutions and sources
    solEc1 = @(t,x,y) x.*(2-x).*y.^3*exp(-0.1*t)/27;
    solEc2 = @(t,x,y) (x-1).^2.*y.^2*exp(-0.1*t)/9;
    % % with Matlab:
    % syms c1(t,x,y)
    % c1 = x*(2-x)*y^3*exp(-0.1*t)/27;
    % syms c2(t,x,y)
    % c2 = (x-1)^2*y^2*exp(-0.1*t)/9;
    % F1= diff(c1,t)-diff(c1,y)-0.1*diff(diff(c1,x),x)-0.1*diff(diff(c1,y),y)+1*c1*c2^2;
    % F2= diff(c2,t)-diff(c2,y)-0.1*diff(diff(c2,x),x)-0.1*diff(diff(c2,y),y)+2*c1*c2^2;
    F1 = @(t,x,y) (y.^3.*exp(-t./10))./135 + (x.*y.*exp(-t./10).*(x - 2))./45 + (x.*y.^2*exp(-t/10).*(x - 2))/9 ...
                  + (x.*y.^3.*exp(-t./10).*(x - 2))./270 - (x.*y.^7*exp(-t./5).*exp(-t./10).*(x - 1).^4.*(x - 2))./2187;
    F2 = @(t,x,y) - (y.^2.*exp(-t./10))./45 - (exp(-t./10).*(x - 1).^2)./45 - (2.*y.*exp(-t./10).*(x - 1).^2)./9 ...
                  - (y.^2.*exp(-t./10).*(x - 1).^2)./90 - (2.*x.*y.^7.*exp(-t./5).*exp(-t./10).*(x - 1).^4.*(x - 2))./2187;
    %%   Initial Conditions
    N=10^24; N1=N; N2=N1;
    c10=solEc1(0,X,Y)*N; c1=c10;
    c20=solEc2(0,X,Y)*N; c2=c20;
    %% solution
    t=0;
    while t<=T
        t=t+dt;        
        fg1=F1(t,X,Y)*N1;
        fg2=F2(t,X,Y)*N2;
        %% Transport step
        nspec=2; %  molecular_species
        for i=1:nspec
            switch i
                case 1
                    c10E=solEc1(t,X,Y)*N1;
                    [c1t]=BGRW_2D(c10,c10E,fg1,I,J,dx,dy,dt,Vx,Vy,D1,D2); % transport solver
                case 2
                    c20E=solEc2(t,X,Y)*N2;
                    [c2t]=BGRW_2D(c20,c20E,fg2,I,J,dx,dy,dt,Vx,Vy,D1,D2); % transport solver
            end
        end
        %% reaction
        Kr1=1; Kr2=2;
        [c1,c2]=reaction(c1t,c2t,dt,Kr1,Kr2,I,J,N1,N2);
        c10=c1; c20=c2;
        tE=t;
    end
    %% Results
    fprintf('The space step is : %0.4e \n',dx) ;
    fprintf('The time step is : %0.4e \n',dt) ;
    fprintf('The total time is : %0.4e \n',tE) ;
    
    c1E=solEc1(tE,X,Y); c1=c1/N1;
    L2_c1(Itest) = ( dx * dy )^(1/2) *norm(c1-c1E);  %
    c2E=solEc2(tE,X,Y); c2=c2/N2;
    L2_c2(Itest) = ( dx * dy )^(1/2) *norm(c2-c2E);  %
    if Itest >1
        eoc_c1(Itest-1)=log10(L2_c1(Itest-1)/L2_c1(Itest))/log10(2);
        eoc_c2(Itest-1)=log10(L2_c2(Itest-1)/L2_c2(Itest))/log10(2);
    end
    Itest = Itest+1; 
end
fprintf('L2_c1  : %0.4e \n',L2_c1)
fprintf('EOC_c1 : %0.4e \n',eoc_c1)
fprintf('L2_c2  : %0.4e \n',L2_c2)
fprintf('EOC_c2 : %0.4e \n',eoc_c2)

figure
cE1=solEc1(t,X,Y);
surf(c1-cE1)
xlabel('x'); ylabel('y'); zlabel('c1-cE1');
figure
cE2=solEc2(t,X,Y);
surf(c2-cE2)
xlabel('x'); ylabel('y'); zlabel('c2-cE2');

toc

%% main_2D_ReactTransp_test
% The space step is : 2.0000e-01 
% The time step is : 5.0000e-02 
% The total time is : 1.0000e+00 
% The space step is : 1.0000e-01 
% The time step is : 1.2500e-02 
% The total time is : 1.0125e+00 
% The space step is : 5.0000e-02 
% The time step is : 3.1250e-03 
% The total time is : 1.0000e+00 
% The space step is : 2.5000e-02 
% The time step is : 7.8125e-04 
% The total time is : 1.0008e+00 
% The space step is : 1.2500e-02 
% The time step is : 1.9531e-04 
% The total time is : 1.0000e+00 
% The space step is : 6.2500e-03 
% The time step is : 4.8828e-05 
% The total time is : 1.0000e+00 
% L2_c1  : 3.5266e-03 
% L2_c1  : 8.6394e-04 
% L2_c1  : 2.1245e-04 
% L2_c1  : 5.2980e-05 
% L2_c1  : 1.3224e-05 
% L2_c1  : 3.3044e-06
%                     EOC_c1 : 2.0293e+00 
%                     EOC_c1 : 2.0238e+00 
%                     EOC_c1 : 2.0036e+00 
%                     EOC_c1 : 2.0023e+00
%                     EOC_c1 : 2.0006e+00
% L2_c2  : 4.9122e-03
% L2_c2  : 1.1199e-03 
% L2_c2  : 2.6653e-04 
% L2_c2  : 6.4950e-05 
% L2_c2  : 1.6042e-05 
% L2_c2  : 3.9862e-06
%                     EOC_c2 : 2.1331e+00
%                     EOC_c2 : 2.0709e+00 
%                     EOC_c2 : 2.0369e+00 
%                     EOC_c2 : 2.0175e+00 
%                     EOC_c2 : 2.0088e+00
% Elapsed time is 6.700850 seconds (for it=[10 20 40 80]).
% Elapsed time is 32 min (for it=[10 20 40 80 160 320]).
