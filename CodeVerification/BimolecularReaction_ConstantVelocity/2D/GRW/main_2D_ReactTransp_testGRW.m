%% 2D unbiased-GRW code for constant velocity and nolinear bimolecular reactions

clear all; 
close all
tic

Itest = 1;
L2_c1 = zeros(4,1); eoc_c1=zeros(3,1);
L2_c2 = zeros(4,1); eoc_c2=zeros(3,1);
for it = [10 20 40 80 160 320] 
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
    D1=0.0001; D2=D1; tht=1;
    Vy=-1*ones(J,I); Vx=zeros(J,I); 
    T = 1;
    stepU=1;
    U_MEAN=-1;
    dt=stepU*dx/abs(U_MEAN);
    %% Exact solutions and sources
    solEc1 = @(t,x,y) x.*(2-x).*y.^3*exp(-0.1*t)/27;
    solEc2 = @(t,x,y) (x-1).^2.*y.^2*exp(-0.1*t)/9;
    % % with Matlab:
    % syms c1(t,x,y)
    % c1 = x*(2-x)*y^3*exp(-0.1*t)/27;
    % syms c2(t,x,y)
    % c2 = (x-1)^2*y^2*exp(-0.1*t)/9;
    % F1=
    % diff(c1,t)-diff(c1,y)-0.0001*diff(diff(c1,x),x)-0.001*diff(diff(c1,y),y)+1*c1*c2^2;  % -0.1*diff = orig
    % F2= diff(c2,t)-diff(c2,y)-0.0001*diff(diff(c2,x),x)-0.001*diff(diff(c2,y),y)+2*c1*c2^2;
    F1 = @(t,x,y) (y.^3.*exp(-t./10))./135/1000 + (x.*y.*exp(-t./10).*(x - 2))./45/1000 + (x.*y.^2*exp(-t/10).*(x - 2))/9 ...
                  + (x.*y.^3.*exp(-t./10).*(x - 2))./270 - (x.*y.^7*exp(-t./5).*exp(-t./10).*(x - 1).^4.*(x - 2))./2187;
    F2 = @(t,x,y) - (y.^2.*exp(-t./10))./45/1000 - (exp(-t./10).*(x - 1).^2)./45/1000 - (2.*y.*exp(-t./10).*(x - 1).^2)./9 ...
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
                    [c1t]=GRW_2D(c10,c10E,fg1,tht,I,J,dx,dy,dt,Vx,Vy,D1,D2,stepU);
                case 2
                    c20E=solEc2(t,X,Y)*N2;
                    [c2t]=GRW_2D(c20,c20E,fg2,tht,I,J,dx,dy,dt,Vx,Vy,D1,D2,stepU);
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
    L2_c1(Itest) = ( dx * dy )^(1/2) *norm(c1-c1E);  
    c2E=solEc2(tE,X,Y); c2=c2/N2;
    L2_c2(Itest) = ( dx * dy )^(1/2) *norm(c2-c2E);  
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
mesh(c1-cE1)
xlabel('x'); ylabel('y'); zlabel('c1-cE1');
figure
cE2=solEc2(t,X,Y);
mesh(c2-cE2)
xlabel('x'); ylabel('y'); zlabel('c2-cE2');

toc

%% main_2D_ReactTransp_testGRW
% The space step is : 2.0000e-01 
% The time step is : 2.0000e-01 
% The total time is : 1.2000e+00 
% The space step is : 1.0000e-01 
% The time step is : 1.0000e-01 
% The total time is : 1.1000e+00 
% The space step is : 5.0000e-02 
% The time step is : 5.0000e-02 
% The total time is : 1.0000e+00 
% The space step is : 2.5000e-02 
% The time step is : 2.5000e-02 
% The total time is : 1.0000e+00 
% The space step is : 1.2500e-02 
% The time step is : 1.2500e-02 
% The total time is : 1.0125e+00 
% The space step is : 6.2500e-03 
% The time step is : 6.2500e-03 
% The total time is : 1.0062e+00 
% L2_c1  : 1.6882e-01 
% L2_c1  : 8.4264e-02 
% L2_c1  : 4.1300e-02 
% L2_c1  : 2.0701e-02 
% L2_c1  : 1.0370e-02 
% L2_c1  : 5.1366e-03 
%                     EOC_c1 : 1.0025e+00 
%                     EOC_c1 : 1.0288e+00 
%                     EOC_c1 : 9.9645e-01 
%                     EOC_c1 : 9.9720e-01 
%                     EOC_c1 : 1.0136e+00 
% L2_c2  : 7.7225e-02 
% L2_c2  : 3.0105e-02 
% L2_c2  : 1.2747e-02 
% L2_c2  : 5.9591e-03 
% L2_c2  : 2.9070e-03 
% L2_c2  : 1.4198e-03 
%                     EOC_c2 : 1.3591e+00 
%                     EOC_c2 : 1.2399e+00 
%                     EOC_c2 : 1.0969e+00 
%                     EOC_c2 : 1.0356e+00 
%                     EOC_c2 : 1.0339e+00 
% Elapsed time is 12.081652 seconds.