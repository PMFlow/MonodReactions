%% 1D noniterative BGRW scheme for constant velocity and nolinear bimolecular reactions

clear all; close all
tic
Itest = 1;
L2_c1 = zeros(4,1); eoc_c1=zeros(3,1);
L2_c2 = zeros(4,1); eoc_c2=zeros(3,1);
for it = [15 30 60 120 240 480] 
    %% Grid Initialization
    I = it+1;
    a=0; b=3;
    dx = (b-a)/(I-1);
    x = [a:dx:b];
    %% Parameters
    D1=0.1; 
    q=-1*ones(1,I);
    T = 1;
    Dfactor=2;
    dt=Dfactor*(2*D1/dx^2); dt=1./dt;
    %% Exact solutions and sources
    solEc1 = @(t,x) x.^3*exp(-0.1*t)/27;
    solEc2 = @(t,x) x.^2*exp(-0.1*t)/9;
    % % with Matlab:
    % syms c1(t,x)
    % c1 = x^3*exp(-0.1*t)/27;
    % syms c2(t,x)
    % c2 = x^2*exp(-0.1*t)/9;
    % F1= diff(c1,t)-diff(c1,x)-0.1*diff(diff(c1,x),x)+1*c1*c2^2;
    % F2= diff(c2,t)-diff(c2,x)-0.1*diff(diff(c2,x),x)+2*c1*c2^2;
    F1 = @(t,x) (exp(-t/5)*exp(-t/10).*x.^7)/2187 - (exp(-t/10).*x.^3)/270 - (exp(-t/10).*x.^2)/9 - (exp(-t/10).*x)/45;
    F2 = @(t,x) (2*exp(-t/5)*exp(-t/10).*x.^7)/2187 - (exp(-t/10).*x.^2)/90 - (2*exp(-t/10).*x)/9 - exp(-t/10)/45;
    %% Initial Conditions
    N=10^24;
    c10=solEc1(0,x)*N; c1=c10;
    c20=solEc2(0,x)*N; c2=c20;
    %% Solution
    t=0;
    c1a=c10; c2a=c20; N1=N; N2=N1;
    while t<=T
        t=t+dt;
        fg1=F1(t,x)*N1;
        fg2=F2(t,x)*N2;
        %% Transport step
        nspec=2; %  molecular_species
        for i=1:nspec           
            switch i
                case 1
                    c10E=solEc1(t,x)*N1;
                    [c1t]=BGRW_1D(c10,c10E,fg1,I,dx,dt,q,D1); % transport solver
                case 2
                    c20E=solEc2(t,x)*N2;
                    [c2t]=BGRW_1D(c20,c20E,fg2,I,dx,dt,q,D1); % transport solver
            end
        end
        %% reaction
        Kr1=1; Kr2=2;
        [c1,c2]=reaction(c1t,c2t,dt,Kr1,Kr2,I,N1,N2);
        c10=c1; c20=c2;
        tE=t;
    end
    %% Results
    fprintf('The space step is : %0.2e \n',dx) ;
    fprintf('The time step is : %0.2e \n',dt) ;
    fprintf('The total time is : %0.2e \n',tE) ;
    c1E=solEc1(tE,x); c1=c1/N1;
    L2_c1(Itest) = ( dx )^(1/2) *norm(c1-c1E); 
    c2E=solEc2(tE,x); c2=c2/N2;
    L2_c2(Itest) = ( dx )^(1/2) *norm(c2-c2E); 
    if Itest >1
        eoc_c1(Itest-1)=log10(L2_c1(Itest-1)/L2_c1(Itest))/log10(2);
        eoc_c2(Itest-1)=log10(L2_c2(Itest-1)/L2_c2(Itest))/log10(2);
    end
    Itest = Itest+1;
end
fprintf('L2_c1  : %0.2e \n',L2_c1)
fprintf('EOC_c1 : %0.2e \n',eoc_c1)
fprintf('L2_c2  : %0.2e \n',L2_c2)
fprintf('EOC_c2 : %0.2e \n',eoc_c2)
fprintf('grid_Peclet = %0.4e \n',dx*mean(abs(q))/D1);

plot_Figs_React(x,c1,c1E,c2,c2E)
toc
%% main_ReactTransp_test
% The space step is : 2.00e-01 
% The time step is : 1.00e-01 
% The total time is : 1.10e+00 
% The space step is : 1.00e-01 
% The time step is : 2.50e-02 
% The total time is : 1.00e+00 
% The space step is : 5.00e-02 
% The time step is : 6.25e-03 
% The total time is : 1.01e+00 
% The space step is : 2.50e-02 
% The time step is : 1.56e-03 
% The total time is : 1.00e+00 
% The space step is : 1.25e-02 
% The time step is : 3.91e-04 
% The total time is : 1.00e+00 
% The space step is : 6.25e-03 
% The time step is : 9.77e-05 
% The total time is : 1.00e+00 
% L2_c1  : 1.98e-02 
% L2_c1  : 4.87e-03 
% L2_c1  : 1.20e-03 
% L2_c1  : 3.01e-04 
% L2_c1  : 7.52e-05 
% L2_c1  : 1.88e-05 
%                  EOC_c1 : 2.02e+00 
%                  EOC_c1 : 2.01e+00 
%                  EOC_c1 : 2.00e+00 
%                  EOC_c1 : 2.00e+00 
%                  EOC_c1 : 2.00e+00 
% L2_c2  : 4.79e-02 
% L2_c2  : 1.17e-02 
% L2_c2  : 2.90e-03 
% L2_c2  : 7.23e-04 
% L2_c2  : 1.81e-04 
% L2_c2  : 4.52e-05 
%                  EOC_c2 : 2.04e+00 
%                  EOC_c2 : 2.01e+00 
%                  EOC_c2 : 2.00e+00 
%                  EOC_c2 : 2.00e+00 
%                  EOC_c2 : 2.00e+00 
% 
% Peclet = 2, 1, 0.5, 0.25, 0.125, 0.0625

% Elapsed time is 2.979036 seconds.
