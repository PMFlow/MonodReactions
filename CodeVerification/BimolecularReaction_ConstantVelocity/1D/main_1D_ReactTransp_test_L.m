%% 1D iterative BGRW scheme for constant velocity and nolinear bimolecular reactions

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
    past=T/3;
    S= 10000; 
    Tolerance = 1e-6;
    Lc=1; 
    Dfactor=2; 
    dtc=Dfactor*(2*D1/Lc/dx^2); dtc=1./dtc; dt=dtc;
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
    c10=solEc1(0,x)*N; c1=c10; c10t=c10;
    c20=solEc2(0,x)*N; c2=c20; c20t=c20;
    c10E=c10; c20E=c20;
    %% Solution
    tgraf=0; t=0; kt=1; test_kt=0;
    epsc1=zeros(1,S); epsc2=zeros(1,S); convc1=zeros(3,S); convc2=zeros(3,S);
    c1a=c10/N; c2a=c20/N; N1=N; N2=N1;
    while t<=T
        t=t+dt;        
        fg1=F1(t,x)*N1;
        fg2=F2(t,x)*N2;       
        for s=1:S % iterations of the L-scheme
            %% Transport and reaction steps
            nspec=2; %  molecular_species
            Kr1=1; Kr2=2;
            for i=1:nspec
                switch i
                    case 1
                        dc1=(c10t-c1);
                        c10E=solEc1(t,x)*N1;
                        [c1t]=BGRW_1D_L(c10,c10E,dc1,fg1,I,dx,dt,Lc,q,D1); % transport solver
                        [c1,~]=reaction(c1t,c20,dt,Kr1,Kr2,I,N1,N2); % reaction
                        c10=c1;
                    case 2
                        dc2=(c20t-c2);
                        c20E=solEc2(t,x)*N2;
                        [c2t]=BGRW_1D_L(c20,c20E,dc2,fg2,I,dx,dt,Lc,q,D1); % transport solver
                        [~,c2]=reaction(c10,c2t,dt,Kr1,Kr2,I,N1,N2); % reaction
                        c20=c2;
                end
            end
            %% Convergence criterion
            tol_epsc1=norm(c1/N1-c1a);
            tol_epsc2=norm(c2/N2-c2a);
            if kt*past>=t && kt*past<t+dt && t<=T
                epsc1(s)=tol_epsc1;
                epsc2(s)=tol_epsc2;
            end
            if max(tol_epsc1,tol_epsc2) <= Tolerance
                break
            end
            c1a=c1/N1; c2a=c2/N2;
        end
        if  kt*past>=t && kt*past<t+dt && t<=T
            tgraf=tgraf+1;
            rndt=kt*past;
            str=['t=',num2str(rndt)];
            strvect(tgraf,1:length(str))=str;
            convc1(kt,:)=epsc1;
            convc2(kt,:)=epsc2;
            kt=kt+1;
            test_kt=1;
        else
            test_kt=0;
        end
        tE=t;
        c10t=c10; c20t=c20;        
    end
    if  kt==4
        save('convc1');
        save('convc2');
    end
    %% Results
    kt_plot=kt-1;
    plot_conv_React(kt_plot,S,strvect)
    fprintf('The space step is : %0.2e \n',dx) ;
    fprintf('The time step is : %0.2e \n',dt) ;
    fprintf('The total time is : %0.2e \n',tE) ;    
    c1E=solEc1(tE,x); c1=c1/N1;
    L2_c1(Itest) = ( dx )^(1/2) *norm(c1-c1E);  %
    c2E=solEc2(tE,x); c2=c2/N2;
    L2_c2(Itest) = ( dx )^(1/2) *norm(c2-c2E);  %
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

%% main_ReactTransp_test_L
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
% L2_c1  : 1.20e-02 
% L2_c1  : 2.32e-03 
% L2_c1  : 5.56e-04 
% L2_c1  : 1.36e-04 
% L2_c1  : 3.40e-05 
% L2_c1  : 8.48e-06 
                    % EOC_c1 : 2.37e+00 
                    % EOC_c1 : 2.06e+00 
                    % EOC_c1 : 2.03e+00 
                    % EOC_c1 : 2.00e+00 
                    % EOC_c1 : 2.00e+00 
% L2_c2  : 4.57e-02 
% L2_c2  : 1.12e-02 
% L2_c2  : 2.81e-03 
% L2_c2  : 7.02e-04 
% L2_c2  : 1.76e-04 
% L2_c2  : 4.39e-05 
                    % EOC_c2 : 2.03e+00 
                    % EOC_c2 : 2.00e+00 
                    % EOC_c2 : 2.00e+00 
                    % EOC_c2 : 2.00e+00 
                    % EOC_c2 : 2.00e+00 
% grid_Peclet = 6.2500e-02 
% Current plot held
% Elapsed time is 6.208576 seconds.
