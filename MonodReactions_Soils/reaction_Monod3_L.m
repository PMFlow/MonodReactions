function [c1,c2,c3]=reaction_Monod3_L(c01,c02,c03,dt,tht,Lc,N1,N2)

a1=1; a2=3; M1=2; M2=0.2; Y=0.09; Kd=0.05; maxm=5; 

maxc3=1;
c01=c01/N1; c02=c02/N2;

m=maxm*(c01./(M1+c01)).*(c02./(M2+c02)).*c03;
c1=c01-dt*a1*m.*tht/Lc;
c2=c02-dt*a2*m.*tht/Lc;
c3=c03+dt*(Y*m.*(1-c03./maxc3)-Kd*c03);

c1=c1*N1; c2=c2*N2;