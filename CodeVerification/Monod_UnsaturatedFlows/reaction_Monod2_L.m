function [c1,c2]=reaction_Monod2_L(c01,c02,dt,tht,N1,N2,Lc)

a1=1; a2=3; M1=2; M2=0.2;
c01=c01/N1; c02=c02/N2;

m=1e-3*(c01./(M1+c01)).*(c02./(M2+c02));
c1=c01-dt*a1*m.*tht/Lc;
c2=c02-dt*a2*m.*tht/Lc;

c1=c1*N1; c2=c2*N2;