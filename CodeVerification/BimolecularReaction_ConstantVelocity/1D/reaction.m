function [c1,c2]=reaction(c01,c02,dt,Kr1,Kr2,I,N1,N2)

c01=c01/N1; c02=c02/N2;

c1(1)=c01(1); c1(I)=c01(I);
c2(1)=c02(1); c2(I)=c02(I);

c1(2:I-1)=c01(2:I-1)-dt*Kr1*c01(2:I-1).*c02(2:I-1).^2;
c2(2:I-1)=c02(2:I-1)-dt*Kr2*c01(2:I-1).*c02(2:I-1).^2;

c1=c1*N1; c2=c2*N2;