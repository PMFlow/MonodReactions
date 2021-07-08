function [Vx,Vy] = velocity(I,J,dx,dy,p0,D,DK)
Vx=zeros(J,I); Vy=zeros(J,I);
%       V_interior of \Omega        
        Vx(2:J-1,2:I-1)=-D.*((p0(2:J-1,3:I)-p0(2:J-1,1:I-2))/(2*dx)); 
        Vy(2:J-1,2:I-1)=-D.*((p0(3:J,2:I-1)-p0(1:J-2,2:I-1))/(2*dy)+1); 
%       V_normal to boundary approximated by finite-differences
        Vx(:,1)=-DK(:,1).*(p0(:,2)-p0(:,1))/dx; Vx(:,I)=-DK(:,I).*(p0(:,I)-p0(:,I-1))/dx; 
        Vx(1,2:I-1)=-DK(1,2:I-1).*(p0(1,3:I)-p0(1,1:I-2))/(2*dy); 
        Vx(1,1)=-DK(1,1).*(p0(1,2)-p0(1,1))/dy; 
        Vx(1,I)=-DK(1,I-1).*(p0(1,I)-p0(1,I-1))/dy; 
        Vx(J,2:I-1)=-DK(J,2:I-1).*(p0(J,3:I)-p0(J,1:I-2))/(2*dy); 
        Vx(J,1)=-DK(J,1).*(p0(J,2)-p0(J,1))/dy; 
        Vx(J,I)=-DK(J,I-1).*(p0(J,I)-p0(J,I-1))/dy; 
        Vy(1,:)=-DK(1,:).*((p0(2,:)-p0(1,:))/dy+1); Vy(J,:)=-DK(J,:).*((p0(J,:)-p0(J-1,:))/dy+1);
        Vy(2:J-1,1)=-DK(2:J-1,1).*((p0(3:J,1)-p0(1:J-2,1))/(2*dy)+1);
        Vy(1,1)=-DK(1,1).*((p0(2,1)-p0(1,1))/dy+1);
        Vy(J,1)=-DK(J-1,1).*((p0(J,1)-p0(J-1,1))/dy+1);
        Vy(2:J-1,I)=-DK(2:J-1,I).*((p0(3:J,I)-p0(1:J-2,I))/(2*dy)+1);
        Vy(1,I)=-DK(1,I).*((p0(2,I)-p0(1,I))/dy+1);
        Vy(J,I)=-DK(J-1,I).*((p0(J,I)-p0(J-1,I))/dy+1);
%       V_normal to boundary extended from interior of \Omega
%         Vx(:,1)=Vx(:,2); Vx(:,I)=Vx(:,I-1); 
%         Vx(1,:)=Vx(2,:); Vx(J,:)=Vx(J-1,:); 
%         Vy(:,1)=Vy(:,2); Vy(:,I)=Vy(:,I-1);
%         Vy(1,:)=Vy(2,:); Vy(J,:)=Vy(J-1,:); 
