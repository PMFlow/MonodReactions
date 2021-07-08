%% Contours of flow solutions
function subplot_surf_ptht(x,y,p,tht,Vx,Vy)

z=y;
levels=16;

figure;

subplot(1,3,1)
contourf(x,z,p,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$\psi(x,z,t)$','Interpreter','latex'); 

subplot(1,3,2)
contourf(x,z,tht,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$\theta(x,z,t)$','Interpreter','latex'); 


subplot(1,3,3)
contourf(x,z,sqrt(Vx.^2+Vy.^2),levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$V(x,z,t)$','Interpreter','latex'); 
