%% Contours - solutions of the Monod model
function subplot_surf_c(x,y,c1,c2,c3)

z=y;
levels=16;

figure;

subplot(1,3,1)
contourf(x,z,c1,levels,'EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$c_1(x,z,t)$','Interpreter','latex'); 

subplot(1,3,2)
contourf(x,z,c2,levels,'EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$c_2(x,z,t)$','Interpreter','latex'); 

subplot(1,3,3)
contourf(x,z,c3,levels,'EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
title('$c_3(x,z,t)$','Interpreter','latex'); 
