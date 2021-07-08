%% Contours - solutions of the Monod model
function subplot_surf(x,y,c1,c2,c3)

y1=0; y2=20; 

figure;

subplot(3,1,1)
surf(y,x,c1','EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
xlim([0 80]); ylim([y1 y2]); caxis([0 2]);
title('$c_1(x,y,t)$','Interpreter','latex'); 

subplot(3,1,2)
surf(y,x,c2','EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
xlim([0 80]); ylim([y1 y2]); caxis([0 5]);
title('$c_2(x,y,t)$','Interpreter','latex'); 

subplot(3,1,3)
surf(y,x,c3','EdgeColor','none'); colorbar; view(0,90)
set(gca, 'XTick', []); set(gca, 'YTick', []);
xlim([0 80]); ylim([y1 y2]);
title('$c_3(x,y,t)$','Interpreter','latex'); 
