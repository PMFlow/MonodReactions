%% Contours of solutions to the coupled flow-transport in the (x,z) plane

close all; clear all;

N=10^24;
y1=1; y2=20; 
kp=1;

figure;
x0=2;
y0=5;
width=30;
height=6;
set(gcf,'units','centimeters','position',[x0,y0,width,height])

subplotLayout = [2,3];
ax = gobjects(fliplr(subplotLayout));
for i = 1:prod(subplotLayout)
    ax(i) = subplot(subplotLayout(1),subplotLayout(2),i); 
end
ax = ax'; 
titles = {{'$BGRW$','$c_2(x,y,t)$'}, {'$GRW$','$c_2(x,y,t)$'}};  
centerSubHandles = ax(:,ceil(subplotLayout(2)/2));

for i=1:2
    switch i
        case 1
            load bgrwT300.mat
        case 2
            load grwT300.mat
    end

subplot(2,3,kp)
surf(y,x,c1'/N,'EdgeColor','none'); colorbar; view(0,90)
xlim([0 80]); ylim([y1 y2]); caxis([0 2]);
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('$c_1(x,y,t)$','Interpreter','latex'); 
kp=kp+1;

subplot(2,3,kp)
surf(y,x,c2'/N,'EdgeColor','none'); colorbar; view(0,90)
xlim([0 80]); ylim([y1 y2]); caxis([0 5]);
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('$c_2(x,y,t)$','Interpreter','latex'); 
kp=kp+1;

subplot(2,3,kp)
surf(y,x,c3','EdgeColor','none'); colorbar; view(0,90)
xlim([0 80]); ylim([y1 y2]); 
xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
title('$c_3(x,y,t)$','Interpreter','latex'); 
kp=kp+1;

end
for j = 1:numel(centerSubHandles)
    title(centerSubHandles(j), titles{j},'Interpreter','latex')
end