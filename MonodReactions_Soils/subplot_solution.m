%% Contours of solutions to Monod model
% function subplot_solution(x,y,c1,c2,c3)

levels=18;
% itest=1 - loam; itest=2 - clay

for itest=1:2
    kp=1;
    
    figure()
    subplotLayout = [3,5];
    ax = gobjects(fliplr(subplotLayout));
    for i = 1:prod(subplotLayout)
        ax(i) = subplot(subplotLayout(1),subplotLayout(2),i);
    end
    ax = ax';
    % list all titles in order (top to bottom)
    titles = {{'$T=1$','$c_3(x,z,t)$'}, {'$T=3$','$c_3(x,z,t)$'}, {'$T=5$','$c_3(x,z,t)$'}};
    centerSubHandles = ax(:,ceil(subplotLayout(2)/2));
    for i=1:3
        switch i
            case 1
                if itest==2
                    load clayT1
                else
                    load loamT1
                end
            case 2
                if itest==2
                    load clayT3
                else
                    load loamT3
                end
            case 3
                if itest==2
                    load clayT5
                else
                    load loamT5
                end
        end
        
        z=y;
        subplot(3,5,kp)
        contourf(x,z,c1/N,levels,'EdgeColor','none'); colorbar; view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []); caxis([0 2]);
        title('$c_1(x,z,t)$','Interpreter','latex');
        kp=kp+1;
        subplot(3,5,kp)
        contourf(x,z,c2/N,levels,'EdgeColor','none'); colorbar; view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []); caxis([0 5]);
        title('$c_2(x,z,t)$','Interpreter','latex');
        kp=kp+1;
        subplot(3,5,kp)
        contourf(x,z,c3,levels,'EdgeColor','none'); colorbar; view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []);
        kp=kp+1;
        subplot(3,5,kp)
        contourf(x,z,p,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []);
        title('$\psi(x,z,t)$','Interpreter','latex');
        kp=kp+1;
        subplot(3,5,kp)
        contourf(x,z,tht,levels,'EdgeColor','none'); colorbar; colormap(parula); view(0,90)
        set(gca, 'XTick', []); set(gca, 'YTick', []);
        title('$\theta(x,z,t)$','Interpreter','latex');
        kp=kp+1;
    end
    for j = 1:numel(centerSubHandles)
        title(centerSubHandles(j), titles{j},'Interpreter','latex')
    end
end