function plot_fig_Monod_Sat(x,y,c1,c2,c3,Vx,Vy)

figure;
mesh(x,y,c1); 
xlabel('$y$','Interpreter','latex'); ylabel('$x$','Interpreter','latex');
zlabel('$c_1(x,y,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;

figure;
mesh(x,y,c2); 
xlabel('$y$','Interpreter','latex'); ylabel('$x$','Interpreter','latex');
zlabel('$c_2(x,y,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;

figure;
mesh(x,y,c3); 
xlabel('$y$','Interpreter','latex'); ylabel('$x$','Interpreter','latex');
zlabel('$c_3(x,y,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;

figure
mesh(x,y,Vy)
xlabel('$y$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); 
zlabel('$q_x(x,y)$','Interpreter','latex'); view(115,15);

figure
mesh(x,y,Vx)
xlabel('$y$','Interpreter','latex'); ylabel('$x$','Interpreter','latex'); 
zlabel('$q_y(x,y)$','Interpreter','latex'); view(115,15);
