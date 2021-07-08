function plot_Figs_React(x,c1,c1E,c2,c2E)

dx=10;
figure;
plot(x(1:dx:end),c1(1:dx:end),'*r',x,c1E,'b'); 
xlabel('$z$','Interpreter','latex');
ylabel('$c_1(z,t)$','Interpreter','latex'); 
legend('numeric','exact');

figure;
plot(x(1:dx:end),c2(1:dx:end),'*r',x,c2E,'b'); 
xlabel('$z$','Interpreter','latex');
ylabel('$c_2(z,t)$','Interpreter','latex'); 
legend('numeric','exact');

figure
plot(x,c1-c1E)
hold
plot(x,c2-c2E)
xlabel('$z$','Interpreter','latex');
ylabel('$c_1(z,t)$','Interpreter','latex'); 
legend('$c_1-c_{1E}$','$c_2-c_{2E}$','Interpreter','latex');
