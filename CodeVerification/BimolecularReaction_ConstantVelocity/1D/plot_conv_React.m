function plot_conv_React(kt_plot,S,strvect)

load convc1
iS=1:S;
figure; hold all;
for k=1:kt_plot
P(k)=plot(iS,convc1(k,:)); 
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|c_{1}^{s} - c_{1}^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');

load convc2
figure; hold all;
for k=1:kt_plot
P(k)=plot(iS,convc2(k,:));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|c_{2}^{s} - c_{2}^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');

