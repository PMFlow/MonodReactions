function plot_conv(kt_plot,S,strvect)

load convf
iS=1:S;
dk=10; % 50;
figure; hold all;
for k=1:kt_plot
    P(k)=plot(iS(1:dk:end),convf(k,1:dk:end));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');


