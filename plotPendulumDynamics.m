numSamp = 66;
hold on
phi0 = [0;pi/4];
phi = zeros(2,numSamp);
phi(:,1) = phi0; 
for ii = 2:length(phi)
    phi(:,ii) = pendulumDynamics(phi(:,ii-1),0.1);
end

plot(phi(1,:),phi(2,:),'b','linewidth',3)
plot(phi(1,1:2:end),phi(2,1:2:end),'s','color',[0.4940 0.1840 0.5560],'markerSize',5,'linewidth',2)

numSamp = 78;
hold on
phi0 = [0.5;pi/2];
phi = zeros(2,numSamp);
phi(:,1) = phi0; 
for ii = 2:length(phi)
    phi(:,ii) = pendulumDynamics(phi(:,ii-1),0.1);
end

plot(phi(1,:),phi(2,:),'b','linewidth',3)
plot(phi(1,1:2:end),phi(2,1:2:end),'.','color',[0.8500 0.3250 0.0980],'markerSize',20)
xlabel('$\phi_1$','interpreter','latex','fontsize',30)
ylabel('$\phi_2$','fontsize',30,'interpreter','latex')
set(gcf, 'Position',  [100, 100, 800, 700])
set(gca,'FontSize',20)