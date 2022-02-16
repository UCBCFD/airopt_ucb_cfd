% Plot Repaneling
f=figure;
% f.Position = [10,10,800,800];
% - Adjust figure
% set(gcf,'Color','White')
aspect_ratio = 1;
H = zeros(4001,2);

i = 25;

M = readmatrix(sprintf('./BaseShapes/%d.txt', i));
subplot(2,1,1)
plot(M(:,1),M(:,2),'LineWidth',2);
% hold on
set(gca,'XTick',[],'YTick',[],'XAxisLocation','top',...
    'fontsize',14,'ycolor',[1 1 1],'XColor',[1 1 1],'color','none')
axis off
box off
xlabel('off')
ylabel('off')
daspect([1 1 1])

M_copy = M;
[M_2] = Panelling(M_copy,50);

sp = subplot(2,1,2);
plot_dir(M(:,1),M(:,2),400,20,10,10,'r');
% set(gca,'XTick',[],'YTick',[],'XAxisLocation','top',...
%     'fontsize',14,'ycolor',[1 1 1],'XColor',[1 1 1],'color','none')
% axis off
% box off
% daspect([1 1 1])
hold on
% subplot(2,1,3)
scatter(M_2(:,1),M_2(:,2),20,[0, 0.4470, 0.7410],'filled')%,'MarkerFaceAlpha',0.6);

set(gca,'XTick',[],'YTick',[],'fontsize',14,'ycolor',[1 1 1],'XColor',[1 1 1],'color','none')
axis off
% box off
daspect([1 1 1])
sp.Position(4) = sp.Position(4) + 0.3;
sp.Position(2) = sp.Position(2) - 0.1;