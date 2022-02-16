% Define Color
% Orange = [0.8500 0.3250 0.0980]	;
Orange = 'red';

% Plot Morphing
f=figure;
f.Position = [10,10,800,800];
i_range = [1,2,3,14,25];
weight = [0.6,-0.2,-0.3,0.5,0.5];
% - Adjust figure
set(gcf,'Color','White')
aspect_ratio = 2;
width = 0.2;
H = zeros(4001,2);
for j = 1:numel(i_range)
    i = i_range(j);
    M = readmatrix(sprintf('./BaseShapes/%d.txt', i));

    % SUBPLOT: Baseline Shapes
    sp = subplot(numel(i_range),3,3*j-2);
    scatter(M(:,1),M(:,2),5,'filled');
    % - Adjust Sizes
    sp.Position(3) = width;
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==1
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'YColor',[1 1 1],'fontsize',14)
        %         xlabel('Baseline Shapes')
    else
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'YColor',[1 1 1],'fontsize',14)
    end
    box off
    %     ylabel(sprintf('Shape #%d',i_range(j)))
    title(sprintf('S_{%g}',i))

    % SUBPLOT: Baseline Coordinates
    sp = subplot(numel(i_range),3,3*j-1);
    scatter(1:4001,M(:,2),5,'filled');
    xlim([0 4000]);
    % - Adjust Sizes
    sp.Position(3) = width;
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==4
        set(gca,'fontsize',14,'XTick',[],'YTick',[],'XAxisLocation','origin','YAxisLocation','origin')
        %         xlabel('Baseline Coordinates & Weighted Coordinates')
    else
        %         set(gca,'XTick',[],'fontsize',14,'XAxisLocation','origin','YAxisLocation','origin')
        set(gca,'fontsize',14,'XTick',[],'YTick',[],'XAxisLocation','origin','YAxisLocation','origin')
    end
    box off
    xlabel('i')
    ylabel('y')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','bottom')
    hXLabel = get(gca,'XLabel');
    set(hXLabel,'rotation',0,'HorizontalAlignment','center','VerticalAlignment','cap')

    % SUBPLOT: Baseline Coordinates
    sp = subplot(numel(i_range),3,3*j);
    s = scatter(1:4001,M(:,2),5,'filled');
    s.MarkerFaceAlpha = 0.004;
    xlim([0 4000]);
    % - Adjust Sizes
    sp.Position(3) = width;
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==4
        set(gca,'fontsize',14,'XTick',[],'YTick',[],'XAxisLocation','origin','YAxisLocation','origin')
        %         xlabel('Baseline Coordinates & Weighted Coordinates')
    else
        %         set(gca,'XTick',[],'YTick',[],'fontsize',14,'XAxisLocation','origin','YAxisLocation','origin')
        set(gca,'XTick',[],'YTick',[],'fontsize',14,'XAxisLocation','origin','YAxisLocation','origin')
    end
    box off

    M(:,1) = M(:,1).*weight(j);
    M(:,2) = M(:,2).*weight(j);
    hold on
    % SUBPLOT: Weighted Baseline Coordinates
    scatter(1:4001,M(:,2),5,Orange,'filled');
    hold off

    title(sprintf('w_{%g}={%g}',[i,weight(j)]))
    xlabel('i')
    ylabel('y')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','bottom')
    hXLabel = get(gca,'XLabel');
    set(hXLabel,'rotation',0,'HorizontalAlignment','center','VerticalAlignment','cap')

    H = H + M;
end

H = H./sum(weight);
f=figure;
f.Position = [10,10,800,800];
set(gcf,'Color','White')
% sgtitle('Morphed Shape')
% SUBPLOT: Morphed Shape
sp = subplot(2,1,2);
scatter(H(:,1),H(:,2),5,Orange,'filled');
% - Adjust Sizes
% sp.Position(2) = 0.6;
sp.Position(3) = 0.8;
sp.Position(4) = sp.Position(3)*(max(H(:,2))-min(H(:,2)))*aspect_ratio;
% set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','ycolor',[1 1 1],'xcolor',[1 1 1])
set(gca,'XTick',[],'YTick',[],'ycolor',[1 1 1],'xcolor',[1 1 1],...
    'fontsize',14)
% xlabel('Morphed Shape')
title('Morphed Shape')
box off

% SUBPLOT: Morphed Coordinates
sp = subplot(2,1,1);
scatter(1:4001,H(:,2),5,Orange,'filled');
xlim([0 4000]);
% - Adjust Sizes
% sp.Position(2) = 0.6;
sp.Position(3) = 0.8;
sp.Position(4) = sp.Position(3)*(max(H(:,2))-min(H(:,2)))*aspect_ratio;
% set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','ycolor',[1 1 1],'xcolor',[1 1 1])
% set(gca,'XTick',[],'YTick',[],'ycolor',[1 1 1],'xcolor',[1 1 1],...
%     'fontsize',14)
set(gca,'fontsize',14,'XTick',[],'YTick',[],'XAxisLocation','origin','YAxisLocation','origin')
% xlabel('Morphed Coordinates')
title('Morphed Coordinates')
box off
xlabel('i')
ylabel('y')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','bottom')
hXLabel = get(gca,'XLabel');
set(hXLabel,'rotation',0,'HorizontalAlignment','center','VerticalAlignment','cap')