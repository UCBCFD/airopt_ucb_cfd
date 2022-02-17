% Plot Morphing
figure
i_range = [1,5,10,20];
weight = [0.3,-0.2,0.4,0.4];
% - Adjust figure
set(gcf,'Color','White')
aspect_ratio = 2;
width = 0.18;
H = zeros(4001,2);
for j = 1:numel(i_range)
    i = i_range(j);
    M = readmatrix(sprintf('./BaseShapes/%d.txt', i));

    % SUBPLOT: Baseline Shapes
    %sp = subplot(numel(i_range)*2,2,4*j-3);
    sp = subplot(numel(i_range),4,4*j-3);
    plot(M(:,1),M(:,2));
    % - Adjust Sizes
    sp.Position(3) = width;
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==1
        set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','xcolor',[1 1 1])
        xlabel('Baseline Shapes')
    else
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1])%,'ycolor',[1 1 1])
    end
    box off
    ylabel(sprintf('Shape #%d',i_range(j)))

    % SUBPLOT: Baseline Coordinates
    %sp = subplot(numel(i_range)*2,2,4*j-1);
    sp = subplot(numel(i_range),4,4*j-2);    
    plot(M(:,2));
    xlim([0 4000]);
    % - Adjust Sizes
    sp.Position(3) = width;
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==1
        set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','ycolor',[1 1 1],'xcolor',[1 1 1])
        xlabel('Baseline Coordinates')
    else
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'ycolor',[1 1 1])
    end
    box off

    M(:,1) = M(:,1).*weight(j);
    M(:,2) = M(:,2).*weight(j);

    % SUBPLOT: Weighted Baseline Shapes
    %sp = subplot(numel(i_range)*2,2,4*j-2);
    sp = subplot(numel(i_range),4,4*j-1);    
    plot(M(:,1),M(:,2));
    % - Adjust Sizes
    sp.Position(3) = width*abs(weight(j));
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==1
        set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','xcolor',[1 1 1])%,'ycolor',[1 1 1])
        xlabel('Weighted Baseline Shapes')
    else
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1])%,'ycolor',[1 1 1])
    end
    ylabel(sprintf('Weight: %g',weight(j)))
    box off

    % SUBPLOT: Weighted Baseline Coordinates
    %sp = subplot(numel(i_range)*2,2,4*j-0);
    sp = subplot(numel(i_range),4,4*j-0);    
    plot(M(:,2));    
    xlim([0 4000]);
    % - Adjust Sizes
    sp.Position(3) = width*abs(weight(j));
    sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    if j==1
        set(gca,'XTick',[],'YTick',[],'XAxisLocation','top','ycolor',[1 1 1],'xcolor',[1 1 1])
        xlabel('Weighted Baseline Coordinates')
    else
        set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'ycolor',[1 1 1])
    end
    box off
    
    H = H + M;
end


% figure
% set(gcf,'Color','White')
% sgtitle('Morphed Shape')
% % SUBPLOT: Baseline Shapes
% sp = subplot(2,1,1);
% plot(H(:,1),H(:,2));
% % - Adjust Sizes
% sp.Position(3) = 0.8;
% sp.Position(4) = sp.Position(3)*(max(H(:,2))-min(H(:,2)))*aspect_ratio;
% set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'ycolor',[1 1 1])
% box off
% 
% % SUBPLOT: Baseline Coordinates
% sp = subplot(2,1,2);
% plot(H(:,2));
% xlim([0 4000]);
% % - Adjust Sizes
% sp.Position(3) = 0.8;
% sp.Position(4) = sp.Position(3)*(max(H(:,2))-min(H(:,2)))*aspect_ratio;
% set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'ycolor',[1 1 1])
% box off