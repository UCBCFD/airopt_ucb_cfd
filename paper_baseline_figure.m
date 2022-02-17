% Define Color
Orange = [0.8500 0.3250 0.0980]	;

% Plot Morphing
f=figure;
% f.Position = [10,10,800,800];
i_range = 1:25;
% - Adjust figure
set(gcf,'Color','White')
aspect_ratio = 2;
width = 0.15;
for j = 1:numel(i_range)
    i = i_range(j);
    M = readmatrix(sprintf('./BaseShapes/%d.txt', i));

    % SUBPLOT: Baseline Shapes
    sp = subplot(3,9,j);
    scatter(M(:,1),M(:,2),5,'filled');
    % - Adjust Sizes
%     sp.Position(3) = width;
%     sp.Position(4) = sp.Position(3)*(max(M(:,2))-min(M(:,2)))*aspect_ratio;
    set(gca,'XTick',[],'YTick',[],'xcolor',[1 1 1],'YColor',[1 1 1],'fontsize',14)
    box off
    daspect([1 1 1])
    %     ylabel(sprintf('Shape #%d',i_range(j)))
    %     title(sprintf('S_{%g}',i))

end
