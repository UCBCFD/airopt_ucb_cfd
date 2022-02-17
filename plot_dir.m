function [h] = plot_dir (vX, vY, data_skip, data_span, headWidth, headLength, color)
%function [h1, h2] = plot_dir (vX, vY)
%Plotting x-y variables with direction indicating vector to the next element.
%Example
%   vX = linspace(0,2*pi, 10)';
%   vY = sin (vX);
%   plot_dir(vX, vY);

% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = data_skip:data_skip:(round(lenTime/2)-data_skip);
vSelect0 = [vSelect0,round(lenTime/2)+data_skip:data_skip:(lenTime-data_skip)];
% vSelect0 = data_skip:data_skip:(lenTime-data_skip);
% Indices of tails of arrows
vSelect1 = vSelect0 + data_span;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% Norm
norm = ((vYQ1-vYQ0).^2+(vXQ1-vXQ0).^2).^.5;
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0)./norm ;
vPy = (vYQ1 - vYQ0)./norm ;
% make plot
% h1 = plot (vX, vY, '.-'); hold on;
% add arrows
h = quiver (vXQ0,vYQ0, vPx, vPy, 0);

% Adjust arrowhead
%get the data from regular quiver
U = h.UData;
V = h.VData;
X = h.XData;
Y = h.YData;
hold on;
for ii = 1:length(X)
    LineLength = 0.001;
    ah = annotation('arrow',...
        'headStyle','cback1','HeadLength',headLength,'HeadWidth',headWidth,'Color',color);
    set(ah,'parent',gca);
    %         set(ah,'position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);
    set(ah,'position',[X(ii) Y(ii) LineLength*U(ii) LineLength*V(ii)]);
end

h.Visible = 'off';

% grid on; %hold off
% axis equal

end