function [morphed2]=Smoothing(morphed,angle_threshold,level,display_switch)

% morphed = importdata(filename);
coord = morphed(:,1)+morphed(:,2).*1i;
morphed2 = morphed;
vec = diff(coord);
ang = wrapToPi(diff(angle(vec))+pi); %angle formed by neighboring 3 pts
ang2 = abs(rad2deg(ang));
edge = find(ang2<angle_threshold)+1; %collocation pt is offset by 1

for i=1:size(edge)
    lowbnd = max(edge(i)-10,1);
    upbnd = min(edge(i)+10,size(coord,1));
   subcoord = morphed(lowbnd:upbnd,:);
   for j = 1:level
       subcoord = smoothdata(subcoord,'gaussian',upbnd-lowbnd+1);
       subcoord = smoothdata(subcoord,'sgolay',upbnd-lowbnd+1);
   end
   morphed2(lowbnd+2:upbnd-2,:) = subcoord(3:end-2,:);
end

if nargin < 4
    display_switch = 0;
else
    display_switch = 1;
end
if display_switch == 1
%     scatter(subcoord_s.x,subcoord_s.y,'.')
%     hold on
%     scatter(subcoord.x,subcoord.y,'.')
    plot(morphed(:,1),morphed(:,2))
    hold on
    plot(morphed2(:,1),morphed2(:,2))
end

end