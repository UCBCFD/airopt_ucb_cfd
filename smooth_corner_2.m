function [morphed_total]=smooth_corner_2(filename,angle_threshold)

% morphed_total = importdata(filename,' ',1);
% morphed_total = morphed_total.data;
morphed_total = importdata(filename);

morphed=morphed_total(10:end-10,:);
coord = morphed(:,1)+morphed(:,2).*1i;
morphed2 = morphed;
vec = diff(coord);
ang = wrapToPi(diff(angle(vec))+pi); %angle formed by neighboring 3 pts
ang2 = abs(rad2deg(ang));
edge = find(ang2<angle_threshold)+1; %collocation pt is offset by 1

for i=1:size(edge)
    lowbnd = max(edge(i)-20,1);
    upbnd = min(edge(i)+20,size(coord,1));
   subcoord = morphed(lowbnd:upbnd,:);
   subcoord_s = smoothdata(subcoord,'sgolay',upbnd-lowbnd+1);
   morphed2(lowbnd+2:upbnd-2,:) = subcoord_s(3:end-2,:);
end

% scatter(subcoord_s.x,subcoord_s.y,'.')
% hold on
% scatter(subcoord.x,subcoord.y,'.')
% plot(morphed_total(:,1),morphed_total(:,2))
% hold on
morphed_total(10:end-10,:)=morphed2;
% plot(morphed_total(:,1),morphed_total(:,2),'b')
% hold on
morphed_total(:,1)=smooth(morphed_total(:,1),7);
morphed_total(:,2)=smooth(morphed_total(:,2),7);
% plot(morphed_total(:,1),morphed_total(:,2),'r')

fid = fopen('morphed_smooth.txt','w');
for i = 1:size(morphed_total,1)
fprintf(fid,'%g %g \n',morphed_total(i,1),morphed_total(i,2));
end
fclose(fid);

end