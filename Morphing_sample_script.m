%% MATLAB morphing sample code run
% S Lee 10/19/21
% Airfoil profiles recorded in 1-25.txt
% As for the model specs, see the google spreadsheet shared by Jinge
% [M, T] = Morphing(val1, val2, ...)
% M: Coord. (x,y). Follows the Selig data convention (CCW from the trail)
% T: Flag param. Indicates whether or not intersect. happened more than 5
% val1, val2: weight params. (-1~1) for AF1(1.txt), AF2(2.txt), ...
clc; clear;

[M, T] = Morphing(  0.3,  -1.0,  0.3,  0.3,  0.2, ...
                    0. ,   0.2 ,  0. ,  0. ,  0. , ...
                    0. ,   1. ,  0. ,  0.1 ,  0.3 , ...
                    0.4 ,   1. ,  0.2 ,  0. ,  0. , ...
                    0. ,   0. ,  1. ,  0.4 ,  0. ); % Change weights here

figure(1);
plot(M(:,1),M(:,2),'-');
xlim([0 1]); ylim([-0.2 0.2]);
title('Morphing based on Haris` code');

% % Compare with uncorrected morphing below ... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampleweight = [0.3 -1.0 0.3 0.3 0.2];                                    %
% w_tot = 0;                                                                %
% for i=1:length(sampleweight)                                              %
%    if i==1                                                                %
%        M2=sampleweight(i)*readmatrix(sprintf('./BaseShapes/%d.txt', i));               %
%    else                                                                   %
%        M2=M2+sampleweight(i)*readmatrix(sprintf('./BaseShapes/%d.txt', i));            %
%    end                                                                    %
%    w_tot=w_tot+sampleweight(i);                                           %
% end                                                                       %
% M2=M2./w_tot;                                                             %
%                                                                           %
% figure(2);                                                                %
% plot(M2(:,1),M2(:,2),'-');                                                %
% xlim([0 1]); ylim([-0.2 0.2]);                                            %
% title('Morphing, no intersection correction');                            %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save morphed shapes:
% 1. Down-sample
if numel(M(:,1))>495 %495 is the max num of panel nodes allowed in XFOIL 
    M_undersampled = M(1:ceil(numel(M(:,1))/495):end,:);
else
    M_undersampled = M;
end
% 2. Save morphed shape
fid = fopen('./XFOIL/morphed.txt', 'w');
for i = 1:numel(M_undersampled(:,1))
    fprintf(fid,'%15.12f %15.12f\n', M_undersampled(i,:));
end
fclose(fid);

% If you run XFOIL with this morph...
if (T == true) % If morphing succeeds..

%     [CLD_max,alpha_stall,alpha_CLD_max] = xfoil_read('load XFOIL/morphed.txt');
    [CLD_max,alpha_stall] = xfoil_fmax('load XFOIL/morphed.txt');

else % If morphing fails...

    % Set low values
    CLD_max = 0;
    alpha_stall = 0;

end