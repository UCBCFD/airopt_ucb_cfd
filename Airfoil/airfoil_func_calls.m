function [CLD_max,alpha_stall]=airfoil_func_calls(varargin)

error_flag = -1;
[M, T] = Morphing(varargin{:}); % Change weights here

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

% Smoothing
angle_threshold = 30;
smooth_level = 1;
[M_smooth]=Smoothing(M,angle_threshold,smooth_level);

% If smoothing fails..
if (sum(isnan(M_smooth),'all') > 0)
    % Set low values
    CLD_max = 0;
    alpha_stall = 0;
    error_flag = 66;
    return
end

% Paneling
panel_default = 200;
panel_level = 0;
panel_step = 50;
counter_break=0;
while error_flag < 0 && panel_level < 2 && counter_break<5
    counter_break=counter_break+1;
    panel_current = panel_default + panel_level*panel_step;
    [M]=Panelling(M_smooth,panel_current);
    
    % Saving
    % 1. Down-sample % Since Paneling always gives less than 400 pts, this
    % portion is currently dummy.
    if numel(M(:,1))>495 %495 is the max num of panel nodes allowed in XFOIL
        M_undersampled = M(1:ceil(numel(M(:,1))/495):end,:);
    else
        M_undersampled = M;
    end
    % 2. Normalize (i.e. scale the chord length to 1.00)
    le = find(M_undersampled(:,1)==min(M_undersampled(:,1))); le = le(1);
    M_undersampled(:,1) = M_undersampled(:,1) - M_undersampled(le,1);
    M_undersampled=M_undersampled / M_undersampled(1,1);
    % 3. Save morphed shape
%     system('del /f XFOIL\morphed_repanel.txt > NUL');
    system('rm ./XFOIL/morphed_repanel.txt &> /dev/null');
    fid = fopen('./XFOIL/morphed_repanel.txt', 'w');
    for i = 1:numel(M_undersampled(:,1))
        fprintf(fid,'%15.12f %15.12f\n', M_undersampled(i,:));
    end
    fclose(fid);
    
    % Find CLD_max & alpha_stall
    % If you run XFOIL with this morph...
    if (T == true) % If morphing succeeds..
        [CLD_max,CLD_alpha_max,~,alpha_stall,error_flag] = xfoil_scan();
        alpha_stall=alpha_stall-CLD_alpha_max;
        if error_flag > 0
            CLD_max = 0;
            alpha_stall = 0;
        end
    else % If morphing fails...
        % Set low values
        CLD_max = 0;
        alpha_stall = 0;
        error_flag = 66;
    end
    
    % Increase panels if necessary
%     if error_flag == -321 || error_flag == -322
    if error_flag < 0
        panel_level = panel_level + 1;
    end
end
if panel_level >= 2 % Can't refine the panel towards convergence
    CLD_max = 0; 
    % The airfoil itself might still be a good one though. Not sure how to
    % deal with this right now
    alpha_stall = 0;
    error_flag = 77;
end
if counter_break==5
  disp(counter_break)
end
% ERROR FLAG:
% 66: failed morphing
% 77: maximum panel level reached and still no convergence (note: when too
%     many panel pts, XFOIL can also lead to numerical artifact)
% 3XX: error flag from xfoil_scan
