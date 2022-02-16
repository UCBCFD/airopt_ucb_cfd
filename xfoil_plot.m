function [] = xfoil_plot(input)
% input = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

% Generate result and update morphed_repanel.txt
input = num2cell(input);
tstart = tic;
[CLD_max,alpha_stall,error_flag]=airfoil_func_calls(input{:});
tend = toc(tstart);
fprintf('%g %g %g %g \n', [CLD_max,alpha_stall,tend,error_flag]);

% Plot xfoil results
alpha = 0;
alpha_end = 40;
alpha_step = 0.25;
[CLD,CL,alpha_list,~] = xfoil_loop(alpha,alpha_end,alpha_step,3,[]);
figure
subplot(2,1,1);
scatter(alpha_list,CL,'.')
title('CL vs alpha')
subplot(2,1,2);
scatter(alpha_list,CLD,'.')
title('CLD vs alpha')

end