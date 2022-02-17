input = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
% fid = fopen('test.txt','r');
for i = 1:numel(input)
    input{i} = 1;
    tstart = tic;
    [CLD_max,alpha_stall,error_flag]=airfoil_func_calls(input{:});
    tend = toc(tstart);
    input{i} = 0;
    fprintf('%g %g %g %g %g \n', [i,CLD_max,alpha_stall,tend,error_flag]);
end
% fclose(fid);