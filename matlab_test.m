fid = fopen('./XFOIL/xfoilinput.txt', 'w');
fprintf(fid,'%haha 15.12f\n', 3.1415926);
fclose(fid);

system('echo ls')