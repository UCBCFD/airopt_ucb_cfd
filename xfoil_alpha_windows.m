function [output,alpha_return] = xfoil_alpha_windows(alpha, A, flag, alpha_end, alpha_step)

% UPDATE XFOIL'S INPUT
if ~exist('alpha_end','var')
    A{13} = sprintf('alfa %g',alpha);
else
    A{13} = sprintf('aseq %g %g %g',[alpha,alpha_end,alpha_step]);
end
fid = fopen('./XFOIL/xfoilInput.input', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break 
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);

% RUN XFOIL
system('del /f XFOIL\xfoilSave.txt > NUL');
system('del /f XFOIL\xfoilDump.txt > NUL');

% system('xfoil.exe < XFOIL\xfoilInput.input')
results=evalc("system('xfoil.exe < XFOIL\xfoilInput.input > NUL')");

% OPEN XFOIL OUTPUT
fid = fopen('./XFOIL/xfoilSave.txt');
if fid == -1
    output = 0;
    return
end

% READ XFOIL OUTPUT
header = textscan(fid,'%s',7,'HeaderLines',10);
header = char(header{1,1});
frewind(fid)
data_raw = textscan(fid,'%f %f %f %f %f %f %f','HeaderLines',12);
data_raw = cell2mat(data_raw);
fclose(fid);

% CREATE DATA STRUCTURE
for i = 1:7
    name = strtrim(convertCharsToStrings(header(i,:)));
	data.(name) = data_raw(:,i);
end

% max(CL/CD)
if isempty(data.CL) || isempty(data.alpha)
    CLD = 0;
    CL = 0;
%     alpha_return = 0;
    alpha_return = [];
else
    CLD = data.CL./data.CD;
    CL = data.CL;
    alpha_return = data.alpha;
end

if strcmpi(flag,'CLD')
    output = CLD;
elseif strcmpi(flag,'CL')
    output = CL;
elseif strcmpi(flag,'all')
    output = [CLD(:),CL(:)];
else
    output = nan;
end


end