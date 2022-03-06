n_processors=8;
% SO Optimums
gen_intial=readmatrix('CurrentPopulation.txt');
lb=-1*ones(1,25);
ub=ones(1,25);

parpool(n_processors,'IdleTimeout', Inf)
options = optimoptions('gamultiobj','DistanceMeasureFcn',...
    {@distancecrowding,'phenotype'},'ConstraintTolerance',1e-4,...
    'MaxGenerations',5000,'PopulationSize',372,'UseParallel',true,...
    'UseVectorized', false, 'InitialPopulationMatrix', gen_intial,...
    'Display','iter','OutputFcn',@gaoutfun, 'FunctionTolerance',10^-8);
[x,fval,exitflag,output,population,scores] = gamultiobj(@(x) airfoil_wrapper(x),25,[],[],[],[],lb,ub,[],options);
poolobj = gcp('nocreate');
delete(poolobj);

save('GAOutput')
