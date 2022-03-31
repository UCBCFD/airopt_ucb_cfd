n_processors=8;
% SO Optimums
gen_intial=zeros(2,25);
gen_intial_opt=zeros(2,1);
lb=-1*ones(1,25);
ub=ones(1,25);
parpool(n_processors,'IdleTimeout', Inf)
% options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false);
options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false, 'Display','iter','OutputFcn',@gaoutfun);
[gen_intial(1,:),gen_intial_opt(1)]= ga(@(x) airfoil_wrapper_CLD(x),25,[],[],[],[],lb,ub,[],options);
poolobj = gcp('nocreate');
delete(poolobj);

parpool(n_processors,'IdleTimeout', Inf)
% options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false);
options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false, 'Display','iter','OutputFcn',@gaoutfun);
[gen_intial(2,:),gen_intial_opt(2)]= ga(@(x) airfoil_wrapper_AS(x),25,[],[],[],[],lb,ub,[],options);
poolobj = gcp('nocreate');
delete(poolobj);

gen_intial=vertcat(gen_intial,eye(25));

parpool(n_processors,'IdleTimeout', Inf)
options = optimoptions('gamultiobj','DistanceMeasureFcn',{@distancecrowding,'phenotype'},'ConstraintTolerance',1e-4,'MaxGenerations',5000,'PopulationSize',376,'UseParallel',true, 'UseVectorized', false, 'InitialPopulationMatrix', gen_intial,'Display','iter','OutputFcn',@gaoutfun, 'FunctionTolerance',10^-8);
[x,fval,exitflag,output,population,scores] = gamultiobj(@(x) airfoil_wrapper(x),25,[],[],[],[],lb,ub,[],options);
poolobj = gcp('nocreate');
delete(poolobj);

save('GAOutput')