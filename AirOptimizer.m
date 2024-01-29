n_processors=124;
maxNumCompThreads(1);

lb=-1*ones(1,25);
ub=ones(1,25);

parpool(n_processors,'IdleTimeout', Inf)

if not(isfile('GAOutput_SO.mat'))
	gen_intial=zeros(2,25);
	gen_intial_opt=zeros(2,1);	

	% CLD MAX SINGLE OBJECTIVE OPTIMIZATION
	options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false, 'Display','iter','OutputFcn',@gaoutfun);
	[gen_intial(1,:),gen_intial_opt(1)]= ga(@(x) airfoil_wrapper_CLD(x),25,[],[],[],[],lb,ub,[],options);

	% DELTA ALPHA SINGLE OBJECTIVE OPTIMIZATION
	options = optimoptions('ga','ConstraintTolerance',1e-6,'MaxGenerations',100,'PopulationSize',128,'UseParallel',true, 'UseVectorized', false, 'Display','iter','OutputFcn',@gaoutfun);
	[gen_intial(2,:),gen_intial_opt(2)]= ga(@(x) airfoil_wrapper_AS(x),25,[],[],[],[],lb,ub,[],options);

	gen_intial=vertcat(gen_intial,eye(25));

	save('GAOutput_SO')
else

	if not(isfile('GAOutput_MO.mat'))
		SO=load('GAOutput_SO.mat');
		gen_intial=SO.gen_intial;
		gen_no=0;
	else
		MO=load('GAOutput_MO.mat');
		gen_intial=MO.population;
		gen_no=MO.gen_no;
	end

	options = optimoptions('gamultiobj','DistanceMeasureFcn',{@distancecrowding,'phenotype'},'ConstraintTolerance',1e-4,'MaxGenerations',500,'PopulationSize',376,'UseParallel',true, 'UseVectorized', false, 'InitialPopulationMatrix', gen_intial,'Display','iter','OutputFcn',@gaoutfun, 'FunctionTolerance',1e-8, 'MaxTime', 172000);
	[x,fval,exitflag,output,population,scores] = gamultiobj(@(x) airfoil_wrapper(x),25,[],[],[],[],lb,ub,[],options);
	gen_no=gen_no+output.generations;

	save('GAOutput_MO')
end

poolobj = gcp('nocreate');
delete(poolobj);
