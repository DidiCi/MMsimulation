clear all;
close all;
addpath('../Functions') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of variables
nvars=7;
%Lower bounds
LB=[1 0 0 0 0 0 0];
%Upper bounds
UB=[200 400 1 1 1 500 1];
%Number of simulations to do
n_round=3;
%Path of output
path_out='Condition1/';

%I load the experimental data obtained in 'Data extraction'
%Indicate the right path
load('../Data_extraction/Condition1/allexDcut.mat');
load('../Data_extraction/Condition1/allnum_pieces.mat');
load('../Data_extraction/Condition1/alllength_pieces.mat');
load('../Data_extraction/Condition1/file.mat');
load('../Data_extraction/Condition1/exresult.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Custom options for the genetic algorithm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Check the manual for more information%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matlabpool open local 4 (for parallel computing if using older versions)

%Positions of integer variables
% IntCon=[1 2 4 8];

%Set options
%ga ignores or overwrites some options for mixed integer optimization problems
options = gaoptimset(@ga);

% Handle to the function that creates the initial population
% @gacreationuniform | @gacreationlinearfeasible
% options.CreationFcn=;

% Handle to the function that the algorithm uses to create crossover children
% @crossoverheuristic | {@crossoverscattered} | @crossoverintermediate | @crossoversinglepoint | @crossovertwopoint | @crossoverarithmetic
% options.CrossoverFcn=

%The fraction of the population at the next generation, not including elite
%children, that is created by the crossover function (default=0.8)
options.CrossoverFraction=0.6;

%Positive integer specifying how many individuals in the current generation
%are guaranteed to survive to the next generation ({0.05*(default population size)} for mixed integer problems)
options.EliteCount=3;

%Scalar. If the fitness function attains the value of FitnessLimit, the
%algorithm halts. (default=-Inf)
options.FitnessLimit=0;

%Hibrid faction
% hybridopts = optimset('Display','iter','MaxFunEvals',50);
% options.HybridFcn={@fmincon,hybridopts};

%Handle to the function that scales the values of the fitness function
% @fitscalingshiftlinear | @fitscalingprop | @fitscalingtop | default:@fitscalingrank
% options.FitnessScalingFcn=@fitscalingprop;

%Positive integer specifying the maximum number of iterations before the
%algorithm halts (default=100)
options.Generations=50;


%Size of the population (default={min(max(10*nvars,40),100)} for mixed integer problems)
options.PopulationSize=[20 20 20];

% Initial population used to seed the genetic algorithm; can be partial
% for i=1:nvars
% population(:,i)=random('unif',LB(i),UB(i),options.PopulationSize,1);
% end
% options.InitialPopulation=population;


%Direction of migration
% 'both' | defalut:'forward'
% options.MigrationDirection;

% Scalar between 0 and 1 specifying the fraction of individuals in each subpopulation that migrates to a different subpopulation
% default=0.2
options.MigrationFraction=0.1;

% Positive integer specifying the number of generations that take place between migrations of individuals between subpopulations
% default=20
options.MigrationInterval=5;

% Handle to the function that produces mutation children
% {@mutationuniform,rate] | @mutationadaptfeasible | {@mutationgaussian}
options.MutationFcn={@mutationuniform,0.2};

%Array of handles to functions that plot data computed by the algorithm
%@gaplotbestf | @gaplotbestindiv | @gaplotdistance | @gaplotexpectation | @gaplotgeneology
% | @gaplotmaxconstr | @gaplotrange | @gaplotselection | @gaplotscorediversity | @gaplotscores | @gaplotstopping 
options.PlotFcns={@gaplotbestf, @gaplotbestindiv, @gaplotscores, @gaplotdistance};

%Positive integer specifying the number of generations between consecutive
%calls to the plot functions (default=1)
% options.PlotInterval;

% Matrix or vector specifying the range of the individuals in the initial population 
options.PopInitRange=[LB;UB];


% Handle to the function that selects parents of crossover and mutation children
% @selectionremainder | @selectionuniform | {@selectionstochunif} |
% @selectionroulette | @selectiontournament 
options.SelectionFcn=@selectionroulette;

%Positive integer. The algorithm stops if there is no improvement in the objective function for StallGenLimit consecutive generations.
% default=50
options.StallGenLimit=200;

% Positive scalar. The algorithm stops if there is no improvement in the objective function for StallTimeLimit seconds.
% {default=Inf}
% options.StallTimeLimit;

% Positive scalar. The algorithm stops after running for TimeLimit seconds.
% {default:Inf}
% options.TimeLimit;

% Positive scalar. TolCon is used to determine the feasibility with respect to nonlinear constraints.
% Positive scalar | {1e-6}
% options.TolCon;

% Positive scalar. The algorithm runs until the cumulative change in the fitness function value over StallGenLimit is less than TolFun.
% {default:1e-6}
 options.TolFun=2;
 
% Compute fitness functions of a population in parallel.
% 'always' | {default='never'}
 options.UseParallel='always';


% String specifying whether the computation of the fitness function is
% vectorized (not possible in my case)
% 'on' | {default='off'} 
% options.Vectorized;


A=[0 0 -1 1 0 0 0];
B=[0];


%Set the random number generetor for reproducibility
rng('shuffle');
for z=1:n_round
tic,
clear garesult
[garesult.x,garesult.fval,garesult.exitflag,garesult.output,garesult.population,garesult.scores] = ga(@(var)simulation_costfunction(var,allexDcut,allnum_pieces,alllength_pieces,file,exresult),nvars,A,B,[],[],LB,UB,[],options);
name=[path_out 'garesult' num2str(z) '.mat'];
save(name,'garesult');
saveas(gcf,[path_out 'ga' num2str(z) '.fig'],'fig');
toc
end

%matlabpool close (for parallel computing if using older versions)
