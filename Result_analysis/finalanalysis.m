%In this code I analyse the results of the genetic algorithm
%I take the best of each garesult for the given model and I calculate
%quantities that will be used to produce the images
clear all;
close all;

%I load the experimental data already treated in 'Data extraction and global analysis'
load('../Data_extraction/Condition1/allexDcut.mat');
load('../Data_extraction/Condition1/allnum_pieces.mat');
load('../Data_extraction/Condition1/alllength_pieces.mat');
load('../Data_extraction/Condition1/file.mat');
load('../Data_extraction/Condition1/exresult.mat');

%Number of simulation in "Simulation" folder
n_simu=100;
path_in=['../Simulation/Condition1/'];
path_out=['Condition1/']; 

for NUM_garesult=1:3
    
tic,
sprintf('Treatment of file n. %i',NUM_garesult)
garesult=[];
load([path_in 'garesult' num2str(NUM_garesult) '.mat']); 

%General variables
unit=allexDcut.(['exDcut' file{1}])(1).unit_block; %Define the number of base pair (bp) for each block of the genome
fraction_max=1; %Fraction of replication at which I block the replication
num_DNA=100; %Number of DNA I simulate 
time_max=100; %Time at which I block the replication 
v=8.3; %velocity of fork in bp/s
length_DNA=100000000/unit; %DNA lenght
timeunit_insec=100000000/(length_DNA*v); %I know velocity and that I replicate a block/step
origin_density=1/2.3 ; %kb^(-1)

%Gaps smaller then thre1(here in bp)are combined in the analysis of fibers
thre1=1000; 

%Eyes smaller then thre2(here in bp)are not considered in the anamysis of
%fiber
thre2=1000;

%Eyes smaller than thre3(here in bp) are considered as new origins.
thre3=3000;
interval=thre3/(2*v); %Detectable initiation events can occur in this interval (in sec)

%I analyse the result with the best score for each simulation
var_simu=garesult.x;
namedna=[path_out 'Dna' num2str(NUM_garesult) '.mat'];
Dna=simulationdna_foranalysis(var_simu,namedna,unit,fraction_max,num_DNA,time_max,length_DNA,timeunit_insec);

%I round the variables that are integers in the simulation
var_simu(1)=round(var_simu(1));
var_simu(2)=round(var_simu(2));
var_simu(6)=round(var_simu(6));


%--------------->Analysis of the time of interest<-----------------------
%#########################################################################
%#########################################################################

[thresult,cost,costpart,reference_thtime]=calculatethresult(Dna,file,exresult,allnum_pieces,alllength_pieces,num_DNA,length_DNA,unit,v,1);
thresult.finalresult.reference_thtime=reference_thtime;

%#########################################################################
%-------------------------------------------------------------------------

save([path_out 'thresult' num2str(NUM_garesult) '.mat'],'thresult');

toc
end

