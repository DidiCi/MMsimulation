
%___________________________DATA_EXTRACTION_____________________________
%
%With this program it is possible to extract the data from the excel files 
%obtained from the analysis of combed fibers with ImageJ. 
%The script calculates several global parameters for each experimental time point
%and allow to group time points based on global replicated fraction, then 
%averaged replication parameters are calculated for the selected group.
%In this way, the script allow to check that the difference in replication level 
%beetween different time points is within a certain limit (set with "percentage")
% or to define the thresholds of fluorescence intensity, gaps size, eye
% length size and initiation tracks size. 
% The outputs globalallexDcut.mat, globalalllength_pieces.mat and globalallnum_pieces.mat
%contain the data of all the time points.
% The outputs allexDcut.mat, alllength_pieces.mat and allnum_pieces.mat contain the
% data of the group of time points selected to be analized.
% At the end, the files exresult.mat, file.mat, allexDcut.mat, alllength_pieces.mat and
% allnum_pieces.mat will be used in the folder "Simulation".
%________________________________________________________________________

clear all;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Variables to modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Path to excel data, can contain mutiple path
path={'../Data_demo/data_example/'};
sample_path='Condition1';

%The function "intensities_extraction.m" extract the data from excel files
intensities=intensities_extraction(path);
load('intensities.mat');

%Threshold of fluorescence intensity, one for each data point
thre_int=[20*ones(1,1)];

%Microscope conversion parameters
Conv_microscope=0.065; %1pixel=0.065micrometer
Convmicro_kb=2;

%General variables
unit=1000; %Define the number of base pair (bp) for each block of the genome
v=8.3; %velocity of fork in bp/s

%Variables for extraction
file={'dataexample_1'}; %Name of the files I want to use
% 
pag=1; %Pag from where to start to take the data; correspond to the page in the excel files with the informations on all the fibers

%Gaps smaller then thre1(here in bp)are combined in the analysis of fibers
thre1=1000; 

%Eyes smaller then thre2(here in bp)are not considered in the analysis of
%fiber
thre2=1000;

%Eyes smaller than thre3(here in bp) are considered as new origins.
thre3=3000;
interval=thre3/(2*v); %Detectable initiation events can occur in this interval (in sec) 

%Percentage to set the groups of time points to analyse
percentage=10;
%Number of bins by which the program divide the fibers according to
%the replicated fraction (I didn't use the function 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_binrep=8;
%Number of bins by which the progam divide the ETED, eye length and gap length
%distributions(I didn't use the function 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_bineyes=20;
%Max value by which the progam divide the ETED, eye length and gap length
%distributions
maxlength_bineyes=80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Intensities_treatment%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I apply the threshold for intensities, gaps and eyes and I convert in Kb
[globalallexDcut,globalallnum_pieces,globalalllength_pieces]=intensities_treatment(intensities,file,unit,thre1,thre2,thre_int,Convmicro_kb,Conv_microscope);

%I load the data to analyse, that were previously saved by the function
%'intensities_treatment'
load('globalallexDcut.mat');
load('globalallnum_pieces.mat');
load('globalalllength_pieces.mat');

%%%%%%%%%%%%Calculation of several global parameters for each experimental time point%%%%%%%%%%%%%%%%
for i=1:length(file)   
    
%The fuction calculatestruct.m calculate all the characteristic of the fibers and save them in the structure     
[globalallexDcut.(['exDcut' file{i}])]=calculatestruc(globalallexDcut.(['exDcut' file{i}]),globalallnum_pieces.(['num_pieces' file{i}]),thre1/(globalallexDcut.(['exDcut' file{i}])(1).unit_block),thre2/(globalallexDcut.(['exDcut' file{i}])(1).unit_block));

% I calculate the total replicated fraction, frequency of initiation and fork density for the time points
%"cellfun" apply the given function to each cell in cell array; the
%structure is coverted in cell array with the braket {}.
%The replicated fraction is calculated like the sum of all the 1-blocks of
%all the fibers of a given experimantal time point divided by the sum of
%the lengths of the fibers.
lengthfiber(i)=sum(cellfun(@length,{globalallexDcut.(['exDcut' file{i}]).fiber}));
sumfiber(i)=sum(cellfun(@sum,{globalallexDcut.(['exDcut' file{i}]).fiber}));
exfrac(i)=sumfiber(i)/lengthfiber(i);

%I calculate the new origins on all the fibers of a time point.
%In the vectors containing the lengths of eyes for each fiber, "cellfun" recognize the eye lengths
%less than thre3 (converted in kb) and put a 1 at that positions and 0 at
%the others position. my is a cell array, with a cell for each fiber. 
% 'UniformOutput' in cellfun allow to combine outputs of different size
my=cellfun(@(x) x<=(thre3/(globalallexDcut.(['exDcut' file{i}])(1).unit_block)),{globalallexDcut.(['exDcut' file{i}]).length_eyes},'UniformOutput',0);
neworigins(i)=sum(cellfun(@sum,my));

%I calculate the frequency of initiation (unit=1/kb*sec) and the fork density (1/kb) for each time
%point.
freq_init(i)=neworigins(i)/((lengthfiber(i)-sumfiber(i))*interval);
fork_density(i)=sum(cellfun(@length,{globalallexDcut.(['exDcut' file{i}]).threleft_forkwithside})+cellfun(@length,{globalallexDcut.(['exDcut' file{i}]).threright_forkwithside}))/lengthfiber(i);

end

%For each point I find the points that differ for not more than "percentage" of 
% total fraction of replication

for i=1:length(file)
    groups(i).values=exfrac(exfrac<(exfrac(i)+(percentage/100)) & exfrac>=exfrac(i));
    groups(i).names=file(exfrac<(exfrac(i)+(percentage/100)) & exfrac>=exfrac(i));
    averagereplication(i)=mean(groups(i).values);
end


%############################### PLOTS ###################################
%-------------------------------------------------------------------------
%I plot: the replicated fraction for all the time points; the frequency of
%initiation and the fork density againts the replicated fraction.

figure;
hold all
for i=1:length(file)
scatter(1,exfrac(i),40,'fill');
end
ylabel('Total replicated fraction');
legend(file);
hold off

figure;
hold all
for i=1:length(file)
scatter(exfrac(i),freq_init(i),40,'fill');
end
ylabel('Total frequency of initiation (1/kb*sec)');
xlabel('Replicated fraction');
legend(file);
hold off

figure;
hold all
for i=1:length(file)
scatter(exfrac(i),neworigins(i),40,'fill');
end
ylabel('Total number of origins');
xlabel('Replicated fraction');
legend(file);
hold off

figure;
hold all
for i=1:length(file)
scatter(exfrac(i),fork_density(i),40,'fill');
end
ylabel('Fork density (1/kb)');
xlabel('Replicated fraction');
legend(file);
hold off

%I plot the time points divided in groups that differ for not more than "percentage" of 
% total fraction of replication
figure;
hold all
for i=1:length(file)
scatter(ones(1,length(groups(i).values))*i,groups(i).values,40,'fill');
text(ones(1,length(groups(i).values))*i,groups(i).values,groups(i).names);
end
ylabel('Replicated fraction');
xlabel(['Groups (Rep:' num2str(averagereplication) ')'],'FontSize',8);
hold off
%-------------------------------------------------------------------------
%#########################################################################








%I ask in the command window, what group to analyse
prompt = {'What is the group you want to analyse? '};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'1'};
answer=inputdlg(prompt,dlg_title,num_lines,defaultans);
chosengroup=str2num(answer{1});

%I save in allexDcut.mat, alllength_pieces.mat and allnum_pieces.mat only
%the data of the group that I decided to analyse and I save the new
%structures that I will use to compare with the simulation

allexDcut=globalallexDcut;
for i=find(~(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup)))
    allexDcut=rmfield(allexDcut,['exDcut' file{i}]);
end

allnum_pieces=globalallnum_pieces;
for i=find(~(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup)))
    allnum_pieces=rmfield(allnum_pieces,['num_pieces' file{i}]);
end

alllength_pieces=globalalllength_pieces;
for i=find(~(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup)))
    alllength_pieces=rmfield(alllength_pieces,['length_pieces' file{i}]);
end

file=file(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup));

save([sample_path '/allexDcut.mat'],'allexDcut');
save([sample_path '/allnum_pieces.mat'],'allnum_pieces');
save([sample_path '/alllength_pieces.mat'],'alllength_pieces');
save([sample_path '/file.mat'],'file');

%I save the global replicated fraction, freq. of initiation and fork density only
%for the data of the group that I decided to analyse and I save in exresult.mat
% that I will use to compare with the simulation
exresult.finalresult.tot.fraction_rep=exfrac(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup));
exresult.finalresult.tot.freq_init=freq_init(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup));
exresult.finalresult.tot.fork_density=fork_density(exfrac<(exfrac(chosengroup)+(percentage/100)) & exfrac>=exfrac(chosengroup));
%I save the global replicated fraction, freq. of initiation and fork density 
%of all the data
exresult.finalresult.tot_allpoints.fraction_rep=exfrac;
exresult.finalresult.tot_allpoints.freq_init=freq_init;
exresult.finalresult.tot_allpoints.fork_density=fork_density;

%-------------->Analysis of the points of interest<-------------------

%I calculate separately for each time point frequency of initiation
%and density of forks for each fiber and I group by replicated fraction.
%Then I average the results of the different time points.

for s=1:length(file)
    
exDcut=allexDcut.(['exDcut' file{s}]);
num_pieces=allnum_pieces.(['num_pieces' file{s}]);
length_pieces=alllength_pieces.(['length_pieces' file{s}]);

%Analysis of single time points 
%I divide fiber in accordance with percentage of replication
%The bin are define in a way that edges(k)< x <=edges(k+1)
%The 0 is taken into account in the first bin
%The num_binrep and num_bin eyes parameters should be setted to have a 
%reasonable number of fibers for each bin

centersrep=1/(num_binrep*2):1/(num_binrep):(2*num_binrep-1)/(num_binrep*2);
histrep=zeros(1,length(centersrep));
histforknotave=zeros(1,length(centersrep));
histlengthnotave=zeros(1,length(centersrep));
neworiginsnotave=zeros(1,length(centersrep));
histunrepnotave=zeros(1,length(centersrep));
k=1;
for i=0:1/(num_binrep):(num_binrep-1)/num_binrep
    temp=([exDcut.fracfiber]>i & [exDcut.fracfiber]<=i+(1/(num_binrep)));
    %I calculate the number of fibers in the bin
    histrep(k)=sum(temp);
    %I sum the number of all the forks for the fibers in the bin 
    histforknotave(k)=sum(cellfun(@length,{exDcut(temp).threleft_forkwithside}))+sum(cellfun(@length,{exDcut(temp).threright_forkwithside}));
    %I sum the length of all the fibers in the bin
    histlengthnotave(k)=sum(cellfun(@length,{exDcut(temp).fiber}));
    %I sum the number of all the new origins for the fibers in the bin
    my=cellfun(@(x) x<=(thre3/(exDcut(1).unit_block)),{exDcut(temp).length_eyes},'UniformOutput',0);
    neworiginsnotave(k)=sum(cellfun(@sum,my));
    %I sum the lengths of all unreplicated parts for the fibers in the bin
    histunrepnotave(k)=histlengthnotave(k)-sum(cellfun(@sum,{exDcut(temp).fiber}));
    k=k+1;
end
%I add the quantities for the 0
temp=([exDcut.fracfiber]==0);
histrep(1)=histrep(1)+sum(temp);
histforknotave(1)=histforknotave(1)+sum(cellfun(@length,{exDcut(temp).threleft_forkwithside}))+sum(cellfun(@length,{exDcut(temp).threright_forkwithside}));
histlengthnotave(1)=histlengthnotave(1)+sum(cellfun(@length,{exDcut(temp).fiber}));
my=cellfun(@(x) x<=(thre3/(exDcut(1).unit_block)),{exDcut(temp).length_eyes},'UniformOutput',0);
neworiginsnotave(1)=neworiginsnotave(1)+sum(cellfun(@sum,my));
histunrepnotave(1)=histunrepnotave(1)+histlengthnotave(1)-sum(cellfun(@sum,{exDcut(temp).fiber}));
%I calculate fork density and frequency of initiation and I put to 0 the NaN
%values; a point is NaN if there are not origins or forks for that bin.
histforkdens_notave=histforknotave./histlengthnotave;
histfreqinit_notave=neworiginsnotave./(histunrepnotave.*interval);
histforkdens_notave(isnan(histforkdens_notave))=0;
histfreqinit_notave(isnan(histfreqinit_notave))=0;
%I save the results in a structure
exresult.(file{s}).fraction_rep.centers=centersrep;
exresult.(file{s}).fraction_rep.histrep=histrep/sum(histrep);
exresult.(file{s}).fork_density.histforkdens_notave=histforkdens_notave;
exresult.(file{s}).freq_init.histfreqinit_notave=histfreqinit_notave;
end

%I calculate separately for each time point ETED, eye length and gap length distribution.
%Then I average the result of the different time points.
for s=1:length(file)
    
exDcut=allexDcut.(['exDcut' file{s}]);
num_pieces=allnum_pieces.(['num_pieces' file{s}]);
length_pieces=alllength_pieces.(['length_pieces' file{s}]);

%Analysis of single time points 
%I divide fiber in accordance with a length range.
%The bin are define in a way that edges(k)< x <=edges(k+1)
%The 0 is taken into account in the first bin
%In this case I consider lentghs beetween 1-100kb.

centerseyes=maxlength_bineyes/(num_bineyes*2):maxlength_bineyes/(num_bineyes):maxlength_bineyes*(2*num_bineyes-1)/(num_bineyes*2);

histeted=zeros(1,length(centerseyes));
histeyes=zeros(1,length(centerseyes));
histgaps=zeros(1,length(centerseyes));

k=1;
%I extract the ETED, eye lengths and gap lengths from all the fibers of the
%given time point.
%There will be no 0 values (because are all lengths) so I don't need to add 
%the 0 values in the first bin like for the replicated fraction
my={exDcut.etedist};
eted=cell2mat(my(~cellfun(@isempty,my))');
my={exDcut.length_eyes};
eyes=cell2mat(my(~cellfun(@isempty,my))');
my={exDcut.gap_length};
gaps=cell2mat(my(~cellfun(@isempty,my))');
for i=0:maxlength_bineyes/(num_bineyes):maxlength_bineyes*(num_bineyes-1)/num_bineyes
    %I calculate the number of ETED in the bin
    histeted(k)=sum(eted>i & eted<=i+(maxlength_bineyes/(num_bineyes)));
     %I calculate the number of eye lengths in the bin
    histeyes(k)=sum(eyes>i & eyes<=i+(maxlength_bineyes/(num_bineyes)));
     %I calculate the number of gap lengths in the bin
    histgaps(k)=sum(gaps>i & gaps<=i+(maxlength_bineyes/(num_bineyes)));
    k=k+1;
end
%I normalize the distributions and I put to 0 the NaN
%values; a point is NaN if there are not origins or forks for that bin.
histeted=histeted./sum(histeted);
histeyes=histeyes./sum(histeyes);
histgaps=histgaps./sum(histgaps);
histeted(isnan(histeted))=0;
histeyes(isnan(histeyes))=0;
histgaps(isnan(histgaps))=0;
%I save the results in a structure
exresult.(file{s}).eted.centers=centerseyes;
exresult.(file{s}).eted.histeted=histeted;
exresult.(file{s}).length_eyes.histeyes=histeyes;
exresult.(file{s}).length_gaps.histgaps=histgaps;
end

%I save the total replicated fraction for the points
for g=1:length(file)
exresult.(file{g}).fraction_rep.totfraction_rep=groups(chosengroup).values(g);
end

%I find min and max for each fraction of replication
for s=1:length(exresult.(file{1}).fraction_rep.centers)
%I put the values of replication for all the files in a single array
for i=1:length(file)
    histrepcenter(i)=exresult.(file{i}).fraction_rep.histrep(s);
end
%I find the min and max of the values
histrepmax(s)=max(histrepcenter);
histrepmin(s)=min(histrepcenter);
%I put the values of fork density for all the files in a single array
for i=1:length(file)
    histfreqinit_notavecenter(i)=exresult.(file{i}).freq_init.histfreqinit_notave(s);
end
%I find the min and max of the values
histfreqinit_notavemax(s)=max(histfreqinit_notavecenter);
histfreqinit_notavemin(s)=min(histfreqinit_notavecenter);
%I put the values of initiation frequency for all the files in a single array
for i=1:length(file)
    histforkdens_notavecenter(i)=exresult.(file{i}).fork_density.histforkdens_notave(s);
end
%I find the min and max of the values
histforkdens_notavemax(s)=max(histforkdens_notavecenter);
histforkdens_notavemin(s)=min(histforkdens_notavecenter);

end
%I do the same for ETED, eye length and gap length distributions
for s=1:length(exresult.(file{1}).eted.centers)
%I put the values of eted distribution, for all the files in a single array
for i=1:length(file)
    histetedcenter(i)=exresult.(file{i}).eted.histeted(s);
end
%I find the min and max beetween the values
histetedmax(s)=max(histetedcenter);
histetedmin(s)=min(histetedcenter);
%I put the values of eye lengths distribution for all the files in a single array
for i=1:length(file)
    histeyescenter(i)=exresult.(file{i}).length_eyes.histeyes(s);
end
%I find the min and max of the values
histeyesmax(s)=max(histeyescenter);
histeyesmin(s)=min(histeyescenter);
%I put the values of gaps lengths distribution for all the files in a single array
for i=1:length(file)
    histgapscenter(i)=exresult.(file{i}).length_gaps.histgaps(s);
end
%I find the min and max of the values
histgapsmax(s)=max(histgapscenter);
histgapsmin(s)=min(histgapscenter);

end
%I save all the results in the structure
exresult.finalresult.fraction_rep.mean=(histrepmax+histrepmin)/2;
exresult.finalresult.fraction_rep.error=(histrepmax-histrepmin)/2;
exresult.finalresult.fork_density.mean=(histforkdens_notavemax+histforkdens_notavemin)/2;
exresult.finalresult.fork_density.error=(histforkdens_notavemax-histforkdens_notavemin)/2;
exresult.finalresult.freq_init.mean=(histfreqinit_notavemax+histfreqinit_notavemin)/2;
exresult.finalresult.freq_init.error=(histfreqinit_notavemax-histfreqinit_notavemin)/2;
exresult.finalresult.eted.mean=(histetedmax+histetedmin)/2;
exresult.finalresult.eted.error=(histetedmax-histetedmin)/2;
exresult.finalresult.length_eyes.mean=(histeyesmax+histeyesmin)/2;
exresult.finalresult.length_eyes.error=(histeyesmax-histeyesmin)/2;
exresult.finalresult.length_gaps.mean=(histgapsmax+histgapsmin)/2;
exresult.finalresult.length_gaps.error=(histgapsmax-histgapsmin)/2;
exresult.finalresult.num_binrep=num_binrep;
exresult.finalresult.num_bineyes=num_bineyes;
exresult.finalresult.maxlength_bineyes=maxlength_bineyes;


%############################### PLOTS ###################################
%-------------------------------------------------------------------------
%I plot for replicated fraction, fork density and frequency of initiation
%the plots of all the points togheter and a figure with the plot with
%average and errors bar

figure;
hold all
xlabel('Fraction of replication');
ylabel('Normalized number of fibers');
for s=1:length(file)
plot(centersrep,exresult.(file{s}).fraction_rep.histrep,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('Fraction of replication ');
ylabel('Forks density (1/kb)');
for s=1:length(file)
plot(centersrep,exresult.(file{s}).fork_density.histforkdens_notave,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('Fraction of replication');
ylabel('I(f)(1/kb*sec)');
for s=1:length(file)
plot(centersrep,exresult.(file{s}).freq_init.histfreqinit_notave,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('Fraction of replication','FontSize',14);
ylabel('Normalized number of fibers','FontSize',12);
errorbar(centersrep,exresult.finalresult.fraction_rep.mean,exresult.finalresult.fraction_rep.error,'r');
legend('Experiment');
hold off

figure;
hold all
xlabel('Fraction of replication','FontSize',14);
ylabel('Average forks density (1/kb)','FontSize',12);
errorbar(centersrep,exresult.finalresult.fork_density.mean,exresult.finalresult.fork_density.error,'r');
legend('Experiment');
hold off

figure;
hold all
xlabel('Fraction of replication','FontSize',14);
ylabel('Average I(f)(1/kb*sec)','FontSize',14);
errorbar(centersrep,exresult.finalresult.freq_init.mean,exresult.finalresult.freq_init.error,'r');
legend('Experiment');
hold off

figure;
hold all
xlabel('ETED (kb)','FontSize',14);
ylabel('Normalized ETED distribution','FontSize',12);
for s=1:length(file)
plot(centerseyes,exresult.(file{s}).eted.histeted,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('Eyes length (kb) ','FontSize',14);
ylabel('Normalized eyes length distribution','FontSize',12);
for s=1:length(file)
plot(centerseyes,exresult.(file{s}).length_eyes.histeyes,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('Gaps length (kb) ','FontSize',14);
ylabel('Normalized gaps length distribution','FontSize',14);
for s=1:length(file)
plot(centerseyes,exresult.(file{s}).length_gaps.histgaps,'--s');
end
legend(file);
hold off

figure;
hold all
xlabel('ETED (kb)','FontSize',14);
ylabel('Normalized ETED distribution','FontSize',12);
errorbar(centerseyes,exresult.finalresult.eted.mean,exresult.finalresult.eted.error,'r');
legend('Experiment');
hold off

figure;
hold all
xlabel('Eyes length (kb) ','FontSize',14);
ylabel('Normalized eyes length distribution','FontSize',12);
errorbar(centerseyes,exresult.finalresult.length_eyes.mean,exresult.finalresult.length_eyes.error,'r');
legend('Experiment');
hold off

figure;
hold all
xlabel('Gaps length (kb) ','FontSize',14);
ylabel('Normalized gaps length distribution','FontSize',14);
errorbar(centerseyes,exresult.finalresult.length_gaps.mean,exresult.finalresult.length_gaps.error,'r');
legend('Experiment');
hold off
%--------------------------------------------------------------------------
%##########################################################################

save([sample_path '/exresult.mat'],'exresult');
