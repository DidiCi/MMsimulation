function [thresult,cost,costpart,reference_thtime]=calculatethresult(Dna,file,exresult,allnum_pieces,alllength_pieces,num_DNA,length_DNA,unit,v,flag_totSphase)

%Put flag_totSphase==1 to calculate all S phase quantities

%Gaps smaller then thre1(here in bp)are combined in the analysis of fibers
thre1=1000; 

%Eyes smaller then thre2(here in bp)are not considered in the anamysis of
%fiber
thre2=1000;

%Eyes smaller than thre3(here in bp) are considered as new origins.
thre3=3000;
interval=thre3/(2*v); %Detectable initiation events can occur in this interval (in sec) 

%Number of bins by which the progam divide the fibers according to
%the replicated fraction (I didn't use the fuction 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_binrep=exresult.finalresult.num_binrep;

%Number of bins by which the progam divide the ETED, eye length and gap length
%distributions(I didn't use the fuction 'hist' because there were
%inconsistencies in the treatment of the fibers)
num_bineyes=exresult.finalresult.num_bineyes;

%Max value by which the progam divide the ETED, eye length and gap length
%distributions
maxlength_bineyes=exresult.finalresult.maxlength_bineyes;

%--------------->Cut and calculation cost function<-----------------------
%#########################################################################
%#########################################################################

%I calculate the minimum and maximum time time at which the genomes have finished the
%simulation
mintime=min(cellfun(@length,{Dna.time}));
finetime=max(cellfun(@length,{Dna.time}));
timeline=1:finetime;

%I calculate the average of the replicated fraction of all the genome at
%each time
for i=1:mintime
average_thtotfracrep(i)=mean(cellfun(@(x) x(i).fraction_rep,{Dna.time}));
end

%Average fraction of replication of reference from experiment
reference_exfrac=mean(exresult.finalresult.tot.fraction_rep);

%I find the time of the simulated S phase closest to the value of
%replication of experiment
%At the end of the analysis 'reference_time' is the time at which the
%simulated genomes are cutted
reference_thtime=find(abs(average_thtotfracrep-reference_exfrac)==min(abs(average_thtotfracrep-reference_exfrac)), 1, 'last' );
%I reproduce the combing experiment on the genome
%I cut the simulated genomes with the same distribution of the experimet
%with less fibers.
for h=1:length(file)
    num_piecesperfile(h,1)=allnum_pieces.(['num_pieces' file{h}]);
end
reference_exdistr=find(num_piecesperfile==min(num_piecesperfile));
num_pieces=allnum_pieces.(['num_pieces' file{reference_exdistr}]);
length_pieces=alllength_pieces.(['length_pieces' file{reference_exdistr}]);

%I cut all Dna at a given time, I calculate the global fraction of replication 
%of the cutted fibers for all the genome and I do the mean; I stop the
%cutting when I find an average fraction of replication higher of the experimental
%reference and one lower; then I choose the lower one that is also the
%closer to the experimental reference.

%Until I reach the last time of the simulated S phase or I until I don't
%find an average fraction of replication higher of the experimental
%reference and one lower, I continue to cut at another time point of the
%simulation
%I start by cutting at the 'reference_thtime'
Drecon=zeros(length_DNA,1);
for u=1:num_DNA
%With the function cut, I extract the lengths to cut the genome with the
%same distribution of the experiment, then I take the positions on the
%genome random
pos_pieces=cutsame(num_pieces,length_pieces,length_DNA,unit);
Drecon=reconstruc(Dna(u).time(reference_thtime),length_DNA);
for i=1:num_pieces
    %The pieces are saved in the Dna structure at the relative Dna number
    %at the time cosidered
    Dna(u).time(reference_thtime).thDcut(i).fiber=Drecon(pos_pieces{i});
end
%I calculate the total fraction of replication of the cutted fibers      
 temptotfraction_rep(u)=sum(cellfun(@sum,{Dna(u).time(reference_thtime).thDcut.fiber}))/sum(cellfun(@length,{Dna(u).time(reference_thtime).thDcut.fiber}));
end

%I do the mean of the total replicated fraction of all the genomes
thfrac_choose(1,reference_thtime)=mean(temptotfraction_rep);

%I compare with 'reference_exfrac' and I decide to cut for lower or higher times 
%If thfrac_choose=reference_exfrac the reference_thtime remain just the
%same
if thfrac_choose(1,reference_thtime)<reference_exfrac 
    while reference_thtime~=mintime && thfrac_choose(1,reference_thtime)<reference_exfrac 
    reference_thtime=reference_thtime+1;
    %I cut at the new reference_thtime
    for u=1:num_DNA
    pos_pieces=cutsame(num_pieces,length_pieces,length_DNA,unit);
    Drecon=reconstruc(Dna(u).time(reference_thtime),length_DNA);
    for i=1:num_pieces
    Dna(u).time(reference_thtime).thDcut(i).fiber=Drecon(pos_pieces{i});
    end   
    temptotfraction_rep(u)=sum(cellfun(@sum,{Dna(u).time(reference_thtime).thDcut.fiber}))/sum(cellfun(@length,{Dna(u).time(reference_thtime).thDcut.fiber}));
    end
    thfrac_choose(1,reference_thtime)=mean(temptotfraction_rep);
    end
    %If reference_thtime=mintime the reference time will be the last one
    %considered, if reference_thtime~=mintime will be the oe before
    if reference_thtime~=mintime
    reference_thtime=reference_thtime-1;
    end
elseif thfrac_choose(1,reference_thtime)>reference_exfrac 
    %In this case the reference time will be the last one considered
    while reference_thtime~=1 && thfrac_choose(1,reference_thtime)>reference_exfrac 
    reference_thtime=reference_thtime-1;
    %I cut at the new reference_thtime
    for u=1:num_DNA
    pos_pieces=cutsame(num_pieces,length_pieces,length_DNA,unit);
    Drecon=reconstruc(Dna(u).time(reference_thtime),length_DNA);
    for i=1:num_pieces
    Dna(u).time(reference_thtime).thDcut(i).fiber=Drecon(pos_pieces{i});
    end   
    temptotfraction_rep(u)=sum(cellfun(@sum,{Dna(u).time(reference_thtime).thDcut.fiber}))/sum(cellfun(@length,{Dna(u).time(reference_thtime).thDcut.fiber}));
    end
    thfrac_choose(1,reference_thtime)=mean(temptotfraction_rep);
    end
end

%-------------->Analysis of the time of interest<-------------------

%I calculate separately for each time point frequency of inititiation
%and density of forks for each fiber and I sort by replicated fraction.
%Then I average the result of the different time points.

for s=1:num_DNA
    
%Calculate all I need in theoretical fibers
[Dna(s).time(reference_thtime).thDcut]=calculatestruc(Dna(s).time(reference_thtime).thDcut,num_pieces,thre1/unit,thre2/unit);   
thDcut=Dna(s).time(reference_thtime).thDcut;

%Analysis of single time points 
%I divide fiber in accordance with percentage of replication
%The bin are define in a way that edges(k)< x <=edges(k+1)
%The 0 is taken into account in the first bin
%The num_bin should be setted to have at least one fiber for each central bins

centersrep=1/(num_binrep*2):1/(num_binrep):(2*num_binrep-1)/(num_binrep*2);
histrep=zeros(1,length(centersrep));
histforknotave=zeros(1,length(centersrep));
histlengthnotave=zeros(1,length(centersrep));
neworiginsnotave=zeros(1,length(centersrep));
histunrepnotave=zeros(1,length(centersrep));
k=1;
for i=0:1/(num_binrep):(num_binrep-1)/num_binrep
    temp=([thDcut.fracfiber]>i & [thDcut.fracfiber]<=i+(1/(num_binrep)));
    %I calculate the number of fibers inn the bin
    histrep(k)=sum(temp);
    %I sum the number of all the forks for the fibers divided for
    %replicated fraction
    histforknotave(k)=sum(cellfun(@length,{thDcut(temp).threleft_forkwithside}))+sum(cellfun(@length,{thDcut(temp).threright_forkwithside}));
    %I sum the length of all the fibers divided for replicated fraction
    histlengthnotave(k)=sum(cellfun(@length,{thDcut(temp).fiber}));
    %I sum the number of all the new origins for the fibers divided for
    %replicated fraction
    my=cellfun(@(x) x<=(thre3/unit),{thDcut(temp).length_eyes},'UniformOutput',0);
    neworiginsnotave(k)=sum(cellfun(@sum,my));
    %I sum the lengths of all unreplicated parts for the fibers divided for
    %replicated fraction
    histunrepnotave(k)=histlengthnotave(k)-sum(cellfun(@sum,{thDcut(temp).fiber}));
    k=k+1;
end
%I add the quantities for the 0
temp=([thDcut.fracfiber]==0);
histrep(1)=histrep(1)+sum(temp);
histforknotave(1)=histforknotave(1)+sum(cellfun(@length,{thDcut(temp).threleft_forkwithside}))+sum(cellfun(@length,{thDcut(temp).threright_forkwithside}));
histlengthnotave(1)=histlengthnotave(1)+sum(cellfun(@length,{thDcut(temp).fiber}));
my=cellfun(@(x) x<=(thre3/unit),{thDcut(temp).length_eyes},'UniformOutput',0);
neworiginsnotave(1)=neworiginsnotave(1)+sum(cellfun(@sum,my));
histunrepnotave(1)=histunrepnotave(1)+histlengthnotave(1)-sum(cellfun(@sum,{thDcut(temp).fiber}));
%I calculate fork density and frequency of initiation and I put to 0 the NaN
%values; a point is NaN if there are not origins or forks for that bin.
histforkdens_notave=histforknotave./histlengthnotave;
histfreqinit_notave=neworiginsnotave./(histunrepnotave.*interval);
histforkdens_notave(isnan(histforkdens_notave))=0;
histfreqinit_notave(isnan(histfreqinit_notave))=0;
%I save the results in a structure
thresult.(['Dna' num2str(s)]).fraction_rep.centers=centersrep;
thresult.(['Dna' num2str(s)]).fraction_rep.histrep=histrep/sum(histrep);
thresult.(['Dna' num2str(s)]).fork_density.histforkdens_notave=histforkdens_notave;
thresult.(['Dna' num2str(s)]).freq_init.histfreqinit_notave=histfreqinit_notave;
end

%I calculate separately for each time point ETED, eye length and gap length distribution.
%Then I average the result of the different time points.
for s=1:num_DNA
    
thDcut=Dna(s).time(reference_thtime).thDcut;

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
my={thDcut.etedist};
eted=cell2mat(my(~cellfun(@isempty,my))');
my={thDcut.length_eyes};
eyes=cell2mat(my(~cellfun(@isempty,my))');
my={thDcut.gap_length};
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
thresult.(['Dna' num2str(s)]).eted.centers=centerseyes;
thresult.(['Dna' num2str(s)]).eted.histeted=histeted;
thresult.(['Dna' num2str(s)]).length_eyes.histeyes=histeyes;
thresult.(['Dna' num2str(s)]).length_gaps.histgaps=histgaps;
end

%I find mean and standard deviation for each fraction of replication
histrepcenter=zeros(1,num_DNA);
histfreqinit_notavecenter=zeros(1,num_DNA);
histforkdens_notavecenter=zeros(1,num_DNA);
histetedcenter=zeros(1,num_DNA);
histeyescenter=zeros(1,num_DNA);
histgapscenter=zeros(1,num_DNA);
histrepmean=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histrepstd=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histfreqinit_notavemean=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histfreqinit_notavestd=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histforkdens_notavemean=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histforkdens_notavestd=zeros(1,length(thresult.(['Dna' num2str(1)]).fraction_rep.centers));
histetedmean=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));
histetedstd=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));
histeyesmean=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));
histeyesstd=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));
histgapsmean=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));
histgapsstd=zeros(1,length(thresult.(['Dna' num2str(1)]).eted.centers));

for s=1:length(thresult.(['Dna' num2str(1)]).fraction_rep.centers)
%I put the values of replication for all the files in a single array
for i=1:num_DNA
    histrepcenter(i)=thresult.(['Dna' num2str(i)]).fraction_rep.histrep(s);
end
%I calculate mean and std
histrepmean(s)=mean(histrepcenter);
histrepstd(s)=std(histrepcenter)/sqrt(num_DNA);
%I put the values of fork density for all the files in a single array
for i=1:num_DNA
    histfreqinit_notavecenter(i)=thresult.(['Dna' num2str(i)]).freq_init.histfreqinit_notave(s);
end
histfreqinit_notavemean(s)=mean(histfreqinit_notavecenter);
histfreqinit_notavestd(s)=std(histfreqinit_notavecenter)/sqrt(num_DNA);
%I put the values of initiation frequency for all the files in a single array
for i=1:num_DNA
    histforkdens_notavecenter(i)=thresult.(['Dna' num2str(i)]).fork_density.histforkdens_notave(s);
end
histforkdens_notavemean(s)=mean(histforkdens_notavecenter);
histforkdens_notavestd(s)=std(histforkdens_notavecenter)/sqrt(num_DNA);

end
%I do the same for ETED, eye length and gap length distributions
for s=1:length(thresult.(['Dna' num2str(1)]).eted.centers)
%I put the values of eted distribution, for all the files in a single array
for i=1:num_DNA
    histetedcenter(i)=thresult.(['Dna' num2str(i)]).eted.histeted(s);
end
%I calculate mean and std
histetedmean(s)=mean(histetedcenter);
histetedstd(s)=std(histetedcenter)/sqrt(num_DNA);
%I put the values of eye lengths distribution for all the files in a single array
for i=1:num_DNA
    histeyescenter(i)=thresult.(['Dna' num2str(i)]).length_eyes.histeyes(s);
end
histeyesmean(s)=mean(histeyescenter);
histeyesstd(s)=std(histeyescenter)/sqrt(num_DNA);
%I put the values of gaps lengths distribution for all the files in a single array
for i=1:num_DNA
    histgapscenter(i)=thresult.(['Dna' num2str(i)]).length_gaps.histgaps(s);
end
histgapsmean(s)=mean(histgapscenter);
histgapsstd(s)=std(histgapscenter)/sqrt(num_DNA);

end
%I save all the results in the structure
thresult.finalresult.fraction_rep.mean=histrepmean;
thresult.finalresult.fraction_rep.std=histrepstd;
thresult.finalresult.fork_density.mean=histforkdens_notavemean;
thresult.finalresult.fork_density.std=histforkdens_notavestd;
thresult.finalresult.freq_init.mean=histfreqinit_notavemean;
thresult.finalresult.freq_init.std=histfreqinit_notavestd;
thresult.finalresult.eted.mean=histetedmean;
thresult.finalresult.eted.std=histetedstd;
thresult.finalresult.length_eyes.mean=histeyesmean;
thresult.finalresult.length_eyes.std=histeyesstd;
thresult.finalresult.length_gaps.mean=histgapsmean;
thresult.finalresult.length_gaps.std=histgapsstd;
thresult.finalresult.num_binrep=num_binrep;
thresult.finalresult.num_bineyes=num_bineyes;
thresult.finalresult.maxlength_bineyes=maxlength_bineyes;


%--------------->Calculation total quantities<-----------------------
%I calculate total frequency of initiation, fraction of replication and
%fork density at the given time
for s=1:num_DNA   
    
% I calculate the total replicated fraction, frequency of initiation and fork density for the time points
%"cellfun" apply the given function to each cell in cell array; the
%structure is coverted in cell array with the braket {}.
%The replicated fraction is calculated like the sum of all tha 1-blocks of
%all the fibers of a given experimantal time point divided by the sum of
%the lengths of the fibers.
lengthfiber(s)=sum(cellfun(@length,{Dna(s).time(reference_thtime).thDcut.fiber}));
sumfiber(s)=sum(cellfun(@sum,{Dna(s).time(reference_thtime).thDcut.fiber}));
totthfrac(s)=sumfiber(s)/lengthfiber(s);

%I calculate the new origins on all the fibers of a time point.
%In the vectors containing the lengths of eyes for each fiber, "cellfun" recognize the eye lengths
%less than thre3 (converted in kb) and put a 1 at that positions and 0 at
%the others position. my is a cell array, with a cell for each fiber. 
% 'UniformOutput' in cell fun allow to combine outputs of different size
my=cellfun(@(x) x<=(thre3/unit),{Dna(s).time(reference_thtime).thDcut.length_eyes},'UniformOutput',0);
neworigins(s)=sum(cellfun(@sum,my));

%I calculate the frequency of initiation (unit=1/kb*sec) and the fork density (1/kb) for each time
%point.
freq_init(s)=neworigins(s)/((lengthfiber(s)-sumfiber(s))*interval);
fork_density(s)=sum(cellfun(@length,{Dna(s).time(reference_thtime).thDcut.threleft_forkwithside})+cellfun(@length,{Dna(s).time(reference_thtime).thDcut.threright_forkwithside}))/lengthfiber(s);

end

thresult.finalresult.tot.fraction_rep=totthfrac;
thresult.finalresult.tot.freq_init=freq_init;
thresult.finalresult.tot.fork_density=fork_density;


%--------------->Calculation cost function<-----------------------
%I define the cost function as the normalized sum of squared errors(normalized to the mean of the experimental data that are not zero,
%  since if I consider also the zero this give more importance to the spacial graphs)
 meanone=exresult.finalresult.fraction_rep.mean;
 meantwo=exresult.finalresult.fork_density.mean;
 meanthree=exresult.finalresult.freq_init.mean;
 meanfour=exresult.finalresult.eted.mean;
 meanfourth=thresult.finalresult.eted.mean;
 meanfive=exresult.finalresult.length_eyes.mean;
 meansix=exresult.finalresult.length_gaps.mean;
 %Cost for time point quantitities
  costpart.residual.one=sum((thresult.finalresult.fraction_rep.mean(meanone~=0)-exresult.finalresult.fraction_rep.mean(meanone~=0)).^2/(mean(meanone(meanone~=0)).^2));
 costpart.residual.two=sum((thresult.finalresult.fork_density.mean(meantwo~=0)-exresult.finalresult.fork_density.mean(meantwo~=0)).^2/(mean(meantwo(meantwo~=0)).^2));
 costpart.residual.three=sum((thresult.finalresult.freq_init.mean(meanthree~=0)-exresult.finalresult.freq_init.mean(meanthree~=0)).^2/(mean(meanthree(meanthree~=0)).^2));
  costpart.residual.four=sum((log(thresult.finalresult.eted.mean(meanfour~=0 & meanfourth~=0))-log(exresult.finalresult.eted.mean(meanfour~=0 & meanfourth~=0))).^2/(mean(log(meanfour(meanfour~=0 & meanfourth~=0))).^2));
 costpart.residual.five=sum((thresult.finalresult.length_eyes.mean(meanfive~=0)-exresult.finalresult.length_eyes.mean(meanfive~=0)).^2/(mean(meanfive(meanfive~=0)).^2));
 costpart.residual.six=sum((thresult.finalresult.length_gaps.mean(meansix~=0)-exresult.finalresult.length_gaps.mean(meansix~=0)).^2/(mean(meansix(meansix~=0)).^2));
 %Cost for total quantities time point
 costpart.residual.seven=(mean(thresult.finalresult.tot.fraction_rep)-mean(exresult.finalresult.tot.fraction_rep))^2/mean(exresult.finalresult.tot.fraction_rep)^2;
 costpart.residual.eight=(mean(thresult.finalresult.tot.fork_density)-mean(exresult.finalresult.tot.fork_density))^2/mean(exresult.finalresult.tot.fork_density)^2;
 costpart.residual.nine=(mean(thresult.finalresult.tot.freq_init)-mean(exresult.finalresult.tot.freq_init))^2/mean(exresult.finalresult.tot.freq_init)^2;

cost.residual=costpart.residual.one+costpart.residual.two+costpart.residual.three+costpart.residual.four+costpart.residual.five+costpart.residual.six+costpart.residual.seven+...
                costpart.residual.eight+costpart.residual.nine;

%I define the cost function as the correlation coefficient
costpart.corrcoef.one=corrcoef(thresult.finalresult.fraction_rep.mean,exresult.finalresult.fraction_rep.mean);
 costpart.corrcoef.two=corrcoef(thresult.finalresult.fork_density.mean,exresult.finalresult.fork_density.mean);
 costpart.corrcoef.three=corrcoef(thresult.finalresult.freq_init.mean,exresult.finalresult.freq_init.mean);
 costpart.corrcoef.four=corrcoef(thresult.finalresult.eted.mean,exresult.finalresult.eted.mean);
 costpart.corrcoef.five=corrcoef(thresult.finalresult.length_eyes.mean,exresult.finalresult.length_eyes.mean);
 costpart.corrcoef.six=corrcoef(thresult.finalresult.length_gaps.mean,exresult.finalresult.length_gaps.mean);
 
cost.corrcoef=-costpart.corrcoef.one(1,2)-costpart.corrcoef.two(1,2)-costpart.corrcoef.three(1,2)-costpart.corrcoef.four(1,2)-costpart.corrcoef.five(1,2)-costpart.corrcoef.six(1,2);


end