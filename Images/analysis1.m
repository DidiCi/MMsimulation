close all
clear all

path='../Result_analysis/'; 

format long 
cost_func='NMSE';
y_ax={'N_0','J (s^{-1})','P_{in}','P_{out}','P_{local}','d (kb)','\theta'};
sample_folder='Condition1';
num_sample=3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%% Plot comparison experiment and simulation %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
Color_simu='r'; 
Color_exp='k'; 
Mark_size=2.5;
Line_width_single=0.3;
Line_width=1;
Font_Size=10;
load(['../Data_extraction/' sample_folder '/exresult']);
load(['../Data_extraction/' sample_folder '/file']);

scrsz = get(0,'ScreenSize');
figure('Position',[10 10 scrsz(3)*9/10 scrsz(4)*9/10])

subplot(2,3,1);
hold all
axis([0 1 0 0.85])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized number of fibers','FontSize',Font_Size,'FontName','Arial');
%Save results of all simulations in one matrix
for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    fraction_rep(i,:)=thresult.finalresult.fraction_rep.mean;
end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fraction_rep.mean,exresult.finalresult.fraction_rep.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fraction_rep.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(fraction_rep),std(fraction_rep)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(fraction_rep);
exp_values=exresult.finalresult.fraction_rep.mean;
gof_rep=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,2);
hold all
axis([0 1 0 0.14])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
ylabel('Averaged forks density (1/kb)','FontSize',Font_Size,'FontName','Arial');
for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    fork_density(i,:)=thresult.finalresult.fork_density.mean;
end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fork_density.mean,exresult.finalresult.fork_density.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fork_density.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(fork_density),std(fork_density)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(fork_density);
exp_values=exresult.finalresult.fork_density.mean;
gof_fork=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,3);
hold all
axis([0 1 0 0.0004])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
 ylabel('Averaged I(f)(1/kb*sec)','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    freq_init(i,:)=thresult.finalresult.freq_init.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.freq_init.mean,exresult.finalresult.freq_init.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.freq_init.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(freq_init),std(freq_init)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(freq_init);
exp_values=exresult.finalresult.freq_init.mean;
gof_freq=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,4);
hold all
axis([0 80 -Inf 0.26])
xlabel('ETED (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized ETED distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    eted(i,:)=thresult.finalresult.eted.mean;
 end
 if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.eted.mean,exresult.finalresult.eted.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
 else
     g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.eted.mean,8,Color_exp,'linewidth',Line_width);
 end
k=errorbar(exresult.(file{1}).eted.centers,mean(eted),std(eted)/sqrt(num_sample),'color','r','linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(eted);
exp_values=exresult.finalresult.eted.mean;
gof_eted=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,5);
hold all
axis([0 80 -Inf 0.6])
xlabel('Gaps length (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized gaps length distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    length_eyes(i,:)=thresult.finalresult.length_eyes.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.length_eyes.mean,exresult.finalresult.length_eyes.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.length_eyes.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(exresult.(file{1}).eted.centers,mean(length_eyes),std(length_eyes)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(length_eyes);
exp_values=exresult.finalresult.length_eyes.mean;
gof_eyes=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,6)
hold all
axis([0 80 -Inf 0.7])
xlabel('Eye length (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized eye length distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    length_gaps(i,:)=thresult.finalresult.length_gaps.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.length_gaps.mean,exresult.finalresult.length_gaps.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.length_gaps.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(exresult.(file{1}).eted.centers,mean(length_gaps),std(length_gaps)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(length_gaps);
exp_values=exresult.finalresult.length_gaps.mean;
gof_gaps=goodnessOfFit(th_values',exp_values',cost_func);
legend([g(1),k(1)],{'Experiment','Simulation'},'FontSize',8,'FontName','Arial');
hold off

gof_tot=(gof_rep+gof_fork+gof_freq+gof_eted+gof_eyes+gof_gaps)/6;

im_paper(['Results/' sample_folder])

xlswrite(['Results/Stat_' sample_folder '.xls'],{sample_folder,'Goodness of Fit (average curves)'},'A1:B1');
xlswrite(['Results/Stat_' sample_folder '.xls'],{'norm. N of fibers';'fork density';'I(dt)';'ETED';'Gap length';'Eye length'},'A2:A7');
xlswrite(['Results/Stat_' sample_folder '.xls'],[gof_rep;gof_fork;gof_freq;gof_eted;gof_eyes;gof_gaps],'B2:B7');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%% Plot comparison experiment and simulation %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% with single plots %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear fraction_rep fork_density freq_init eted length_eyes length_gaps
Color_simu='r'; 
Color_exp='k'; 
Mark_size=2.5;
Line_width_single=0.3;
Line_width=1;
Font_Size=10;

scrsz = get(0,'ScreenSize');
figure('Position',[10 10 scrsz(3)*9/10 scrsz(4)*9/10])

subplot(2,3,1);
hold all
axis([0 1 0 0.85])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized number of fibers','FontSize',Font_Size,'FontName','Arial');
%Save results of all simulations in one matrix
for i=1:num_sample
    th_values=thresult.finalresult.fraction_rep.mean;
    exp_values=exresult.finalresult.fraction_rep.mean;
    gof_rep_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).fraction_rep.centers,thresult.finalresult.fraction_rep.mean,'color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    fraction_rep(i,:)=thresult.finalresult.fraction_rep.mean;
end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fraction_rep.mean,exresult.finalresult.fraction_rep.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fraction_rep.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(fraction_rep),std(fraction_rep)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(fraction_rep);
exp_values=exresult.finalresult.fraction_rep.mean;
gof_rep=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,2);
hold all
axis([0 1 0 0.14])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
ylabel('Averaged forks density (1/kb)','FontSize',Font_Size,'FontName','Arial');
for i=1:num_sample
    th_values=thresult.finalresult.fork_density.mean;
    exp_values=exresult.finalresult.fork_density.mean;
    gof_fork_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).fraction_rep.centers,thresult.finalresult.fork_density.mean,'Color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    fork_density(i,:)=thresult.finalresult.fork_density.mean;
end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fork_density.mean,exresult.finalresult.fork_density.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.fork_density.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(fork_density),std(fork_density)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(fork_density);
exp_values=exresult.finalresult.fork_density.mean;
gof_fork=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,3);
hold all
axis([0 1 0 0.0004])
xlabel('Fraction of replication','FontSize',Font_Size,'FontName','Arial');
 ylabel('Averaged I(f)(1/kb*sec)','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    th_values=thresult.finalresult.freq_init.mean;
    exp_values=exresult.finalresult.freq_init.mean;
    gof_freq_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).fraction_rep.centers,thresult.finalresult.freq_init.mean,'Color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    freq_init(i,:)=thresult.finalresult.freq_init.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.freq_init.mean,exresult.finalresult.freq_init.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).fraction_rep.centers,exresult.finalresult.freq_init.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(thresult.(['Dna' num2str(1)]).fraction_rep.centers,mean(freq_init),std(freq_init)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(freq_init);
exp_values=exresult.finalresult.freq_init.mean;
gof_freq=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,4);
hold all
axis([0 80 -Inf 0.26])
xlabel('ETED (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized ETED distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    th_values=thresult.finalresult.eted.mean;
    exp_values=exresult.finalresult.eted.mean;
    gof_eted_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).eted.centers,thresult.finalresult.eted.mean,'Color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    eted(i,:)=thresult.finalresult.eted.mean;
 end
 if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.eted.mean,exresult.finalresult.eted.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
 else
     g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.eted.mean,8,Color_exp,'linewidth',Line_width);
 end
k=errorbar(exresult.(file{1}).eted.centers,mean(eted),std(eted)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(eted);
exp_values=exresult.finalresult.eted.mean;
gof_eted=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,5);
hold all
axis([0 80 -Inf 0.6])
xlabel('Gaps length (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized gaps length distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    th_values=thresult.finalresult.length_eyes.mean;
    exp_values=exresult.finalresult.length_eyes.mean;
    gof_eyes_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).eted.centers,thresult.finalresult.length_eyes.mean,'Color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    length_eyes(i,:)=thresult.finalresult.length_eyes.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.length_eyes.mean,exresult.finalresult.length_eyes.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.length_eyes.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(exresult.(file{1}).eted.centers,mean(length_eyes),std(length_eyes)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(length_eyes);
exp_values=exresult.finalresult.length_eyes.mean;
gof_eyes=goodnessOfFit(th_values',exp_values',cost_func);
hold off

subplot(2,3,6)
hold all
axis([0 80 -Inf 0.7])
xlabel('Eye length (kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Normalized eye length distribution','FontSize',Font_Size,'FontName','Arial');
 for i=1:num_sample
    th_values=thresult.finalresult.length_gaps.mean;
    exp_values=exresult.finalresult.length_gaps.mean;
    gof_gaps_single(i)=goodnessOfFit(th_values',exp_values',cost_func);
    h=plot(thresult.(['Dna' num2str(1)]).eted.centers,thresult.finalresult.length_gaps.mean,'Color',Color_simu,'LineWidth',Line_width_single);
        
    load([path  sample_folder '/thresult' num2str(i)]);
    length_gaps(i,:)=thresult.finalresult.length_gaps.mean;
 end
if length(exresult.finalresult.tot.fraction_rep)>1
    g=errorbar(exresult.(file{1}).eted.centers,exresult.finalresult.length_gaps.mean,exresult.finalresult.length_gaps.error,'Color',Color_exp,'LineStyle','none','Marker','o','MarkerSize',Mark_size,'linewidth',Line_width);
else
    g=scatter(exresult.(file{1}).eted.centers,exresult.finalresult.length_gaps.mean,8,Color_exp,'linewidth',Line_width);
end
k=errorbar(exresult.(file{1}).eted.centers,mean(length_gaps),std(length_gaps)/sqrt(num_sample),'color',Color_simu,'linewidth',Line_width,'LineStyle','--');
%gof
th_values=mean(length_gaps);
exp_values=exresult.finalresult.length_gaps.mean;
gof_gaps=goodnessOfFit(th_values',exp_values',cost_func);
legend([g(1),k(1)],{'Experiment','Simulation'},'FontSize',8,'FontName','Arial');
hold off

gof_tot=(gof_rep+gof_fork+gof_freq+gof_eted+gof_eyes+gof_gaps)/6;
gof_tot_single=mean([gof_rep_single;gof_fork_single;gof_freq_single;gof_eted_single;gof_eyes_single;gof_gaps_single]);

xlswrite(['Results/Stat_' sample_folder '.xls'],{'Goodness of Fit (single curves)',''},'F1:G1');
xlswrite(['Results/Stat_' sample_folder '.xls'],{'norm. N of fibers','fork density','I(dt)','ETED','Gap length','Eye length'},'F2:K2');
xlswrite(['Results/Stat_' sample_folder '.xls'],[gof_rep_single',gof_fork_single',gof_freq_single',gof_eted_single',gof_eyes_single',gof_gaps_single'],['F3:K' num2str(2+num_sample)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%% Plot total quantities %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Position',[10 10 scrsz(3)*9/10 scrsz(4)*9/10])
subplot(1,3,1);

hold on
xlabel('Total replicated fraction','FontSize',Font_Size,'FontName','Arial');
ylabel('Fraction');
axis([0 Inf 0 Inf])
th_totrep=[];
for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    th_totrep=[th_totrep,mean(thresult.finalresult.tot.fraction_rep)];
end
SP=mean(exresult.finalresult.tot.fraction_rep); 
[b,a]=hist(th_totrep,15);
bar(a,b/num_sample,'facecolor',Color_simu);
plot([SP SP],[0 max(b/num_sample)+1/num_sample],'k','linewidth',2)
legend('Simulation','Experiment');
hold off

subplot(1,3,2);
hold on
xlabel('Total freq. of initiation (1/(kb*sec))','FontSize',Font_Size,'FontName','Arial');
ylabel('Fraction');
axis([-Inf Inf 0 Inf])
th_totfreq=[];
for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    th_totfreq=[th_totfreq,mean(thresult.finalresult.tot.freq_init)];
end
[b,a]=hist(th_totfreq,15);
bar(a,b/num_sample,'facecolor',Color_simu);
SP=mean(exresult.finalresult.tot.freq_init); 
plot([SP SP],[0 max(b/num_sample)+1/num_sample],'k','linewidth',2)
hold off

subplot(1,3,3);
hold on
xlabel('Total fork density (1/kb)','FontSize',Font_Size,'FontName','Arial');
ylabel('Fraction');
axis([-Inf Inf 0 Inf])
th_totfork=[];
for i=1:num_sample
    load([path  sample_folder '/thresult' num2str(i)]);
    th_totfork=[th_totfork,mean(thresult.finalresult.tot.fork_density)];
end
SP=mean(exresult.finalresult.tot.fork_density); 
[b,a]=hist(th_totfork,15);
bar(a,b/num_sample,'facecolor',Color_simu);
plot([SP SP],[0 max(b/num_sample)+1/num_sample],'k','linewidth',2)
hold off


