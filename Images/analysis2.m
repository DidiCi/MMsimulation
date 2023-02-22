close all
clear all
addpath('../Functions') 

path='../Simulation/'; 

format long
y_ax={'N_0','J (s^{-1})','P_{in}','P_{out}','P_{local}','d (kb)','\theta'};
sample_folder={'Condition1' 'Condition2'};
title_samples={'Ctrl','Depl'};
num_sample=[3 3];

for i=1:num_sample(1)
    load([path sample_folder{1} '/garesult' num2str(i) '.mat'] );
    sample1_var(i,:)=garesult.x;
    sample1_cost(i)=garesult.fval;
end


for i=1:num_sample(2)
    load([path sample_folder{2} '/garesult' num2str(i) '.mat']);
    sample2_var(i,:)=garesult.x;
    sample2_cost(i)=garesult.fval;
end


%Boxplots  
scrsz = get(0,'ScreenSize');
figure('Position',[5 5 scrsz(3)*15/10 scrsz(4)*4/10])

linewidth_num=1.5;

for i=1:7
    subplot(1,7,i)
    hold on
    errorbar(1,mean(sample1_var(:,i)),std(sample1_var(:,i),1)/sqrt(num_sample(1)),'Marker','o','Color',[.6 .6 .6],'LineWidth',linewidth_num);
    errorbar(2,mean(sample2_var(:,i)),std(sample2_var(:,i),1)/sqrt(num_sample(2)),'Marker','o','Color','black','LineWidth',linewidth_num);
    axis([0 3 -Inf Inf])
    ylabel(y_ax{i},'FontSize',14,'FontName','Arial');

    up=max([mean(sample1_var(:,i))+std(sample1_var(:,i),1)/sqrt(num_sample(1)),mean(sample2_var(:,i))+std(sample2_var(:,i),1)/sqrt(num_sample(2))]);
    down=min([mean(sample1_var(:,i))-std(sample1_var(:,i),1)/sqrt(num_sample(1)),mean(sample2_var(:,i))-std(sample2_var(:,i),1)/sqrt(num_sample(2))]);
    axis([0 3 down-0.1*down up+0.1*up])
    set(gca,'XTick',[1,2])
    set(gca,'XTickLabel',title_samples)
    
    [p1(i),h1(i)] =ranksum(sample1_var(:,i),sample2_var(:,i));
    [h2(i),p2(i)] = kstest2(sample1_var(:,i),sample2_var(:,i));
    chi_2(i)=(mean(sample1_var(:,i))-mean(sample2_var(:,i)))^2/((std(sample1_var(:,i),1)/sqrt(num_sample(1)))^2+(std(sample2_var(:,i),1)/sqrt(num_sample(2)))^2);
    mean_values1(i)=mean(sample1_var(:,i));
    mean_values2(i)=mean(sample2_var(:,i));
    std_values1(i)=std(sample1_var(:,i),1)/sqrt(num_sample(1));
    std_values2(i)=std(sample2_var(:,i),1)/sqrt(num_sample(2));
    hold off
end

xlswrite('Results/Stat_boxplot.xls',{'N0','J (s^-1)','Pin','Pout','Plocal','d(kb)','theta'},'B1:H1');
xlswrite('Results/Stat_boxplot.xls',{'Mean cond1';'Std cond1';'Mean cond2';'Std cond2';'chi2';'ranksum h';'ranksum p';'ks h';'ks p'},'A2:A10');
xlswrite('Results/Stat_boxplot.xls',mean_values1,'B2:H2');
xlswrite('Results/Stat_boxplot.xls',std_values1,'B3:H3');
xlswrite('Results/Stat_boxplot.xls',mean_values2,'B4:H4');
xlswrite('Results/Stat_boxplot.xls',std_values2,'B5:H5');
xlswrite('Results/Stat_boxplot.xls',chi_2,'B6:H6');
xlswrite('Results/Stat_boxplot.xls',h1,'B7:H7');
xlswrite('Results/Stat_boxplot.xls',p1,'B8:H8');
xlswrite('Results/Stat_boxplot.xls',h2,'B9:H9');
xlswrite('Results/Stat_boxplot.xls',p2,'B10:H10');

im_paper2(['Results/boxplot' sample_folder{1} sample_folder{2}])
