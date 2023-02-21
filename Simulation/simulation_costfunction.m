function costfunction=simulation_costfunction(var,allexDcut,allnum_pieces,alllength_pieces,file,exresult)

%With this function I calculate the cost for the particular set of
%variables, to be used in the genetic algorithm. 
%MM5 model described in Ciardo et al. 2021 (https://doi.org/10.3390/genes12081224)
% is used here to simulate the replication process.

%no automatic integer setting in the genetic algorithm
%I convert the variable that are integers in integers
limit_factorinit=round(var(1));
ratej=round(var(2));
prob_globin=var(3);
prob_globout=var(4);
prob_loc=var(5);
dist_box=round(var(6));
teta=var(7);


% INFO ON REPLICATION STEPS
%->I elongate before initiate to not elongate the origins added in the
%  same step and to not initiate new origins in the positions where I'm
%  elongating at the same step.
%->The liming factor released from the merges at the given step is used
%  to initiate at the next step;

%General variables
unit=allexDcut.(['exDcut' file{1}])(1).unit_block; %Define the number of base pair (bp) for each block of the genome
fraction_max=0.25; %Fraction of replication at which I block the replication
num_DNA=15; %Number of DNA I simulate 
time_max=30; %Time at which I block the replication 
v=8.3; %velocity of fork in bp/s
length_DNA=100000000/unit; %DNA lenght
timeunit_insec=100000000/(length_DNA*v); %I know velocity and that I replicate a block/step
origin_density=1/2.3 ; %kb^(-1)


% ---------------------> Simulation of replication<-----------------------
%#########################################################################
%#########################################################################

for f=1:num_DNA
time=0; %Put to zero and increase in the fist line of while
        %Otherwise at the last step it increase and then stop
free_part=0; %Particle of limiting factor free at every step
fraction_rep=0; %Fraction  of replicated DNA
limit_factor=limit_factorinit; %I initialize the limiting factor to the initial quantity
pos_init=[]; %Matrix that contain positions with probability of initiation less than prob_init
left_fork=[]; %Matrix with position of left forks
right_fork=[]; %Matrix with position of right forks
point_fusion=[]; %Matrix with position of the poit of fusion

D=zeros(length_DNA,1); %Strand of DNA

%I define potential origins
potential_origins_position=ceil(length_DNA*rand(round(length_DNA*origin_density),1));
potential_origins=zeros(length_DNA,1);
potential_origins(potential_origins_position)=1;

%I set the different probability of initiation in the two regions of the
%genome defined by teta (the first regions represent domains; the second represent the
%region outside domains)
prob_init=prob_globin.*ones(length_DNA,1);
prob_init(round(teta*length_DNA)+1:end)=prob_globout;

while fraction_rep<fraction_max && time<time_max
    
   time=time+1;
   
   %Elongation and position forks
   %I find the new forks by finding the sequences 0 1 or 1 0 in D 
   %I elongate before initiate to not elongate the origins added in the
   %same step and to not initiate new origins in the positions where I'm
   %elongating at the same time.
    D(left_fork-1)=1;
    D(right_fork+1)=1;
    
    %Points of merge due to the elongation at this step
    if ~isempty(right_fork) && ~isempty(left_fork) %If the two vector are not empty (without this condition have problem)
    if left_fork(1)==2 %Control on left side
        free_part=free_part+0.5; end
    if right_fork(length(right_fork))==length_DNA-1  %Cotrol on right side
        free_part=free_part+0.5; end
    if left_fork(1)<=right_fork(1)  %If the first fork is left fork...
        %temp contain the positions in right fork array where the right fork 
        %differ from left forks (except the fist one) for one or zero
        %The point of merge are the positions in the DNA that are contained
        %in these positions
        temp=(left_fork(2:end,1)-right_fork(1:(length(left_fork)-1),1))<=3; 
        point_fusion=right_fork(temp);
    else %If the first fork is a right fork...
        %temp contain the positions in right fork array where the right fork 
        %differ from left forks(included the fist one) for one or zero
        %The point of merge are the positions in the DNA that are contained
        %in these positions
        temp=(left_fork(1:end,1)-right_fork(1:length(left_fork),1))<=3; 
        point_fusion=right_fork(temp);
    end
    end
    
  %Particle free for the merges at this step are used for the initiations
  %at the next one; they are added to the limiting factor at the end.
   free_part=free_part+length(point_fusion);
    
  %Clustering
   %Action of boxes
    prob=zeros(length_DNA,1);
    if dist_box>0
    for i=1:length(left_fork)
        action_box=[(left_fork(i)-1-dist_box),left_fork(i)-2]; %I put -2 because the code elongates before; I depart from the position before the fork
        if action_box(1)<1  action_box(1)=1; end
        if action_box(2)>length_DNA  action_box(2)=length_DNA; end
        prob(action_box(1):action_box(2))=prob_loc; 
    end
    for i=1:length(right_fork)
        action_box=[right_fork(i)+2,(right_fork(i)+1+dist_box)]; %I put +2 because the code elongates before; I depart from the position after the fork
        if action_box(1)<1  action_box(1)=1; end
        if action_box(2)>length_DNA  action_box(2)=length_DNA; end
        prob(action_box(1):action_box(2))=prob_loc; 
    end
    end
   
   
   %Initiation
   %We consider the potential origins in vast excess; 
   %all blocks can initiate replication
   %Find the positions where the limiting factor encounter the strand by comparing 
   %the randomly extracted probability with a chosen probability
   %If we have more chosen positions than the number of particle free, then
   %we choose randomly in the array 'pos_init' the right number of
   %positions
   
   pos_init=find(rand(length_DNA,1)<(prob+prob_init)); %With prob I take into account the action of P_box
   
   if (length(pos_init)>limit_factor && limit_factor>0)
           [~,a]=sort(rand(1,length(pos_init))); %In 'a' we store the positions of sorted elements
           pos_init=pos_init(a(1:limit_factor)); %I take in 'a' the fist positions in number equal to the free particles
   elseif limit_factor<=0
           pos_init=[];
   end
   
   %I initiate the new origins
   for j=1:length(pos_init) 
      if  (D(pos_init(j))==0 && potential_origins(pos_init(j))==1) %If the position in DNA is not already replicated and the position is a potential origin, I initiate origin
          D(pos_init(j))=1;
          limit_factor=limit_factor-1;
          %We must free a particle when we put an origin near a fork; 
          %we have a problem also when the origin is at the side
          %I do two 'if' because a new origin could be close both to the side
          %and to a fork
          if pos_init(j)==1 free_part=free_part+0.5; 
          elseif D(pos_init(j)-1)==1 free_part=free_part+1; end
          if pos_init(j)==length_DNA free_part=free_part+0.5;
          elseif D(pos_init(j)+1)==1 free_part=free_part+1; end
      end
   end
   
  %I add here free particles because I don't want to use them for
  %initiate new origins at the same step
  %Add the free particle to limiting particle
  limit_factor=limit_factor+floor(free_part);
  
  %Put to zero free particle for the next step
  free_part=0;
  
  %I add new particles of limit_factor with a rate of ratej/s
  limit_factor=limit_factor+ratej*round(timeunit_insec);
  
  %Find forks
  left_fork=(findstr(D',[0 1])+1)';
  right_fork=(findstr(D',[1 0]))';
  Dna(f).time(time).left_fork=left_fork;
  Dna(f).time(time).right_fork=right_fork;
  
  %Fraction DNA replicated
  fraction_rep=sum(D)/length_DNA;
  Dna(f).time(time).fraction_rep=fraction_rep;
    
end
end
%#########################################################################
%#########################################################################
%-------------------------------------------------------------------------



%--------------->Analysis of the time of interest<-----------------------
%#########################################################################
%#########################################################################

[thresult,cost,costpart,reference_thtime]=calculatethresult(Dna,file,exresult,allnum_pieces,alllength_pieces,num_DNA,length_DNA,unit,v,0);
costfunction=cost.residual;
%#########################################################################
%#########################################################################
%-------------------------------------------------------------------------

end

