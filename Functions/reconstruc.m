function [Drecon]=reconstruc(Dna,length_DNA)

if length(Dna.left_fork)>1 || length(Dna.right_fork)>1
if Dna.left_fork(1)<=Dna.right_fork(1) %If the first fork is left fork...
          Drecon(1:Dna.left_fork(1)-1,1)=zeros(Dna.left_fork(1)-1,1);
          for k=1:length(Dna.right_fork)-1
          Drecon(Dna.left_fork(k):Dna.right_fork(k),1)=ones(Dna.right_fork(k)-Dna.left_fork(k)+1,1);
          Drecon(Dna.right_fork(k)+1:Dna.left_fork(k+1)-1,1)=zeros(Dna.left_fork(k+1)-1-Dna.right_fork(k),1);
          end
          if Dna.left_fork(end)<=Dna.right_fork(end) %If the last fork is right fork...
          Drecon(Dna.left_fork(end):Dna.right_fork(end),1)=ones(Dna.right_fork(end)-Dna.left_fork(end)+1,1);  
          Drecon(Dna.right_fork(end)+1:length_DNA,1)=zeros(length_DNA-Dna.right_fork(end),1);
          else %If the last fork is left fork...
          Drecon(Dna.left_fork(end-1):Dna.right_fork(end),1)=ones(Dna.right_fork(end)-Dna.left_fork(end-1)+1,1);  
          Drecon(Dna.right_fork(end)+1:Dna.left_fork(end)-1,1)=zeros(Dna.left_fork(end)-1-Dna.right_fork(end),1);
          Drecon(Dna.left_fork(end):length_DNA,1)=ones(length_DNA-Dna.left_fork(end)+1,1);
          end
      else %If the first fork is right fork...
          Drecon(1:Dna.right_fork(1),1)=ones(Dna.right_fork(1),1);
          Drecon(Dna.right_fork(1)+1:Dna.left_fork(1)-1,1)=zeros(Dna.left_fork(1)-1-Dna.right_fork(1),1);
          for k=2:length(Dna.right_fork)-1
          Drecon(Dna.left_fork(k-1):Dna.right_fork(k),1)=ones(Dna.right_fork(k)-Dna.left_fork(k-1)+1,1);
          Drecon(Dna.right_fork(k)+1:Dna.left_fork(k)-1,1)=zeros(Dna.left_fork(k)-1-Dna.right_fork(k),1);
          end
          if Dna.left_fork(end)<=Dna.right_fork(end) %If the last fork is right fork...
          Drecon(Dna.left_fork(end):Dna.right_fork(end),1)=ones(Dna.right_fork(end)-Dna.left_fork(end)+1,1);  
          Drecon(Dna.right_fork(end)+1:length_DNA,1)=zeros(length_DNA-Dna.right_fork(end),1);
          else %If the last fork is left fork...
          Drecon(Dna.left_fork(end-1):Dna.right_fork(end),1)=ones(Dna.right_fork(end)-Dna.left_fork(end-1)+1,1);  
          Drecon(Dna.right_fork(end)+1:Dna.left_fork(end)-1,1)=zeros(Dna.left_fork(end)-1-Dna.right_fork(end),1);
          Drecon(Dna.left_fork(end):length_DNA,1)=ones(length_DNA-Dna.left_fork(end)+1,1);
          end
end
elseif length(Dna.left_fork)==1 && length(Dna.right_fork)==1
      if Dna.left_fork(1)<=Dna.right_fork(1) %If the first fork is left fork...
          Drecon(1:Dna.left_fork(1)-1,1)=zeros(Dna.left_fork(1)-1,1);
          Drecon(Dna.left_fork(1):Dna.right_fork(1),1)=ones(Dna.right_fork(1)-Dna.left_fork(1)+1,1);
          Drecon(Dna.right_fork(end)+1:length_DNA,1)=zeros(length_DNA-Dna.right_fork(end),1);
      else %If the first fork is right fork...
          Drecon(1:Dna.right_fork(1),1)=ones(Dna.right_fork(1),1);
          Drecon(Dna.right_fork(1)+1:Dna.left_fork(1)-1,1)=zeros(Dna.left_fork(1)-1-Dna.right_fork(1),1);
          Drecon(Dna.left_fork(end):length_DNA,1)=ones(length_DNA-Dna.left_fork(end)+1,1);
      end
elseif ~isempty(Dna.left_fork) && isempty(Dna.right_fork)
    Drecon(1:Dna.left_fork(1)-1,1)=zeros(Dna.left_fork(1)-1,1);
    Drecon(Dna.left_fork(1):length_DNA,1)=ones(length_DNA-Dna.left_fork(end)+1,1);
elseif isempty(Dna.left_fork) && ~isempty(Dna.right_fork)
    Drecon(1:Dna.right_fork(1),1)=ones(Dna.right_fork(1),1);
    Drecon(Dna.right_fork(1)+1:length_DNA,1)=zeros(length_DNA-Dna.right_fork(1),1);
elseif Dna.fraction_rep==0
    Drecon=zeros(length_DNA,1);
elseif Dna.fraction_rep==1
    Drecon=ones(length_DNA,1);
end
end