function [Dcut]=calculatestruc(Dcut,num_pieces,thre1,thre2)
%This function calculate all the characteristics of the fibers
%(experimental or theoretical) and save them in the structure exDcut.
%I apply first the threshold on the gaps (thre1) and then the threshold on
%the eyes length (thre2).

%Calculate the replicated fraction for each fiber
for i=1:num_pieces
    Dcut(i).fracfiber=sum(Dcut(i).fiber)/length(Dcut(i).fiber);
end  

%I find all forks 
%The forks are recognize by the sequences [0 1] and [1 0]. The position of
%the fork is cosidered on the 1.
for i=1:num_pieces
    
%With "rough" I indicate all the forks present on the fiber, before have
%applied the thresholding.
Dcut(i).roughleft_fork=(findstr(Dcut(i).fiber',[0 1])+1)';
Dcut(i).roughright_fork=(findstr(Dcut(i).fiber',[1 0]))';
%With "thre" I indicate the forks after thresholding. They are modified at
%each step to apply the thresholding.
Dcut(i).threleft_fork=Dcut(i).roughleft_fork;
Dcut(i).threright_fork=Dcut(i).roughright_fork;

%---------->I combine eyes divided by gaps smaller then thre1<----------
%If on the fiber there is at least an eye...
if ~isempty(Dcut(i).threleft_fork) && ~isempty(Dcut(i).threright_fork)
    %If the fist fork is a left fork
if Dcut(i).threleft_fork(1)<=Dcut(i).threright_fork(1)  
    % In "temp" -1 is to have the right length of gaps
    temp=find((Dcut(i).threleft_fork(2:end,1)-Dcut(i).threright_fork(1:(length(Dcut(i).threleft_fork)-1),1)-1)<thre1); 
    Dcut(i).threleft_fork(temp+1)=[];
    Dcut(i).threright_fork(temp)=[];
else %If the first fork is a right fork
    temp=find((Dcut(i).threleft_fork(1:end,1)-Dcut(i).threright_fork(1:length(Dcut(i).threleft_fork),1)-1)<thre1); 
    Dcut(i).threleft_fork(temp)=[];
    Dcut(i).threright_fork(temp)=[];
end
end
%>----------------------------------------------------------------------<
end


%------------------->I remove eyes smaller then thre2<-------------------
for i=1:num_pieces
%If on the fiber there is at least an eye...
if ~isempty(Dcut(i).threright_fork) && ~isempty(Dcut(i).threleft_fork)
    %If the fist fork is a left fork
    if Dcut(i).threleft_fork(1)<=Dcut(i).threright_fork(1)  
    %In "length_ayes" +1 is to have the right length of eyes
    length_eyes=(Dcut(i).threright_fork(1:end)-Dcut(i).threleft_fork(1:length(Dcut(i).threright_fork)))+1;
    Dcut(i).threleft_fork(length_eyes<thre2)=[];
    Dcut(i).threright_fork(length_eyes<thre2)=[];
else %If the first fork is a right fork
    length_eyes=(Dcut(i).threright_fork(2:end)-Dcut(i).threleft_fork(1:length(Dcut(i).threright_fork)-1))+1;
    Dcut(i).threleft_fork(length_eyes<thre2)=[];
    Dcut(i).threright_fork(find(length_eyes<thre2)+1)=[];
    end

end

end
%>----------------------------------------------------------------------<


%I calculate gap length (between eyes) considering also eyes at side
for i=1:num_pieces
%If there is at least an eye...
if ~isempty(Dcut(i).threleft_fork) && ~isempty(Dcut(i).threright_fork) 
if Dcut(i).threleft_fork(1)<=Dcut(i).threright_fork(1)  
    Dcut(i).gap_length=(Dcut(i).threleft_fork(2:end)-Dcut(i).threright_fork(1:(length(Dcut(i).threleft_fork)-1)))-1;
else 
    Dcut(i).gap_length=(Dcut(i).threleft_fork(1:end)-Dcut(i).threright_fork(1:length(Dcut(i).threleft_fork)))-1;
end
%If there are no eyes, "gap_length" is an empty vector
else
    Dcut(i).gap_length=[];
end
end



%I eliminate eyes at the side
for i=1:num_pieces
%I copy the forks in ".threright_forkwithside" and ".threleft_forkwithside"
%and I eliminate the forks of the eyes at side from ".threright_fork" and ".threleft_fork"
    Dcut(i).threright_forkwithside=Dcut(i).threright_fork;
    Dcut(i).threleft_forkwithside=Dcut(i).threleft_fork;
%If there is just a right fork, is a fork at the side so I eliminate it.
if ~isempty(Dcut(i).threright_fork) && isempty(Dcut(i).threleft_fork)
   Dcut(i).threright_fork(1)=[];
%If there is at least an eye and the fist right fork is before the fist
%left fork, I eliminate the fist right fork.
elseif ~isempty(Dcut(i).threright_fork) && ~isempty(Dcut(i).threleft_fork)
    if  Dcut(i).threright_fork(1)<Dcut(i).threleft_fork(1)
    Dcut(i).threright_fork(1)=[];
    end
end
%If after the fist two checks, the number of left and right fork is not the
%same, means that there is a left fork at side and I eliminate it.
if length(Dcut(i).threright_fork)~=length(Dcut(i).threleft_fork)
    Dcut(i).threleft_fork(length(Dcut(i).threleft_fork))=[];
end
end


%I calculate length of eyes
for i=1:num_pieces
if ~isempty(Dcut(i).threright_fork) && ~isempty(Dcut(i).threleft_fork)
Dcut(i).length_eyes=(Dcut(i).threright_fork(1:end)-Dcut(i).threleft_fork(1:end))+1;
else
 Dcut(i).length_eyes=[];
end
end


%I calculate centre of eyes
for i=1:num_pieces
if ~isempty(Dcut(i).threright_fork) && ~isempty(Dcut(i).threleft_fork)
Dcut(i).centre_eyes=(Dcut(i).threright_fork(1:end)+Dcut(i).threleft_fork(1:end))/2;
else
  Dcut(i).centre_eyes=[];
end
end


%I calculate eye to eye distances and distances between neighborhing eyes
for i=1:num_pieces
Dcut(i).etedist=Dcut(i).centre_eyes(2:2:end)-Dcut(i).centre_eyes(1:2:end-1);
end
end

%I calculate gap length(between eyes) without eyes at side
% for i=1:num_pieces
% Dcut(i).gap_length=(Dcut(i).threleft_fork(2:end)-Dcut(i).threright_fork(1:end-1));
% end