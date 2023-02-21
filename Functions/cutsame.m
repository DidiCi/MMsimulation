function [pos_pieces]=cutsame(num_pieces,length_pieces,length_fiber,unit)
%Extract number of fiber and distribution of lentgh from excel 
%I define the positions on the dna to cut it with the same
%distribution of the experiment
pos_pieces=cell(num_pieces,1);
cut=round(length_pieces*(1000/unit));
%I create a fiber with the numbers from 1 to the length of the genome and
%at each step I put to zero the numbers at the positions where I will cut
%my genome
fant=1:length_fiber;
%Then I really cut the genome
for i=1:num_pieces
    %ceil rounds the elements of A to the nearest integers greater than or
    %equal to A.
    f = ceil(length_fiber.*rand()); 
    while (f+cut(i))>length_fiber || sum(fant(f:f+cut(i)))~=sum(f:f+cut(i))
        f = ceil(length_fiber.*rand()); 
    end
    pos_pieces{i}=fant(f:f+cut(i));
    fant(f:f+cut(i))=0;
end

end