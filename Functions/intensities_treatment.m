function [globalallexDcut,globalallnum_pieces,globalalllength_pieces]=intensities_treatment(intensities,file,unit,thre1,thre2,thre_int,Convmicro_kb,Conv_Zeiss100)

%Convert pixels in replicated or not

for j=1:length(file)
    intensities2=[];
    for i=1:length(intensities{j})
        intensities2=intensities{j}{i}{2};
        intensities2(intensities{j}{i}{2}>=thre_int(j))=1;
        intensities2(intensities{j}{i}{2}<thre_int(j))=0;
        fibers_pixel{j}{i}{2}=intensities2;
    end
end
% I remove gaps at side that are smaller than 3.12pixel and
%I combine eyes divided by gaps smaller then thre1=1kb=3.12pixel
%(1kb/(Convmicro_kb*Conv_Zeiss100)) and 
for j=1:length(file)
    intensities2=[];
    for i=1:length(fibers_pixel{j})
        intensities2=fibers_pixel{j}{i}{2};
        pos=1;
        num=1;
        %Remove first gap if <thre1
        while pos<length(intensities2) && intensities2(pos+1)==intensities2(pos)
                num=num+1;
                pos=pos+1;
        end
         if num<thre1/(unit*Convmicro_kb*Conv_Zeiss100) && intensities2(pos)==0  
             intensities2(pos-num+1:pos)=[];
         end
         pos=length(intensities2);
        num=1;
         %Remove last gap if <thre1
        while pos>1 && intensities2(pos-1)==intensities2(pos)
                num=num+1;
                pos=pos-1;
        end
         if num<thre1/(unit*Convmicro_kb*Conv_Zeiss100) && intensities2(pos)==0  
             intensities2(pos:pos+num-1)=[];
         end
         fibers_pixel_correct{j}{i}{2}=intensities2;
     end
        
    
end


for j=1:length(file)
    intensities2=[];
    for i=1:length(fibers_pixel_correct{j})
        intensities2=fibers_pixel_correct{j}{i}{2};
        pos=1;
        num=1;
        %Since I will never gaps<thre1 at the end because I removed them I
        %don't care that I don't consider the last pixel in this analysis
        %(< instead than <= in the while)
        while pos<length(fibers_pixel_correct{j}{i}{2}) 
            if intensities2(pos+1)==intensities2(pos)
                num=num+1;
                pos=pos+1;
            else
                if num<thre1/(unit*Convmicro_kb*Conv_Zeiss100) && intensities2(pos)==0
                intensities2(pos-num+1:pos)=1;
                end
                pos=pos+1;
                num=1;
            end
        end
        fibers_pixel_correct{j}{i}{2}=intensities2;
    end
end

%I convert eyes smaller than thre2 in unreplicated genome
for j=1:length(file)
    intensities2=[];
    for i=1:length(fibers_pixel_correct{j})
        intensities2=fibers_pixel_correct{j}{i}{2};
        pos=1;
        num=1;
        while pos<length(fibers_pixel_correct{j}{i}{2}) 
            if intensities2(pos+1)==intensities2(pos)
                num=num+1;
                pos=pos+1;
            else
                if num<thre2/(unit*Convmicro_kb*Conv_Zeiss100) && intensities2(pos)==1
                intensities2(pos-num+1:pos)=0;
                end
                pos=pos+1;
                num=1;
            end
        end
        %Fo eyes at the end smallest than thre2
        if num<thre2/(unit*Convmicro_kb*Conv_Zeiss100) && intensities2(pos)==1
            intensities2(pos-num+1:pos)=0;       
        end
        fibers_pixel_correct{j}{i}{2}=intensities2;
    end
end


%Convert fibers in pixels in fibers in kb
for j=1:length(file)
    intensities2=[];
    for i=1:length(fibers_pixel_correct{j})
        intensities3=[];
        intensities2=fibers_pixel_correct{j}{i}{2};
        pos=1;
        num=1;
        while pos<length(fibers_pixel_correct{j}{i}{2})
            if intensities2(pos+1)==intensities2(pos)
                num=num+1;
                pos=pos+1;
            else
                intensities3=[intensities3;intensities2(pos)*ones(round(num*Convmicro_kb*Conv_Zeiss100),1)];
                pos=pos+1;
                num=1;
            end
        end
        %Convertion of last sequence of 0/1
        if intensities2(pos)==intensities2(pos-1)
           intensities3=[intensities3;intensities2(pos)*ones(round(num*Convmicro_kb*Conv_Zeiss100),1)];
        else
           intensities3=[intensities3;intensities2(pos-1)*ones(round(num*Convmicro_kb*Conv_Zeiss100),1)];
           intensities3=[intensities3;intensities2(pos)*ones(round(1*Convmicro_kb*Conv_Zeiss100),1)];
        end
        fibers_kb{j}{i}{2}=intensities3;
    end
end

%I prepare structure for the analysis
for j=1:length(file)
    length_pieces=[];
    exDcut=[];
    for i=1:length(fibers_kb{j})
        exDcut(i).intensity=intensities{j}{i}{2};
        exDcut(i).fiber=fibers_kb{j}{i}{2};
        length_pieces(i)=length(intensities{j}{i}{2})*Convmicro_kb*Conv_Zeiss100;
    end
    [exDcut.unit_block]=deal(unit);
    globalallexDcut.(['exDcut' file{j}])=exDcut
    globalallnum_pieces.(['num_pieces' file{j}])=length(fibers_kb{j})
    globalalllength_pieces.(['length_pieces' file{j}])=length_pieces
end

save 'globalallexDcut.mat' globalallexDcut;
save 'globalallnum_pieces.mat' globalallnum_pieces;
save 'globalalllength_pieces.mat' globalalllength_pieces;

end