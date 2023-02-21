function intensities=intensities_extraction(path)

 %Fiber intensities{file}{fiber}{intensities=2}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Extract intensities%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(path)
sprintf('Treatment of file n. %i',j)
fileID = fopen([path{j} 'Log.txt']);
C = textscan(fileID,'%s');
fclose(fileID);
for i=1:length(C{1})
    fileID = fopen([path{j} C{1}{i} '.xls']);
     C_text = textscan(fileID,'%s',2);
    intensities{j}{i} = textscan(fileID,'%f %f'); %{file}{fiber}
    fclose(fileID);
end
end
save('intensities.mat','intensities')

end