function im_paper2(name)
pos = get(gcf, 'papersize');

width=15; 
height=3; 

set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])

print('-dtiff','-r600',[name  '.tif'])
print('-djpeg','-r600',[name  '.jpeg'])

end