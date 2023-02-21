function im_paper(name)
pos = get(gcf, 'papersize');
set(findobj(gcf,'type','axes'),'TickDir','in','TickLength',[0.02 0.035],'FontSize',10,'FontName','Arial','XMinorTick','on','box','on');

width=10;
height=5;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])

print('-dtiff','-r600',[name  '.tif'])
print('-djpeg','-r600',[name  '.jpeg'])

end