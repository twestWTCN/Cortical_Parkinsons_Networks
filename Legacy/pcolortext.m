function pcolortext(matrix)
%prepare position and size of textboxes
pos=get(gca,'position');
[rows,cols]=size(matrix);
width=pos(3)/(cols-1);
height =pos(4)/(rows-1);
%create textbox annotations
for i=1:cols-1
    for j=rows-1:-1:1
        if matrix(j,i)~=0
            annotation('textbox',[pos(1)+width*(i-1)-0.01,pos(2)+height*(j-1),width,height], ...
                'string',sprintf('%.2f',matrix(j,i)),'LineStyle','none','HorizontalAlignment','center',...
                'VerticalAlignment','middle');
        end
    end
end
