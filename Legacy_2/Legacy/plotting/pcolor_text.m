function pcolor_text(mat,eps)
pos=get(gca,'position');
% mat = squeeze(stat.prob);
amat = zeros(size(mat)+1);
amat(1:end-1,1:end-1) = mat;
[rows,cols]=size(amat);
width=pos(3)/(cols-1);
height =pos(4)/(rows-1);
%create textbox annotations
for i=1:cols-1
    for j=rows-1:-1:1
        if amat(j,i)~=0 & amat(j,i)<eps
        annotation('textbox',[pos(1)+width*(j-1),pos(2)+height*(i-1),width,height], ...
            'string',sprintf('P = %0.3f ',amat(j,i)),'LineStyle','none','HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontWeight','bold','FontSize',8);
        end
    end
end