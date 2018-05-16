f = figure;
N = 1;
implotter = @(N) imshow(M(N).cdata);
% setoptions(h,'XLim',[0,10],'YLim',[0,2]);
b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',N, 'min',1, 'max',70);
bgcolor = f.Color;
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                'String','1','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                'String','70','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                'String','Damping Ratio','BackgroundColor',bgcolor);
            
b.Callback =  @(N,M) imshow(M(N).cdata);