close all
addpath('C:\Users\twest\Documents\Work\MATLAB ADDONS\linspecer')
% Plot Results
cmap  = linspecer(6,'qualitative');
cmap = cmap(end:-1:1,:)

% Plot Fractional Occupancies across Subjects
figure(1)
[~,idx] = sort(FO(:,6));
sFO = FO(idx,:);
sFO = FO;
b = bar(sFO,'stacked');
condind = (idx<27)*0.5;
condind(condind<0.5) = NaN;
hold on;
sc = scatter(1:size(FO,1),condind)
sc.Marker = 'x'; sc.LineWidth = 2; sc.CData = [0 0 0 ];
for st = 1:no_states
    b(st).FaceColor = cmap(st,:);
end
xlabel('Subject')
ylabel('Fractional Occupancy')
legname = strcat({'State: '},int2str((1:6).')).';
legname{7} = 'OFF Recordings';
legend(legname,'Location','SouthEast');ylim([0 1])
set(gcf,'Position',[250         280        1236         590])

% Plot State Spectra
figure(2)
sti = 0;
for st = 1:no_states
    sti = sti +1;
    subplot(2,6,sti)
    a = spectra.state(st).psd;
    x1 = abs(a(:,1,1));
    x2 = abs(a(:,2,2));
    cx1 = spectra.state(st).coh(:,1,2);
    cx2 = spectra.state(st).coh(:,2,1);
    
    f = spectra.state(st).f;
    plot(f,x1,'color',cmap(st,:),'LineWidth',2); hold on
    plot(f,x2,'color',cmap(st,:).*0.5,'LineWidth',2)
    set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log')
    xlim([2 40]); ylim([0 0.1])
    grid on
    xlabel('Frequency'); ylabel('PSD')
    title(sprintf('State %.0f',st))
    legend('SMA','STN')
end

% Plot State Occupancy Series
figure(2)
clear a
sti = 0;
for st = 1:no_states
    sti = sti +1;
    a(st) = subplot(2,6,sti+6)
    x = Gamma(1:T{1},st);
    t = linspace(0,T{1}/fsamp,T{1});
    plot(t,x,'color',cmap(st,:),'LineWidth',2);
    xlabel('Time (s)');
    ylabel('P(x)');
    grid on
    title(sprintf('State %.0f',st))
end
set(gcf,'Position',[250         280        1236         490])
linkaxes(a,'xy')
xlim([80 90]); ylim([0 1])

% Plot State Coherences
figure(40)
clear a
sti = 0;
for st = 1:no_states
    sti = sti +1;
    subplot(2,6,sti)
    a = spectra.state(st).psd;
    x1 = abs(a(:,1,1));
    x2 = abs(a(:,2,2));
    cx1 = spectra.state(st).coh(:,1,2);
    cx2 = spectra.state(st).coh(:,2,1);
    
    f = spectra.state(st).f;
    plot(f,cx1,'color',cmap(st,:),'LineWidth',2); hold on
    plot(f,cx2,'color',cmap(st,:).*0.5,'LineWidth',2)
    
    xlim([2 40]); ylim([0 1])
    grid on
    xlabel('Frequency'); ylabel('Coherence')
    title(sprintf('State %.0f',st))
    legend('SMA','STN')
end
set(gcf,'Position',[250         280        1236         490])

figure(20)
clear a
sti = 0;
for st = 1:no_states
    x = Gamma(1:T{1},st);
    t = linspace(0,T{1}/fsamp,T{1});
    plot(t,x,'color',cmap(st,:),'LineWidth',2);
    hold on
    xlabel('Time (s)');
    ylabel('P(x)');
    grid on
    title(sprintf('State %.0f',st))
end
set(gcf,'Position',[250         280        1236         490])
xlim([100 110]); ylim([0 1])

figure(3)
% A = zeros(no_states+1)
% A(1:end-1,1:end-1) = TP'
% pcolor(A)
G  =digraph(TP);
colormap(cmap)
plot(G,'EdgeLabel',G.Edges.Weight,'NodeCData',1:size(G.Nodes,1),'MarkerSize',10,'EdgeColor','k')
set(gca,'Visible','off')
set(gcf,'color','white')

cmap = linspecer(2);
figure(4)
for st = 1:no_states
    subplot(1,6,st)
    FO_OFF = FO(1:26,st);
    FO_ON = FO(26:52,st);
    histogram(FO_OFF,0:.1:1,'Normalization','probability','FaceColor',cmap(1,:),'FaceAlpha',0.9);
    hold on
    histogram(FO_ON,0:.1:1,'Normalization','probability','FaceColor',cmap(2,:),'FaceAlpha',0.9);
    p =ranksum(FO_OFF,FO_ON);
    if p<0.05
        text(0.05,0.4,sprintf('RS Test: P = %.3f',p)); ylim([0 1])
    end
    xlabel('Fractional Occupancy')
    ylabel('P(X)');
    grid on
    title(sprintf('State %.0f',st))
    if st == no_states
        legend(R.condname)
    end
end
set(gcf,'Position',[6         856        1898         242])

figure(5)
clear x
for st = 1:no_states
    subplot(1,6,st)
    binEdge = logspace(log10(0.01),log10(15),12);
    binMid = binEdge(1:end-1)+((binEdge(2)-binEdge(1))/2);
    x(1,:) = histcounts((LifeTimes{20,st}/fsamp),binEdge);
    x(1,:) = (x(1,:)./(T{20}/fsamp)')*60;
    x(2,:) = histcounts((LifeTimes{46,st}/fsamp),binEdge);
    x(2,:) = (x(2,:)./(T{46}/fsamp)')*60;
    b = bar(log10(binMid),x',1.2); hold on
    %     b(2) = bar(log10(binMid),x(2,:),1); hold on
    b(1).FaceColor = cmap(1,:);
    b(2).FaceColor = cmap(2,:);
    hold on; grid on
    ylim([0 8]); ylabel('Occurence Rate (min^-1)')
    xlabel('Life Time (log s)')
    title(sprintf('State %.0f',st))
    if st == no_states
        legend(R.condname)
    end
end
set(gcf,'Position',[6         856        1898         242])

R.obs.logdetrend = 0;
figure(25)
R.frqz = 4:48; R.condnames{1} = '1'; R.obs.obsstates = [1 2];
sti = 0;
for cond = 1:2
    if cond == 1;
        Ts = get_state_tseries(data,Scor,idx(1:10)');
    else
        Ts = get_state_tseries(data,Scor,idx(end-10:end)');
    end
    for st = 1:size(Ts,2)
        sti = sti+1;
        subplot(2,size(Ts,2),sti)
        if ~isempty(Ts{st})
            [F meannpd] = constructNPDMat_190618({Ts{st}'},{'SMA','STN'},{'SMA','STN'},Hz,8,R);
            plot(F,squeeze(meannpd(:,1,2,2,:)),'r'); hold on
            plot(F,squeeze(meannpd(:,1,2,3,:)),'b')
        end
    end
end
set(gcf,'Position',[250         280        1236         490])

