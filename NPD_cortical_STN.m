function [npdspctrm_out, Hz] = NPD_cortical_STN(vchansave_ON,vchansave_OFF,R,plotop,refid)
close all
for cond = 1:2
    if cond == 1
        vchansave = vchansave_ON;
    else
        vchansave = vchansave_OFF;
    end
    
    for i =1:numel(vchansave)
        data = vchansave(i);
        x = data.trial{1}(1,:);
        x = (x-mean(x))./std(x);
        y = data.trial{1}(refid(cond),:);
        y = (y-mean(y))./std(y);
%         [f13,t13,c13]=sp2a2_R2_mt(x',y',R.pp.cont.full.fs,10,'M1.5');
                [f13,t13,c13]=sp2a2_R2(x',y',R.pp.cont.full.fs,9);
        Hz = f13(:,1);
        %     if plotop ==1
        %     subplot(1,4,1);plot(Hz,f13(:,4),'k'); hold on; xlim([2 48]); ylim([0 0.1]);
        %     subplot(1,4,2);plot(Hz,f13(:,10),'r'); hold on; xlim([2 48]); ylim([0 0.1]);
        %     subplot(1,4,3);plot(Hz,f13(:,11),'b'); hold on; xlim([2 48]); ylim([0 0.1]);
        %     subplot(1,4,4);plot(Hz,f13(:,12),'g'); hold on; xlim([2 48]); ylim([0 0.1]);
        %     end
        npdspctrm{cond,1,1}(:,i) = f13(:,10);
        npdspctrm{cond,1,2}(:,i) = f13(:,11);
        npdspctrm{cond,1,3}(:,i) = f13(:,12);
        npdspctrm{cond,1,4}(:,i) = f13(:,4);
        
        %  f column 1       frequency in Hz
        %  f column 2       Log input/x  spectrum
        %  f column 3       Log output/y spectrum
        %  f column 4       Coherence
        %  f column 5       Phase
        %  f column 6       MMSE whitened spectra for process x, JIN (2.4)
        %  f column 7       MMSE whitened spectra for process y, JIN (2.4)
        %  f column 8       Coherence between whitened processes, JIN (2.5)
        %  f column 9       Phase between whitened processes
        %  f column 10      Zero lag coherence component, JIN (2.19)
        %  f column 11      Forward  coherence component, JIN (2.20)
        %  f column 12      Reverse  coherence component, JIN (2.18)
    end
end
%zero mean
if plotop ==1
    plotNPD(Hz,npdspctrm,R)
end
for cond = 1:2
    npdspctrm_out{cond,1,1}(:,1) = mean(npdspctrm{cond,1,1},2);
    npdspctrm_out{cond,1,2}(:,1) = mean(npdspctrm{cond,1,2},2);
    npdspctrm_out{cond,1,3}(:,1) = mean(npdspctrm{cond,1,3},2);
    npdspctrm_out{cond,1,4}(:,1) = mean(npdspctrm{cond,1,4},2);
end
