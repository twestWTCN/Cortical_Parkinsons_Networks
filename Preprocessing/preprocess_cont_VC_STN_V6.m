function preprocess_cont_VC_STN_V6(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
for sub = 1:length(R.subname)
    for cond = 1:length(R.condname)
        for breg = 1:length(R.bregname)
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\VC_new_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'VC_new')
%                 load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                
                % find ipsi channels
                reflist = find(strncmp(VC_new.label,['STN_' R.siden{side}(1)],5));
                ctxlist = find(strncmp(VC_new.label,['ipsi'],4));
                clear x2data x3data x4data
                xdata = VC_new.trial{1}([ctxlist reflist],:);
                fsamp = VC_new.fsample;
                for i = 1:size(xdata,1)
                    x = xdata(i,:);
                    x = x(floor(fsamp*0.75):end-floor(fsamp*0.75));
                    x = (x-mean(x));%./std(x);
                    x2data(i,:) = x;
                end
                x3data = ft_preproc_bandpassfilter(x2data, fsamp,R.pp.cont.thin.bp, [], 'but', 'twopass', 'reduce');
                for i = 1:size(x3data,1)
                    x = x3data(i,:);
                    x = (x-mean(x)); %./std(x);
                    x4data(i,:) = x;
                end
                vc_clean.label = {VC_new.label{[ctxlist reflist]}};
                vc_clean.fsample = VC_new.fsample;
                vc_clean.trial{1} = x4data;
                vc_clean.time{1} = linspace(0,length(x4data)./fsamp,length(x4data));
                vc_clean.history = [VC_new.history{:} {[mfilename '_' date]}];
                mkdir([R.datapathr R.subname{sub} '\ftdata\cleaned'])
                save([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                clear vc_clean
            end
        end
    end
end

