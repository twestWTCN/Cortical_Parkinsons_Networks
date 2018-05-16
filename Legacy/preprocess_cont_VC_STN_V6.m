function preprocess_cont_VC_STN_V6(R)
sidenrev = {R.siden{end:-1:1}};
for sub = 1:numel(R.subname)
    for cond = 1:2
        for band = [1 3]
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\VC_new_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'VC_new')
                % find ipsi channels
                reflist = find(strncmp(VC_new.label,['STN_' sidenrev{side}(1)],5));
                
                clear x2data x3data x4data
                xdata = VC_new.trial{1}([1 reflist],:);
                fsamp = VC_new.fsample;
                for i = 1:size(xdata,1)
                    x = xdata(i,:);
                    x = x(floor(fsamp*0.75):end-floor(fsamp*0.75));
                    x = (x-mean(x))/std(x);
                    x2data(i,:) = x;
                end
                x3data = ft_preproc_bandpassfilter(x2data, fsamp,R.pp.cont.thin.bp, [], 'but', 'twopass', 'reduce');
                for i = 1:size(x3data,1)
                    x = x3data(i,:);
                    x = (x-mean(x))./std(x);
                    x4data(i,:) = x;
                end
                
                vc_clean.label = {VC_new.label{[1 reflist]}};
                vc_clean.fsample = VC_new.fsample;
                vc_clean.trial{1} = x4data;
                vc_clean.time{1} = linspace(0,length(x4data)./fsamp,length(x4data));
                vc_clean.history = [VC_new.history{:} {[mfilename '_' date]}];
                save([R.datapathr R.subname{sub} '\ftdata\virtualV6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'vc_clean')
                clear vc_clean
            end
        end
    end
end

