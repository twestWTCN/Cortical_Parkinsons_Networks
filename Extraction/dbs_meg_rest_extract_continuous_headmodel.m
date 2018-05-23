function D = dbs_meg_rest_extract_continuous_headmodel(initials, drug)

druglbl = {'off', 'on'};
prefix = '';

keep = 0;
tsss = 1;

try
    [files, seq, root, details] = dbs_subjects(initials, drug);
catch
    D = [];
    return
end

%remove date from directory
ind=strfind(root,'\');
root=root(1:ind(end-1));

%%
if ~exist(root, 'dir')
    if ismember(root(end), {'/', '\'})
        root(end) = [];
    end
    [p, p1] = fileparts(root);
    if ~exist(p, 'dir')
        if ismember(p(end), {'/', '\'})
            p(end) = [];
        end
        [p2, p3] = fileparts(p);
        if ~exist(p2, 'dir')
            if ismember(p2(end), {'/', '\'})
                p2(end) = [];
            end
            [p4, p5] = fileparts(p2);
            res = mkdir(p4, p5);
        end
        res = mkdir(p2, p3);
    end
    res = mkdir(p, p1);
end

cd(root);

% if details.cont_head_loc
%     [alldewar, alldist, allsens, allfid] = dbs_meg_headloc(files);
% else
alldewar = [];
alldist = [];
% end

fD = {};
%%
origsens  = {};
origfid   = [];
fixedsens = {};
fixedfid  = [];
cnt=1;
nf=1;

for f = 1:numel(files)
    if ~isequal('R', seq{f}(1));
        continue
    end
    
    
    % =============  Conversion =============================================
    S = [];
    S.dataset = files{f};
    S.channels = details.chanset;
    S.checkboundary = 0;
    S.saveorigheader = 1;
    S.conditionlabels = seq{f};
    S.mode = 'continuous';
    S.outfile = ['spmeeg' num2str(f) '_' spm_file(S.dataset,'basename')];
    
    
    if details.brainamp
        if details.hamburg
            S.ref2 = 'UADC004';
        end
        
        D = dbs_meg_brainamp_preproc(S);
    else
        D = spm_eeg_convert(S);
    end
    
    % =============  Post conversion refinements
    
    D = sensors(D, 'EEG', []);
    
    origsens =[origsens  {ft_datatype_sens(ft_read_sens(S.dataset, 'senstype', 'meg'), 'amplitude', 'T', 'distance', 'mm')}];
    
    if details.berlin
        origfid  =[origfid ...
            ft_convert_units(ft_read_headshape(spm_file(S.dataset, 'filename', details.markers(f).files{1})))];
        
        D = fiducials(D, origfid(end));
        
        hdr = ft_read_header(S.dataset);
        
        event     = ft_read_event(S.dataset, 'detectflank', 'down', 'trigindx', ...
            spm_match_str({'EEG157', 'EEG159', 'EEG160'}, hdr.label), 'threshold', 5e-3);
        
        save(D);
    else
        origfid  =[origfid ft_convert_units(ft_read_headshape(S.dataset), 'mm')];
        event    = ft_read_event(fullfile(D), 'detectflank', 'both');
    end
    
    
    eventdata = zeros(1, D.nsamples);
    trigchanind = D.indchannel(details.eventchannel);
    
    if ~isempty(event)
        trigind  = find(strcmp(details.eventtype, {event.type}));
        eventdata([event(trigind).sample]) = 1;
    end
    
    D = chanlabels(D, trigchanind, 'event');
    
    D(D.indchannel('event'), :) = eventdata;
    
    save(D);
    
    if details.cont_head_loc
        S = [];
        S.D = D;
        if ~isempty(alldewar)
            S.valid_fid = squeeze(alldewar(:, :, f));
        end
        S.mode = 'mark';
        D = spm_eeg_fix_ctf_headloc(S);
    end
    
    fixedsens = [fixedsens {sensors(D, 'MEG')}];
    fixedfid  = [fixedfid fiducials(D)];
    
    % Montage ========== =============================================
    if ~isempty(details.montage)
        montage = details.montage;
        
        S = [];
        S.D = D;
        
        if details.oxford
            ind = [];
            
            neeg = length(strmatch('EEG', D.chanlabels));
            if neeg < 8
                ind = [ind; strmatch('Cz', montage.labelnew, 'exact')];
            end
            
            if neeg < 5
                ind = [ind; strmatch('HEOG', montage.labelnew, 'exact')];
                ind = [ind; strmatch('VEOG', montage.labelnew, 'exact')];
            end
            
        elseif details.berlin
        elseif details.hamburg
        else
            ind = [];
            if isempty(strmatch('HLC', D.chanlabels))
                ind = [ind; strmatch('HLC', montage.labelnew)];
            end
            
            neeg = length(strmatch('EEG', D.chanlabels));
            if neeg < 16
                ind = [ind; strmatch('Cz', montage.labelnew, 'exact')];
            end
            
            if neeg < 13
                ind = [ind; strmatch('HEOG', montage.labelnew, 'exact')];
                ind = [ind; strmatch('VEOG', montage.labelnew, 'exact')];
            end
            
            if neeg == 0
                ind = [ind; strmatch('EMG', montage.labelnew)];
                ind = [ind; strmatch('LFP', montage.labelnew)];
            end
            
            if ~isempty(ind)
                montage.labelnew(ind) = [];
                montage.tra(ind, :) = [];
            end
        end
        
        S.montage = montage;
        S.keepothers = 0;
        D = spm_eeg_montage(S);
        
        if ~keep, delete(S.D);  end
        
    end
    
    D = chantype(D, D.indchannel(details.chan), 'LFP');
    
    if details.hamburg
        D = chantype(D, D.indchannel(details.eegchan), 'EEG');
    end
    
    save(D);
    
    % Interpolate Missing Data =========================================
    a = D(:,:,:);
    for x = 1:size(a,1)
        if sum(isnan(a(x,:)))> 1
            warning([D.chanlabels{x} ' has ' num2str(sum(isnan(a(x,:)))) ' missing values, interpolating'])
            [F,TF] = fillmissing( a(x,:),'linear','SamplePoints',1:length(a(x,:)));
            a(x,:) = F;
        end
    end
    D(:,:,:) = a;
    clear a
    save(D);
    
    % Downsample =======================================================
    S = [];
    S.D = D; S.fsample_new = 128;
    D = spm_eeg_downsample(S);
    
    % Artefact marking in MEG =========================================
    % only detect and mark flat MEG channels & jumps
    
    if details.berlin
        [ampthresh, flatthresh] = berlin_gain(spm_file(files{f}, 'ext', 'gain.txt'));
    end
    
    S = [];
    S.D = D;
    S.mode = 'mark';
    S.badchanthresh = 0.8;
    if details.oxford
        S.methods(1).channels = {'MEGMAG'};
    else
        S.methods(1).channels = {'MEGGRAD'};
    end
    S.methods(1).fun = 'flat';
    if details.berlin
        S.methods(1).settings.threshold = flatthresh;%***** Flat thresh
    else
        S.methods(1).settings.threshold = 1e-010;
    end
    S.methods(1).settings.seqlength = 10;
    S.methods(2).channels = {'MEGPLANAR'};
    S.methods(2).fun = 'flat';
    S.methods(2).settings.threshold = 0.1;
    S.methods(2).settings.seqlength = 10;
    if details.oxford
        S.methods(3).channels = {'MEGMAG'};
    else
        S.methods(3).channels = {'MEGGRAD'};
    end
    S.methods(3).fun = 'jump';
    if details.oxford
        S.methods(3).settings.threshold = 50000;
    elseif details.berlin
        S.methods(3).settings.threshold = 1e4; %****Jump thresh
    else
        S.methods(3).settings.threshold = 20000;
    end
    S.methods(3).settings.excwin = 2000;
    S.methods(4).channels = {'MEGPLANAR'};
    S.methods(4).fun = 'jump';
    S.methods(4).settings.threshold = 5000;
    S.methods(4).settings.excwin = 2000;
    
    if details.oxford
        %S.methods = S.methods(2:4); % Just check MEGPLANAR for now
    elseif details.berlin
        S.methods(5).channels = {'MEG'};
        S.methods(5).fun = 'threshchan';
        S.methods(5).settings.threshold = ampthresh;
        S.methods(5).settings.excwin = 1000;
    end
    
    D = spm_eeg_artefact(S);
    
    if ~keep
        delete(S.D);
    end
    
    if details.oxford && tsss
        S = [];
        S.D = D;
        S.magscale   = 100;
        S.tsss       = 1;
        S.xspace     = 0;
        D = tsss_spm_enm(S);
        
        D = badchannels(D, D.indchantype('MEGANY'), 0);save(D);
        
        if ~keep, delete(S.D);  end
        
        prefix = 'tsss';
    end
    
    % Filtering =========================================
    S = [];
    S.D = D;
    S.band = 'high';
    S.freq = 4;
    S.type = 'butterworth';
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    
    
    save(D);
    % Truncate
    S = [];
    S.D = D;
    S.timewin = [1000 1000.*(D.time(end)-1)];
    S.prefix = 'cr';
    D = spm_eeg_crop(S);
    save(D);
    % Epoching =========================================
    % removed no epoching performed
    
    % Trial rejection =========================================
    % removed no trial rejection
    
    %     save(D);
    
    
    %%
    
    % =====================================================
    % removed merging of same type blocks
    
    % ====================================================
    % compute headmodel - for source extraction later on
    
    %%
    sens = spm_cat_struct(origsens{:}, fixedsens{:});
    fid  = spm_cat_struct(origfid, fixedfid);
    %%
    %         pnt = [];
    %         for i = 1:numel(fid)
    %             pnt = cat(3, pnt, fid(i).fid.pnt);
    %         end
    %% DONT CHECK VALID SENSORS
    valid = logical(ones(1,numel(fid)));
    %         cont_fid  = permute(pnt, [3 1 2]);
    %
    %         dist = [
    %             squeeze(sqrt(sum((cont_fid(:, 1, :) - cont_fid(:, 2, :)).^2, 3)))';...
    %             squeeze(sqrt(sum((cont_fid(:, 2, :) - cont_fid(:, 3, :)).^2, 3)))';...
    %             squeeze(sqrt(sum((cont_fid(:, 3, :) - cont_fid(:, 1, :)).^2, 3)))' ];
    %
    %         rdist = 10*round([dist alldist]/10);
    %
    %         valid = find(all(abs(dist - repmat(mode(rdist, 2), 1, size(dist, 2)))<10));
    %
    %         if isempty(valid)
    %             if ~isempty(alldist)
    %                 sens = allsens;
    %                 fid  = allfid;
    %
    %                 valid = 1:numel(allsens);
    %             else
    %                 error('No valid senosrs found');
    %             end
    %         end
    
    spm_figure('GetWin','Graphics');clf;
    h = axes;
    [asens afid] = ft_average_sens(sens(valid), 'fiducials', fid(valid), 'feedback', h);
    %%
    D = sensors(D, 'MEG', asens);
    D = fiducials(D, afid);
    
    for i = 1:nf
        if ismember(nf + i, valid)
            D.allsens(i) = sens(nf + i);
            D.allfid(i)  = fid(nf + i);
        elseif valid(i)
            D.allsens(i) = sens(i);
            D.allfid(i)  = fid(i);
        else
            D.allsens(i) = asens;
            D.allfid(i)  = afid;
        end
    end
    %%
    D.initials = initials;
    
    if details.oxford
        fids = D.fiducials;
        [~, sel] = spm_match_str({'Nasion', 'LPA', 'RPA'}, fids.fid.label);
        fids.fid.pnt = fids.fid.pnt(sel, :);
        fids.fid.label = {'nas'; 'lpa'; 'rpa'};
        
        D = fiducials(D, fids);
    end
    
    save(D);
    
    D = dbs_meg_headmodelling(D);
    %%
    D = D.move([prefix initials '_' seq{f}(1) '_' num2str(cnt) '_' druglbl{drug+1}]);
    
    cnt=cnt+1;
end
end
