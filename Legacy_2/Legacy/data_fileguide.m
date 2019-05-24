function [filname,pp_mark,nrep,senscheck] = data_fileguide(sub,OFF)
pp_mark = []; senscheck{1} = 1; senscheck{2} = 1;
switch sub
    case 'CP'
        if OFF == 1
            filname{1} = [];
            disp(['No Data available for ' sub ' OFF ' num2str(OFF)])
        elseif OFF == 0
            filname{1} = 'cp021245a_vladimirDBS_20090130_03.ds';
            pp_mark  = 'r_';
        end
        
    case 'DF'
        if OFF == 1
            filname{1} = 'df140551_vladimirDBS_20090807_08.ds';
        elseif OFF == 0
            filname{1} = 'df140551_vladimirDBS_20090806_08.ds';
            pp_mark  = 'r_';
        end
        
    case 'DP'
        if OFF == 1
            filname{1} = 'dp190452_vladimirDBS_20090828_07.ds';
            senscheck{1} = 0;
        elseif OFF == 0
            filname{1} = 'dp190452_vladimirDBS_20090827_07.ds';
            pp_mark  = 'r_';
        end
        
    case 'DS'
        if OFF == 1
            filname{1} = 'ds150148_vladimirDBS_20080225_07.ds';
            senscheck{1} = 0;
        elseif OFF == 0
            filname{1} = 'ds150148_vladimirDBS_20080226_04.ds';
            pp_mark  = 'r_';
        end
        
    case 'JA'
        if OFF == 1
            filname{1} = 'ja271247_vladimirDBS_20090507_08.ds';
            senscheck{1} = 0;
        elseif OFF == 0
            filname{1} ='ja271247_vladimirDBS_20090508_03.ds';
            pp_mark  = 'r_';
            senscheck{1} = 0;
        end
        
    case 'JB'
        if OFF == 1
            filname{1} = 'jb060560_vladimirDBS_20081117_07.ds';
        elseif OFF == 0
            filname{1} = 'jb060560_vladimirDBS_20081114_05.ds';
            filname{2} = 'jb060560_vladimirDBS_20081114_15.ds';
        end
        
    case 'JN'
        if OFF == 1
            filname{1} = 'jn270257_vladimirDBS_20090731_01.ds';
        elseif OFF == 0
            filname{1} = 'jn270257_vladimirDBS_20090730_02.ds';
            pp_mark  = 'r_';
        end
        
    case 'JP'
        if OFF == 1
            filname{1} = 'jp180257_vladimirDBS_20080223_08.ds';
        elseif OFF == 0
            filname{1} = 'jp180257_vladimirDBS_20080224_06.ds';
            pp_mark  = 'r_';
        end
        
    case 'JS'
        if OFF == 1
            filname{1} = [];
            disp(['No Data available for ' sub ' OFF ' num2str(OFF)])
        elseif OFF == 0
            filname{1} = 'js161040_vladimirDBS_20070728_03.ds';
            pp_mark  = 'r_';
        end
        
    case 'KB'
        if OFF == 1
            filname{1} = 'kb300763_vladimirDBS_20070919_08.ds';
            senscheck{1} = 0;
            filname{2} = 'kb300763_vladimirDBS_20070919_13.ds';
            senscheck{2} = 0;
            pp_mark  = 'r_';
        elseif OFF == 0
            filname{1} = 'kb300763_vladimirDBS_20070917_01.ds';
            filname{2} = 'kb300763_vladimirDBS_20070917_10.ds';
            pp_mark  = 'r_';
        end
        
    case 'LM'
        if OFF == 1
            filname{1} = 'ls050651_vladimirDBS_20071102_05.ds';
            pp_mark  = 'r_';
            senscheck{1} = 0;
        elseif OFF == 0
            filname{1} = 'ls050651_vladimirDBS_20071103_06.ds';
            pp_mark  = 'r_';
        end
        
    case 'MC'
        if OFF == 1
            filname{1} = 'mc140651_vladimirDBS_20090724_07.ds';
            pp_mark  = 'r_';
        elseif OFF == 0
            filname{1} = 'mc140651_vladimirDBS_20090725_03.ds';
            pp_mark  = 'r_';
        end
        
    case 'MW'
        if OFF == 1
            filname{1} = 'mw161052_vladimirDBS_20070330_01.ds';
            pp_mark  = 'r_';
        elseif OFF == 0
            filname{1} = 'mw161052_vladimirDBS_20070331_05.ds';
            pp_mark  = 'r_';
        end
        
    case 'PB'
        if OFF == 1
            filname{1} = [];
            disp(['No Data available for ' sub ' OFF ' num2str(OFF)])
        elseif OFF == 0
            filname{1} = 'pb280543_vladimirDBS_20080207_01.ds';
            pp_mark  = 'r_';
        end
        
    case 'RP'
        if OFF == 1
            filname{1} = [];
            disp(['No Data available for ' sub ' OFF ' num2str(OFF)])
        elseif OFF == 0
            filname{1} = 'rp191142_vladimirDBS_20080806_06.ds';
            pp_mark  = 'r_';
        end
        
    case 'SW'
        if OFF == 1
            filname{1} = 'sw260568_vladimirDBS_20070219_02.ds'; senscheck{1} = 0;
            filname{2} = 'sw260568_vladimirDBS_20070219_14.ds'; senscheck{2} = 0;
            pp_mark  = 't_';
        elseif OFF == 0
            filname{1} = 'sw260568_vladimirDBS_20070216_01.ds';
            pp_mark  = 't_';
        end
        
    case 'WB'
        if OFF == 1
            filname{1} = '100154_vladimirDBS_20080530_08.ds';
            senscheck{1} = 0;
            pp_mark  = 'r_';
            
        elseif OFF == 0
            filname{1} = 'wb100154_vladimirDBS_20080529_02.ds';
            pp_mark  = 'r_';
        end
    case 'LN01'
          if OFF == 1
            filname{1} = '100154_vladimirDBS_20080530_08.ds';% fake filname
            filname{2} = '100154_vladimirDBS_20080530_08.ds';% fake filname
            pp_mark  = 'r_';
            
        elseif OFF == 0
            filname{1} = 'wb100154_vladimirDBS_20080529_02.ds'; % fake filname
            filname{2} = '100154_vladimirDBS_20080530_08.ds'; % fake filname
            pp_mark  = 'r_';
          end    
    case 'LN02'
          if OFF == 1
              filname{1} = '100154_vladimirDBS_20080530_08.ds';% fake filname
              filname{2} = '100154_vladimirDBS_20080530_08.ds';% fake filname
            pp_mark  = 'r_';
            
        elseif OFF == 0
            filname{1} = 'wb100154_vladimirDBS_20080529_02.ds'; % fake filname
            pp_mark  = 'r_';
          end
    case 'LN03'
          if OFF == 1
            filname{1} = '100154_vladimirDBS_20080530_08.ds';% fake filname
            filname{2} = '100154_vladimirDBS_20080530_08.ds';% fake filname
            pp_mark  = 'r_';
            
        elseif OFF == 0
            filname{1} = 'wb100154_vladimirDBS_20080529_02.ds'; % fake filname
            filname{2} = '100154_vladimirDBS_20080530_08.ds'; % fake filname
            pp_mark  = 'r_';
        end              
end
nrep = numel(filname);





