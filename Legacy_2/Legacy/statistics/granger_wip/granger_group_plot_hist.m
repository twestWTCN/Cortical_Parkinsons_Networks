x1 = []; x2 = [];
for sub = 1:numel(subname)
    for side = 1:2
        W = gc_dist_sub_save{sub,side};
        x1 = [x1 W{1}];
        x2 = [x2 W{2}];
    end
end
[h p] =ttest2(x1,x2)
figure(1)
GH(1) = histogram(x1,-1:0.025:1,'Normalization','probability'); hold on
GH(2) = histogram(x2,-1:0.025:1,'Normalization','probability');
xlabel('GC ')

x1 = []; x2 = [];
for sub = 1:numel(subname)
    for side = 1:2
        W = gc_cd_dist_sub_save{sub,side};
        x1 = [x1 W{1}];
        x2 = [x2 W{2}];
    end
end
[h p] =ttest2(x1,x2)
figure(2)
GH(1) = histogram(x1,0:0.025:.5,'Normalization','probability'); hold on
GH(2) = histogram(x2,0:0.025:.5,'Normalization','probability');
xlabel('GC - CD')

% amplitude
x1 = []; x2 = [];
for sub = 1:numel(subname)
    for side = 1:2
        W =  amp_dist_fs{1,side,sub};
        x1 = [x1 W];
        W =  amp_dist_fs{2,side,sub};
        x2 = [x2 W];
    end
end
bank1 = x1;
bank2 = x2;

% Granger
x1 = []; x2 = [];
for sub = 1:numel(subname)
    for side = 1:2
        W =  gc_dist_c_fs{1,side,sub};
        x1 = [x1 W];
        W =  gc_dist_c_fs{2,side,sub};
        x2 = [x2 W];
    end
end
bank1 = [bank1; x1];
bank2 = [bank2; x2];

% SegL
x1 = []; x2 = [];
for sub = 1:numel(subname)
    for side = 1:2
        W =  segL_dist_fs{1,side,sub};
        x1 = [x1 W];
        W =  segL_dist_fs{2,side,sub};
        x2 = [x2 W];
    end
end
% bank1 = [bank1; x1];
% bank2 = [bank2; x2];

X = bank2;
for i = 1:size(X,1)
    a = X(i,:);
    X(:,isnan(a)) = [];
%     if i == 6
%         X(:,X(i,:)<0.05) = [];
%     end
end