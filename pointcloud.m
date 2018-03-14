% script to make point cloud

[x,y,z] = ndgrid(1:source.dim(1),1:source.dim(2),1:source.dim(3));
x = reshape(x,1,[]);y = reshape(y,1,[]); z = reshape(z,1,[]);
scatter3(x,y,z,reshape(source.avg.coh,[],1)*500,reshape(source.avg.coh,[],1),'filled')