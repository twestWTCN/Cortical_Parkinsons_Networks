function [mask] = makeSphereMask(sDim,sCen,rad,plotop)
if nargin<4
    plotop = 0;
end
% Create a logical image volume of a sphere with specified
% diameter, center, and image size.
% First create the image.
imageSizeX = sDim(1);
imageSizeY = sDim(2);
imageSizeZ = sDim(3);
[columnsInImage rowsInImage pagesInImage] = meshgrid(1:imageSizeX, 1:imageSizeY,1:imageSizeZ);

% Next create the sphere in the image.
centerX = sCen(1);
centerY = sCen(2);
centerZ = sCen(3);
radius = rad;
mask = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 + (pagesInImage - centerZ).^2 <= radius.^2;

% sphereVoxels is a 3D "logical" array.
% Now, display it using an isosurface and a patch
if plotop == 1
    hold on
fv = isosurface(mask,0);
patch(fv,'FaceColor',[0 0 .7],'EdgeColor',[0 0 1],'FaceAlpha',0.2);
view(45,45);
axis equal;
title('Binary volume of a sphere');
end