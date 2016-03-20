function new_img=edge_dilation(img,bmap,distThr)
%set the boundry to selected wide
% --distThr is the wide of segmentation boundary
% 


[X Y Z]=size(img);
 [r c]=find(bmap==1);
numBdy=length(r);
%% make mask
offset=makeTable(X,Y,distThr);
numOffsets=length(offset);
% offset_gen=offset(:,2)*X+offset(:,1);
masked_points_idx=zeros(numBdy*numOffsets,1);
for e=1:numOffsets
%     e
    Pmask=[r+offset(e,1) c+offset(e,2)];
    excluded=logical(Pmask(:,1)<1|Pmask(:,1)>X)|...
        logical(Pmask(:,2)<1|Pmask(:,2)>Y);
    %got the edges from points to Pnext
    Pidx=(Pmask(:,2)-1)*X+Pmask(:,1); % set the outlier points to 0
    Pidx(excluded)=0;
    masked_points_idx((e-1)*numBdy+1:e*numBdy)=Pidx;
end
idx_masked_points=unique(masked_points_idx); 
idx_masked_points(1)=[]; %excluded the outlier points
masked_boundary_img=zeros(X,Y);
masked_boundary_img(idx_masked_points)=1;
% bmap=masked_boundary_img;

new_img=img(:,:,1);
new_img(idx_masked_points)=0;
temp=img(:,:,2);
temp(idx_masked_points)=0;
new_img(:,:,2)=temp;
temp=img(:,:,3);
temp(idx_masked_points)=255;
new_img(:,:,3)=temp;
% figure;imshow(new_img);
% imwrite(new_img,strcat(img_name,'.bmp'));
end


%%%%%%%%%%%%%%%%%%%%%%%function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [offset ]=makeTable(X,Y,r)
%Output:
%       offset-- the offset table from the centre point
%       detour-- the 2-path table
%       e-- the num of edges from centre point
%imageSize=size(img);
r_ceil=ceil(r);
ci=[1+r_ceil 1+r_ceil r];
%% Make table offset[]
[xx,yy] = meshgrid((1:X)-ci(1),(1:Y)-ci(2));
mask = (xx.^2 + yy.^2)<=ci(3)^2;
[row,col]=find(mask);
offset=zeros(size(row,1),2);
tempx=row-1-r_ceil;tempy=col-1-r_ceil;
offset(:,1)=tempx;offset(:,2)=tempy;
offset_gen=offset(:,1)*X+offset(:,2);

end
