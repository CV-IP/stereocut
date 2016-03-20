function [SegImage_L,SegImage_R] = post_processing(SegImage_L,SegImage_R,r,c)

 se = strel('ball',7,7);
 
  temp_L = imdilate(SegImage_L,se,'same');
  SegImage_L = imerode(temp_L,se,'same');
 
  temp_R = imdilate(SegImage_R,se,'same');
  SegImage_R = imerode(temp_R,se,'same');
  
  ind = find(SegImage_L>=0.5) ;
  SegImage_L = zeros(r,c) ;
  SegImage_L(ind) =1 ;
  
  ind = find(SegImage_R>=0.5) ;
  SegImage_R = zeros(r,c) ;
  SegImage_R(ind) =1 ;  

  
[arr,re_num] = bwlabel(SegImage_L,4) ;
for u=1:re_num
    ind = find(arr==u);
    num_pixel = length(ind) ;
    
    if num_pixel < 0.01*r*c
        SegImage_L(ind)=0 ;
    end
end

[arr,re_num] = bwlabel((1-SegImage_L),4) ;
for u=1:re_num
    ind = find(arr==u);
    num_pixel = length(ind) ;
    
    if num_pixel < 0.01*r*c
        SegImage_L(ind)=1 ;
    end
end

[arr,re_num] = bwlabel(SegImage_R,4) ;
for u=1:re_num
    ind = find(arr==u);
    num_pixel = length(ind) ;
    
    if num_pixel < 0.01*r*c
        SegImage_R(ind)=0 ;
    end
end
 
[arr,re_num] = bwlabel((1-SegImage_R),4) ;
for u=1:re_num
    ind = find(arr==u);
    num_pixel = length(ind) ;
    
    if num_pixel < 0.01*r*c
        SegImage_R(ind)=1 ;
    end
end

clear arr re_num num_pixel  ind ;