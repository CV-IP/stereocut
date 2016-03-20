function  StereoSaliencyCutMain

string_i = '008' ;
    
%input image pair
name_img_l = strcat('imgs\\bit',string_i,'_l.jpg') ;
name_img_r = strcat('imgs\\bit',string_i,'_r.jpg') ;    
    
I_L = imread(name_img_l);
I_R = imread(name_img_r);
 
%input saliency maps
name_img_l = strcat('saliency\\saliency',string_i,'_l.png') ;
name_img_r = strcat('saliency\\saliency',string_i,'_r.png') ;   

sal_L = imread(name_img_l) ;
sal_R = imread(name_img_r) ;

[r , c] = size(sal_L) ;
numVars = r*c ;

sal_L = double(sal_L(:))/255.0 ;
sal_R = double(sal_R(:))/255.0 ;

Colors_L = reshape(I_L,[r*c,3]);
Colors_R = reshape(I_R,[r*c,3]);

Colors_L = double(Colors_L) ; 
Colors_R = double(Colors_R) ; 

name_img_l = strcat('disp\\bit',string_i,'_dp.jpg') ; 
disp = imread(name_img_l) ;
temp = min(disp(:)) ;
disp = int32(disp/temp) ;       
disp = double(disp) ;

%%%%%%%%%%%%%%%%data term£ºUc%%%%%%%%%%%%%%%%%%%%%%%
%first, get the histograms with weights on l/a/b channel
%then, normalize them
Fhist_r = get_weighted_hist(Colors_L(:,1) , sal_L , Colors_R(:,1) , sal_R ) ;
Fhist_g = get_weighted_hist(Colors_L(:,2) , sal_L , Colors_R(:,2) , sal_R ) ;
Fhist_b = get_weighted_hist(Colors_L(:,3) , sal_L , Colors_R(:,3) , sal_R ) ;

I = ones(numVars,1) ;
Bhist_r = get_weighted_hist(Colors_L(:,1) , I-sal_L , Colors_R(:,1) , I-sal_R ) ;
Bhist_g = get_weighted_hist(Colors_L(:,2) , I-sal_L , Colors_R(:,2) , I-sal_R ) ;
Bhist_b = get_weighted_hist(Colors_L(:,3) , I-sal_L , Colors_R(:,3) , I-sal_R ) ;
%f(sn)=sn

FPr_r_L = Fhist_r(Colors_L(:,1)+1) ;   FPr_g_L = Fhist_g(Colors_L(:,2)+1) ;   FPr_b_L = Fhist_b(Colors_L(:,3)+1) ;
BPr_r_L = Bhist_r(Colors_L(:,1)+1) ;   BPr_g_L = Bhist_g(Colors_L(:,2)+1) ;   BPr_b_L = Bhist_b(Colors_L(:,3)+1) ;
FPr_r_R = Fhist_r(Colors_R(:,1)+1) ;   FPr_g_R = Fhist_g(Colors_R(:,2)+1) ;   FPr_b_R = Fhist_b(Colors_R(:,3)+1) ;
BPr_r_R = Bhist_r(Colors_R(:,1)+1) ;   BPr_g_R = Bhist_g(Colors_R(:,2)+1) ;   BPr_b_R = Bhist_b(Colors_R(:,3)+1) ;

%get Uc
FDist_Uc_L = -log ( FPr_r_L .* FPr_g_L .* FPr_b_L ) ;
BDist_Uc_L = -log ( BPr_r_L .* BPr_g_L .* BPr_b_L ) ;
FDist_Uc_R = -log ( FPr_r_R .* FPr_g_R .* FPr_b_R ) ;
BDist_Uc_R = -log ( BPr_r_R .* BPr_g_R .* BPr_b_R ) ;

clear FPr_r_L FPr_g_L FPr_b_L     BPr_r_L BPr_g_L BPr_b_L   ...
      FPr_r_R FPr_g_R FPr_b_R     BPr_r_R BPr_g_R BPr_b_R  ;
  
  
%%%%%%%%%%%%%%%%%data term£ºUs%%%%%%%%%%%%%%%%%%%%%%%%%
%f(sn)=max(0, sign(sn-¦Ó))    or    f(sn) = (sn)¦Ê 
para_k=1 ;
FDist_Us_L = I - sal_L .^para_k ;
BDist_Us_L = sal_L .^para_k ;
FDist_Us_R = I - sal_R .^para_k ;
BDist_Us_R = sal_R .^para_k ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  The Mex Function for GraphCutSegment %%%%%%%%%%%%
 lambda_s = 4 ;        
 lambda_p = 20 ;  
 lambda_c = 0.000001 ;    
[SegImage_L,SegImage_R] = StereoCut_GMM( r,c , Colors_L, Colors_R , ...
                                         FDist_Uc_L , BDist_Uc_L , FDist_Uc_R , BDist_Uc_R , ...
                                         FDist_Us_L , BDist_Us_L , FDist_Us_R , BDist_Us_R , ...                                         
                                         disp , lambda_s , lambda_p , lambda_c );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%post processing%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SegImage_L,SegImage_R] = post_processing(SegImage_L,SegImage_R,r,c) ;

%%%%%%%%%%%%%%%    Save the segmented image   %%%%%%%%%%%%%%%
SegImage_L = repmat(SegImage_L,[1,1,3]);
SegNewImage_L = uint8(SegImage_L) .* uint8(I_L);

SegImage_R = repmat(SegImage_R,[1,1,3]);
SegNewImage_R = uint8(SegImage_R) .* uint8(I_R);

pad = zeros(r,5,3)+255 ;
result = [SegNewImage_L,pad,SegNewImage_R] ;

name_save = strcat('result\\cut',string_i,'.bmp') ;    
imwrite(result, name_save ,'bmp') ;
 
  
 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Helper Functions declarations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function weighted_hist = get_weighted_hist( channel_L , sal_L , channel_R , sal_R )

numVars = length(channel_L) ;

weighted_hist = zeros(256,1);

for i=1:numVars
    ind = channel_L(i)+1 ;
    weighted_hist(ind) = weighted_hist(ind) + sal_L(i) ;
    ind = channel_R(i)+1 ;
    weighted_hist(ind) = weighted_hist(ind) + sal_R(i) ;    
end

weighted_hist = weighted_hist ./ sum(weighted_hist) ;


        
