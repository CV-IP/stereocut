function  StereoStrokeCutMain

string_i = '191' ;

%% input

%get input image pair
name_img_l = strcat('imgs\\bit',string_i,'_l.jpg') ;
name_img_r = strcat('imgs\\bit',string_i,'_r.jpg') ;    
    
I_L = imread(name_img_l);
I_R = imread(name_img_r);
 
%get Foreground and Background Pixels
name_img_l = strcat('stroke\\stroke',string_i,'_l.jpg') ;
name_img_r = strcat('stroke\\stroke',string_i,'_r.jpg') ;   
 
stroke_L = rgb2gray( imread(name_img_l) );
stroke_R = rgb2gray( imread(name_img_r) );

[r , c] = size(stroke_L) ;

Colors_L = reshape(I_L,[r*c,3]);
Colors_R = reshape(I_R,[r*c,3]);

%'double' is for the k-means operation later
Colors_L = double(Colors_L) ; 
Colors_R = double(Colors_R) ; 


% get disparity map
name_img_l = strcat('disp\\disp',string_i,'.jpg') ; 
disp = imread(name_img_l) ;
temp = min(disp(:)) ;
disp = int32(disp/temp) ;         
disp = double(disp) ;    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding Foreground and Background Labels from pixel seeds
FLabels_L = find( stroke_L == 255 )  ;
BLabels_L = find( stroke_L == 0 )    ;

FLabels_R = find( stroke_R == 255 )  ;
BLabels_R = find( stroke_R == 0 )    ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finding foregound and background color clusters 
FColors_L = Colors_L(FLabels_L,:); 
BColors_L = Colors_L(BLabels_L,:);

FColors_R = Colors_R(FLabels_R,:);  
BColors_R = Colors_R(BLabels_R,:);

FColors = [FColors_L;FColors_R] ;
BColors = [BColors_L;BColors_R] ;

NumFClusters = 5;
NumBClusters = 5;

%%   GMM  %%
[ FCClusters , FCovs , FWeights ] = get_GMM( FColors, NumFClusters ) ;
[ BCClusters , BCovs , BWeights ] = get_GMM( BColors, NumBClusters ) ;

% Foreground and Background edge weights
FDist_L = ClustDistMembership(Colors_L , FCClusters, FCovs, FWeights);
BDist_L = ClustDistMembership(Colors_L , BCClusters, BCovs, BWeights);

FDist_R = ClustDistMembership(Colors_R , FCClusters, FCovs, FWeights);
BDist_R = ClustDistMembership(Colors_R , BCClusters, BCovs, BWeights);

clear FCClusters  FCovs  FWeights  BCClusters  BCovs  BWeights ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  The Mex Function for GraphCutSegment %%%%%%%%%%%%
FLabels_L = FLabels_L-1;
BLabels_L = BLabels_L-1;

FLabels_R = FLabels_R-1;
BLabels_R = BLabels_R-1;


lambda_p = 10 ;       
lambda_c = 0.000001 ;  
[SegImage_L,SegImage_R] = StereoCut_GMM( r,c , Colors_L, FLabels_L, BLabels_L, FDist_L, BDist_L , ...
                                               Colors_R, FLabels_R, BLabels_R, FDist_R, BDist_R , ...
                                         disp , lambda_p , lambda_c );

                                     
%% post processing %%
[SegImage_L,SegImage_R] = post_processing(SegImage_L,SegImage_R,r,c) ;

%%  Display the segmented image   %%
SegImage_L = repmat(SegImage_L,[1,1,3]);
SegNewImage_L = uint8(SegImage_L) .* uint8(I_L);

SegImage_R = repmat(SegImage_R,[1,1,3]);
SegNewImage_R = uint8(SegImage_R) .* uint8(I_R);
%-----------------------------------------

pad = zeros(r,5,3)+255 ;
result = [SegNewImage_L,pad,SegNewImage_R] ;

name_save = strcat('result\\cut ',string_i,'.bmp') ;
imwrite(result, name_save ,'bmp') ;



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Helper Functions declarations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function Dist = ClustDistMembership(MeanColors, CClusters, Covs, Weights)
% CLUSTDISTMEMBERSHIP - Calcuates FG and BG Distances
% Authors - Mohit Gupta, Krishnan Ramnath
% Affiliation - Robotics Institute, CMU, Pittsburgh
% 2006-05-15

NumClusters = size(CClusters,2);
numSP = size(MeanColors,1);  

Dist = zeros(numSP,1);
Ind = zeros(numSP,1);

tmp = zeros(numSP, NumClusters);

for k = 1:NumClusters
    M = CClusters(:,k);  
    CovM = Covs(:,:,k); 
    W = Weights(1,k);    

    V = MeanColors - repmat(M',numSP,1);
    tmp(:,k) = 1000*(W / sqrt(det(CovM))) * exp(-( sum( ((V * inv(CovM)) .* V),2) /2));
    % = W  *    1/sqrt(|CovM|)      * exp[- sum(i) (Ii-Ii')^2 /2*inv(CovM)]  

end

Dist = sum(tmp,2);
        


