function [FCClusters, FCovs , FWeights ] = get_GMM( FColors, NumFClusters )
%%%%%%%%%%% Using Just kmeans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FId FCClusters] = kmeans( FColors, NumFClusters);
Fdim = size(FColors,2);

FCClusters = zeros(Fdim, NumFClusters);
FWeights = zeros(1,NumFClusters);
FCovs = zeros(Fdim, Fdim, NumFClusters);
for k=1:NumFClusters
    relColors = FColors(find(FId==k),:);       
    FCClusters(:,k) = mean(relColors,1)';
    FCovs(:,:,k) = cov(relColors);    
    FWeights(1,k) = length(find(FId==k)) / length(FId); 
end


