

function [AllMetrics, AllMetrics_AreaWise]=GetMetrics(C, Dist, weighted, directed, area_wise)

AllMetrics_AreaWise=[];


if(weighted==1)
    

    if(directed==1)
        
        C_Norm=weight_conversion(C, 'normalize');%Normalize to 0-1 for the weighted clustering coefficient
        CC=clustering_coef_wd(C_Norm);
        BC=betweenness_wei(-log(C));
        
    end
    
    if(directed==0)
        
        C_Norm=weight_conversion(C, 'normalize');%Normalize to 0-1 for the weighted clustering coefficient
        CC=clustering_coef_wu(C_Norm);
        BC=betweenness_wei(-log(C));
        
    end
    
   
    
else
    
    C_Bin=double(C > 0);
    
    if(directed==1)
        CC=clustering_coef_bd(C_Bin);
        BC=betweenness_bin(C_Bin);
    end
    
    if(directed==0)
        CC=clustering_coef_bu(C_Bin);
        BC=betweenness_bin(C_Bin);
    end
    
    BC=BC';
    
end

CC_Mean=mean(CC);

%Get the average shortest path length
if(weighted==1)

%     if(directed==0)
        
        [SPL] = distance_wei_floyd(C,'log');
    
%     end
    
%     if(directed==1)
%         
%         [SPL] = distance_wei_floyd(C,'inv');
%         
%     end

else
    
    [SPL] = distance_bin(C_Bin);
    
end

mean_shortest_path=sum(sum(SPL))/((size(C,1).^2)-size(C,1));

%Get the average in out connectivity similarity

% if(weighted==1)
% 
% [Min,Mout,~]=connectivity_similarity_cosine(C); 
%     
% else
%     
% [Min,Mout,~] = matching_ind(C);
% Min=Min+Min';
% Mout=Mout+Mout';
% 
% end
% 
% mean_Min=sum(sum(Min))/((size(C,1).^2)-size(C,1));
% mean_Mout=sum(sum(Mout))/((size(C,1).^2)-size(C,1));

%Get the average length of connections
ind=find(C~=0);

%Normalize Dist by max distance
Dist=Dist./max(max(Dist));
fiber_length=nanmean(Dist(ind));


%Get modularity
if(weighted==1)
        
    if(directed==1)
        [~,Q]=modularity_dir(C,1);
    end
    
    if(directed==0)
        [~,Q]=modularity_und(C,1);
    end
    
else
    
    if(directed==1)
        [~,Q]=modularity_dir(C_Bin,1);
    end
    
    if(directed==0)
        [~,Q]=modularity_und(C_Bin,1);
    end
    
end

%Get average diffusion efficiency
if(weighted==1)
    
    [G_DE,DiffEff] = diffusion_efficiency(C);

else
    
    [G_DE,DiffEff] = diffusion_efficiency(C_Bin);
    
end

%Get the maximum size of a clique/core - run it on a binarized matrix

    C_Bin=double(C > 0);
    
    [~, ~, ~, CliqueSize, ~]=CorePeriphery(C_Bin);
    max_clique_size=max(CliqueSize);


%Summarize metrics in a vector
AllMetrics=horzcat(CC_Mean, mean_shortest_path, fiber_length, max_clique_size, Q, G_DE);

if(area_wise==1)
    
    
    if(directed==1)
    
      AllMetrics_AreaWise=horzcat(CC, BC, sum(SPL,2)./(size(SPL,2)-1), (sum(SPL,1)./(size(SPL,2)-1))', sum(DiffEff,2)./(size(DiffEff,2)-1), (sum(DiffEff,1)./(size(DiffEff,2)-1))');
    
    end
    
    if(directed==0)
        
      AllMetrics_AreaWise=horzcat(CC, BC, sum(SPL,2)./(size(SPL,2)-1), sum(DiffEff,2)./(size(DiffEff,2)-1));
        
    end
    
    
end





return

