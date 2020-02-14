
%what_to_do:    string indicating the following options:
%               'existence':
%               'core-periphery':
%               'connectome-cyto':
%

function [Stats, Plotting]=CrossSpecies_NetworkMetrics(Species, indexes, weighted, permutations, what_to_do, compare_distributions)

meter=1;

for i=indexes


    fprintf('\nIterating %s...\n',Species(i).name);    


    C=Species(i).C;

    %Decide if it is directed
    if(issymmetric(C))
       directed=0;
    else
       directed=1; 
    end

    Dist=Species(i).Dist;
    ND=Species(i).CorticalType;
    AP=Species(i).AP_Coords;
    LM=Species(i).LM_Coords;
    DV=Species(i).DV_Coords;
    Names=Species(i).Names;
    Delta=Species(i).Delta;

    C(isnan(C))=0;
    Dist(isnan(Dist))=0;

    %Calculate the geometric centrality. To adhere to the intuitive meaning of 
    %centrality, divide by max and reverse the order so that the higher value
    %correspond to areas that are central.
    GC=sum(Dist,2)./(size(Dist,2)-1);
    GC=GC./max(GC);
    GC=1-GC;


    %Perform the basic presence-absence analysis
    if(any(strcmp('existence',what_to_do)))
        
        if(weighted==0)
            C=double(C > 0);
        end

        %Check the AP contribution as well
        %ap=pdist(AP);
        %ap=squareform(ap);

        %use=find(~isnan(Dist) & ~isnan(Delta) & (Dist > 0) & ~isnan(ap));
        use=find(~isnan(Dist) & ~isnan(Delta) & (Dist > 0));
        conn=C(use);
        conn=double(conn>0);
        dist=Dist(use);
        delta=Delta(use);
        %ap=ap(use);
        
        [~, ~, stats] = glmfit(Rescale01(horzcat(dist,abs(delta))),conn,'binomial','link','logit');
        
        if(strcmp('ks', compare_distributions))
            
            [~, pValue_delta, KSstatistic_delta] = kstest2(abs(delta(conn==1)), abs(delta(conn==0)));
            [~, pValue_dist, KSstatistic_dist] = kstest2(abs(dist(conn==1)), abs(dist(conn==0)));
            
            %Keep results in a struct
            Existence_Stats(meter).pval_delta=pValue_delta;
            Existence_Stats(meter).pval_dist=pValue_dist;
            Existence_Stats(meter).stat_delta=KSstatistic_delta;
            Existence_Stats(meter).stat_dist=KSstatistic_dist;
            
        end    
            
        if(strcmp('statistical_energy', compare_distributions))
        
            [pValue_delta, MEstat_delta] = minentest(abs(delta(conn==1)), abs(delta(conn==0)), permutations);
            [pValue_dist, MEstat_dist] = minentest(abs(dist(conn==1)), abs(dist(conn==0)), permutations);
        
            %Keep results in a struct
            Existence_Stats(meter).pval_delta=pValue_delta;
            Existence_Stats(meter).pval_dist=pValue_dist;
            Existence_Stats(meter).stat_delta=MEstat_delta;
            Existence_Stats(meter).stat_dist=MEstat_dist;
            
        end
        
        %Save the needed info in a struct for ploting at the 
        %ebd of the analysis
        Existence_Plotting(meter).conn=conn;
        Existence_Plotting(meter).delta=abs(delta);
        Existence_Plotting(meter).dist=dist;
        Existence_Plotting(meter).Names=Names;
        Existence_Plotting(meter).stats=stats;
        
    end
    
    

    %Perform the core-periphery analysis
    if(any(strcmp('core-periphery',what_to_do)))

        %calculate the core
        fprintf('\n Calculating core and null size...\n');

        [~, InCore, ~, CliqueSize, ~]=CorePeriphery(C);

        %Calculate the pvalue thriugh permutations controlling for
        %degree-distribution (this is needed because with an heterogenous 
        %degree-distribution network WILL exhibit a core. The important issue is 
        %to control for the fact that we observe a core composed of cliques of a 
        %particular size).

        for perm=1:permutations

            [~, ~, ~, CliqueSizeNull, ~]=CorePeriphery(randmio_dir(C,10));
            MaxCliqueSizeNull(perm)=max(CliqueSizeNull);

        end
        
        %Size of maximum clique
        CorePeriphery_Stats(meter).max_clique_size = max(CliqueSize);
        
        %Indexes denioting the core areas (=1) 
        CorePeriphery_Stats(meter).core_indexes = InCore;
        
        %Empirical pvalue based on degree, edge and node matched null
        %network
        CorePeriphery_Stats(meter).max_clique_size_pval = length(find(MaxCliqueSizeNull >= max(CliqueSize))) / (permutations+1);
        
        %[AllMetrics, AllMetrics_AreaWise]=GetMetrics(C, Dist, weighted, 1, 1);
        
        if(weighted==1)
        
            if(~isempty(find((C~=0) & (C~=1))))
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 1, directed, 1);
            else
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 0, directed, 1);
            end
        
        else
            
            [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, weighted, directed, 1);
            
        end
        
        AllMetrics_AreaWise=AllMetrics_AreaWise(:,3:end);
        AllMetrics_AreaWise(:,1)=1./AllMetrics_AreaWise(:,1);
        AllMetrics_AreaWise(:,2)=1./AllMetrics_AreaWise(:,2);
        
        
        %Store the results for plotting
        CorePeriphery_Plotting(meter).InCore=InCore;
        CorePeriphery_Plotting(meter).AllMetrics_AreaWise=AllMetrics_AreaWise;
        CorePeriphery_Plotting(meter).Cyto=ND;
        CorePeriphery_Plotting(meter).AP=AP;
        CorePeriphery_Plotting(meter).DV=DV;
        CorePeriphery_Plotting(meter).LM=LM;
        CorePeriphery_Plotting(meter).GC=GC;
        CorePeriphery_Plotting(meter).Dist=Dist;
        CorePeriphery_Plotting(meter).C=C;
        
        %Plot here the core-periphery difference in terms of in-out
        %efficiency and diffusion efficiency.
        
%         subplot(4,4,meter);
%         
%         %Plot out and in efficiency
%         boxplot(vertcat(AllMetrics_AreaWise(InCore==1, 1), AllMetrics_AreaWise(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
%         boxplot(vertcat(AllMetrics_AreaWise(InCore==1, 2), AllMetrics_AreaWise(InCore==0, 2)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
% 
%         %Plot out and in diffusion efficiency
%         boxplot(vertcat(AllMetrics_AreaWise(InCore==1, 3), AllMetrics_AreaWise(InCore==0, 3)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
%         boxplot(vertcat(AllMetrics_AreaWise(InCore==1, 4), AllMetrics_AreaWise(InCore==0, 4)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));

        %Calculate shortest path length and diffusion efficiency
        % fprintf('\n Calculating shortest path and diffusion efficiency...\n');
        % 
        % if(weighted==1)
        %     
        %     [~,Ediff] = diffusion_efficiency(C);
        %     
        %     DE_Out=sum(Ediff,2)./((size(Ediff,2)-1));
        %     DE_In=sum(Ediff,1)./((size(Ediff,2)-1));
        %     
        %     [SPL,hops,Pmat] = distance_wei_floyd(C,'log');
        %     
        %     SP_Out=sum(SPL,2)./((size(SPL,2)-1));
        %     SP_In=sum(SPL,1)./((size(SPL,2)-1));
        %     
        % else
        %     
        %     C_Bin=C;
        %     C_Bin=double(C_Bin > 0);
        %     
        %     [~,Ediff] = diffusion_efficiency(C_Bin);
        %     
        %     DE_Out=sum(Ediff,2)./((size(Ediff,2)-1));
        %     DE_In=sum(Ediff,1)./((size(Ediff,2)-1));
        %     DE_In=DE_In';
        %     
        %     
        %     [SPL] = distance_bin(C_Bin);
        %     
        %     SP_Out=sum(SPL,2)./((size(SPL,2)-1));
        %     SP_In=sum(SPL,1)./((size(SPL,2)-1));
        %     SP_In=SP_In';
        %     
        % end

        %Plot the summary of the metrics for the core-periphery
        %figure;
        %boxplot(vertcat(DE_Out(InCore==1),DE_Out(InCore==0)),vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
        %boxplot(vertcat(DE_In(InCore==1),DE_In(InCore==0)),vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
        %boxplot(vertcat(GC(InCore==1),GC(InCore==0)),vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
        %[H, pValue, KSstatistic] = kstest2(DE_Out(InCore==1),DE_Out(InCore==0));
        %pause;


    end


    %Compute network metrics and relate to cyto, geo and AP.
    if(any(strcmp('connectome-cyto', what_to_do)))
        
        if(weighted==1)

            if(~isempty(find((C~=0) & (C~=1))))
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 1, directed, 1);
            else
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 0, directed, 1);
            end
        
        else
            
            [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 0, directed, 1);
            
        end

        if(directed==1)
            labels_net_metrics={'CC','BC','SPL_out','SPL_in','DE_out','DE_in'};
        end
        
        if(directed==0)
            labels_net_metrics={'CC','BC','SPL','DE'};
        end
        
        %labels_net_metrics={'SPL_out','SPL_in','DE_out','DE_in'};
        AllMetrics_AreaWise=AllMetrics_AreaWise(:,1:end);
        
        %Invert the shortest paths to convert to efficiency
        %AllMetrics_AreaWise(:,1)=1./AllMetrics_AreaWise(:,1);
        %AllMetrics_AreaWise(:,2)=1./AllMetrics_AreaWise(:,2);

        %AllMetrics_AreaWise(:,3)=1./AllMetrics_AreaWise(:,3);
        %AllMetrics_AreaWise(:,4)=1./AllMetrics_AreaWise(:,4);
        
        AllMetrics_AreaWise_Norm=Rescale01(AllMetrics_AreaWise);

        %[AllMetrics, AllMetrics_AreaWise]=GetMetrics(C, Dist, weighted, 1, 1);

%         AllMetrics_AreaWise_Norm=Rescale01(AllMetrics_AreaWise);
%         use=find((~isnan(AP)) & (~isnan(ND)) & (~isnan(GC)));
%         
%         for m=1:size(AllMetrics_AreaWise_Norm,2)
%             
%              [r, p]=partialcorr(ND(use),AllMetrics_AreaWise_Norm(use,m),horzcat(AP(use)),'type','Spearman');
%              [r, p]=partialcorr(AP(use),AllMetrics_AreaWise_Norm(use,m),horzcat(ND(use)),'type','Spearman');
%              
%             %[r, p]=corr(AP(use),score(use),'type','Spearman');
%             %title(strcat(Species(i).name,' AP rho=',num2str(r),'pval=',num2str(p)));
%             
%         end
        

        [coeff, score, ~, ~, ~, ~] = pca(AllMetrics_AreaWise_Norm);

        pc=1;
        
        figure;
        biplot(coeff(:,1:2), 'scores', score(:,1:2), 'varlabels', labels_net_metrics);
         
        figure;
        scatter(score(:,1),score(:,2));
        text(score(:,1),score(:,2),Names);
        title('PCA space of connectome');
        
        %Correlate 1st component with AP, ND, GC.
        use=find((~isnan(AP)) & (~isnan(ND)) & (~isnan(GC)));

        figure;
        scatter(score(use,pc),ND(use));
        text(score(use,pc),ND(use),Names(use));
        lsline;
        
        %[r, p]=partialcorr(ND(use),score(use,pc),horzcat(AP(use)),'type','Spearman');
        [r, p]=corr(ND(use),score(use,pc),'type','Spearman');
        title(strcat(Species(i).name,' ND rho=',num2str(r),'pval=',num2str(p)));

        figure;
        scatter(score(use,pc),AP(use));
        text(score(use),AP(use),Names(use));
        lsline;
        
        %[r, p]=partialcorr(AP(use),score(use,pc),horzcat(ND(use)),'type','Spearman');
        [r, p]=corr(AP(use),score(use,pc),'type','Spearman');
        title(strcat(Species(i).name,' AP rho=',num2str(r),'pval=',num2str(p)));
        
        figure;
        scatter(score(use,pc),LM(use));
        text(score(use),LM(use),Names(use));
        lsline;
        
        %[r, p]=partialcorr(AP(use),score(use,pc),horzcat(ND(use)),'type','Spearman');
        [r, p]=corr(LM(use),score(use,pc),'type','Spearman');
        title(strcat(Species(i).name,' LM rho=',num2str(r),'pval=',num2str(p)));
        
        figure;
        scatter(score(use,pc),DV(use));
        text(score(use),DV(use),Names(use));
        lsline;
        
        %[r, p]=partialcorr(AP(use),score(use,pc),horzcat(ND(use)),'type','Spearman');
        [r, p]=corr(DV(use),score(use,pc),'type','Spearman');
        title(strcat(Species(i).name,' DV rho=',num2str(r),'pval=',num2str(p)));

%         figure;
%         scatter(score(use),GC(use));
%         text(score(use),GC(use),Names(use));
%         lsline;
%         
%         [r, p]=partialcorr(GC(use),score(use),horzcat(ND(use),AP(use)),'type','Spearman');
%         title(strcat(Species(i).name,' GC rho=',num2str(r),'pval=',num2str(p)));  


%         figure;
%         scatter(score(use),ND(use));
%         text(score(use),ND(use),Names(use));
%         [r, p]=corr(ND(use),score(use),'type','Spearman');
%         title(strcat(Species(i).name,' ND rho=',num2str(r),'pval=',num2str(p)));
% 
%         figure;
%         scatter(score(use),AP(use));
%         text(score(use),AP(use),Names(use));
%         [r, p]=corr(AP(use),score(use),'type','Spearman');
%         title(strcat(Species(i).name,' AP rho=',num2str(r),'pval=',num2str(p)));
% 
%         figure;
%         scatter(score(use),GC(use));
%         text(score(use),GC(use),Names(use));
%         [r, p]=corr(GC(use),score(use),'type','Spearman');
%         title(strcat(Species(i).name,' GC rho=',num2str(r),'pval=',num2str(p)));   


    end

    meter=meter+1;
    
end

%Here the total output across the indicated species will be rendered 
%in figures
if(any(strcmp('existence', what_to_do)))


    for s=1:length(Existence_Plotting)
    
   
        figure('Position',[0 0 500 500]);

        conn= Existence_Plotting(s).conn;
        delta=Existence_Plotting(s).delta;
        dist=Existence_Plotting(s).dist;
        %Existence_Plotting(s).Names=Names;

        scatterhist(Rescale01(dist),Rescale01(delta),'Group',conn,'Direction','out','Kernel','on');

%         scatter(Rescale01(dist(conn==1)),Rescale01(delta(conn==1)),'ro');
%         xlim([-0.1 1.1]);
%         ylim([-0.1 1.1]);
% 
%         figure('Position',[0 0 500 500]);
%         scatter(Rescale01(dist(conn==0)),Rescale01(delta(conn==0)),'bo');
%         xlim([-0.1 1.1]);
%         ylim([-0.1 1.1]);
%     
%     stats=Existence_Plotting(s).stats;
%     
%     subplot(2,2,s);
%    
%     errorbar(horzcat(stats.beta(2),stats.beta(3)),horzcat(stats.se(2),stats.se(3)));
    
    
    end
    

end



%Block of analysis for the core-periphery centric approach
if(any(strcmp('core-periphery',what_to_do)))

        figure('Position',[0 0 500 500]);
       
        %index for submplots
        ind=1;
        
        %Plot all the network centrality metrics for the core and periphery 
        %areas
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            AllMetrics_AreaWise=CorePeriphery_Plotting(s).AllMetrics_AreaWise;
            
            AllMetrics_AreaWise=Rescale01(AllMetrics_AreaWise);

            for kk=1:size(AllMetrics_AreaWise,2)
            
                %Plot out and in efficiency
                subplot(length(CorePeriphery_Plotting),size(AllMetrics_AreaWise,2),ind);
                boxplot(vertcat(AllMetrics_AreaWise(InCore==1, kk), AllMetrics_AreaWise(InCore==0, kk)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
                ind=ind+1;
            
            end

        end
        
        suptitle('Efficiency-Diffusion Efficiency');
        
        
        %Plot the cytology of the core and periphery areas
        figure('Position',[0 0 500 200]);
        
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            Cyto=Rescale01(CorePeriphery_Plotting(s).Cyto);

            %Plot cyto
            subplot(1,4,s);
            boxplot(vertcat(Cyto(InCore==1, 1), Cyto(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
            ylim([-0.1 1.1]);
            
            c1=Cyto(InCore==1);
            c1=c1(~isnan(c1));
            
            c2=Cyto(InCore==0);
            c2=c2(~isnan(c2));
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_cyto_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_cyto_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_cyto_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_cyto_core_periphery = pValue; 
   
            end
            
            
            title(num2str(pValue));
     

        end
        
        suptitle('Cyto');
        
        
        %Plot the Anterior-Posterior coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
        
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            AP=Rescale01(CorePeriphery_Plotting(s).AP);

            %Plot AP
            subplot(1,4,s);
            boxplot(vertcat(AP(InCore==1, 1), AP(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
            %ind=ind+1;
            
            c1=AP(InCore==1);
            c1=c1(~isnan(c1));
            
            c2=AP(InCore==0);
            c2=c2(~isnan(c2));
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_AP_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_AP_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_AP_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_AP_core_periphery = pValue; 
   
            end
           
            title(num2str(pValue));
 

        end
        
        suptitle('AP');
        
        
        
        figure('Position',[0 0 500 200]);
        
        %Plot the geometric centrality for the core and periphery 
        %areas
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            GC=CorePeriphery_Plotting(s).GC;

            %Plot
            subplot(1,4,s);
            boxplot(vertcat(GC(InCore==1, 1), GC(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
            
            c1=GC(InCore==1);
            c1=c1(~isnan(c1));
            
            c2=GC(InCore==0);
            c2=c2(~isnan(c2));
            
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_geo_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_geo_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_geo_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_geo_core_periphery = pValue; 
   
            end
           
            title(num2str(pValue));

        end
        
        suptitle('GC');
        
        %Plot the Dorsal-Ventral coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            DV=Rescale01(CorePeriphery_Plotting(s).DV);

            %Plot cyto
            subplot(1,4,s);
            boxplot(vertcat(DV(InCore==1, 1), DV(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
            
            c1=DV(InCore==1);
            c1=c1(~isnan(c1));
            
            c2=DV(InCore==0);
            c2=c2(~isnan(c2));
            
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_DV_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_DV_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_DV_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_DV_core_periphery = pValue; 
   
            end
            
            title(num2str(pValue));


        end
        
        suptitle('DV');
        
        
        
        %Plot the Lateral-Medial coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
       
        
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            LM=Rescale01(CorePeriphery_Plotting(s).LM);

            %Plot cyto
            subplot(1,4,s);
            boxplot(vertcat(LM(InCore==1, 1), LM(InCore==0, 1)), vertcat(ones(length(find(InCore==1)),1),2.*ones(length(find(InCore==0)),1)));
            
            c1=LM(InCore==1);
            c1=c1(~isnan(c1));
            
            c2=LM(InCore==0);
            c2=c2(~isnan(c2));
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_LM_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_LM_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_LM_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_LM_core_periphery = pValue; 
   
            end
            
            title(num2str(pValue));

        end
        
        suptitle('LM');
       
        %Plot the wiring cost for the core and periphery areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            Dist=CorePeriphery_Plotting(s).Dist;
            C=CorePeriphery_Plotting(s).C;
            
            if(weighted==0)
                C=double(C > 0);
            end
                
            WiringCost=Dist.*C;
            WiringCost_Out=Rescale01(nansum(WiringCost,2));
            WiringCost_In=Rescale01(nansum(WiringCost,1))';
            
            WiringCost_Core=vertcat(WiringCost_Out(InCore==1),WiringCost_In(InCore==1));
            WiringCost_Periphery=vertcat(WiringCost_Out(InCore==0),WiringCost_In(InCore==0));

            c1=WiringCost_Core;
            c2=WiringCost_Periphery;
            
            %Plot wiring cost of core and periphery
            subplot(1,4,s);
            boxplot(vertcat(WiringCost_Core,WiringCost_Periphery), vertcat(ones(length(find(InCore==1)),1), ones(length(find(InCore==1)),1), 2.*ones(length(find(InCore==0)),1), 2.*ones(length(find(InCore==0)),1)));

             if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_wiring_cost_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_wiring_cost_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_wiring_cost_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_wiring_cost_core_periphery = pValue; 
   
            end
            

            title(num2str(pValue));
            

        end
        
        suptitle('WiringCost Core-Periphery');
        
        
        %Plot the projection lengths for the core and periphery areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(CorePeriphery_Plotting)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            Dist=CorePeriphery_Plotting(s).Dist;
            C=CorePeriphery_Plotting(s).C;
                       
            if(weighted==0)
                C=double(C > 0);
            end
            
            Core_Lengths=[];
            Periphery_Lengths=[];
            
            for kk=1:size(C,1)
                
                ind=find((C(kk,:) > 0) & not(isnan(C(kk,:))));
                lengths_out=Dist(kk,ind)';
                
                ind=find((C(:,kk) > 0) & not(isnan(C(:,kk))));
                lengths_in=Dist(ind,kk);
                
                if(InCore(kk)==1)
                    Core_Lengths=vertcat(Core_Lengths,lengths_out,lengths_in);
                else
                    Periphery_Lengths=vertcat(Periphery_Lengths,lengths_out,lengths_in);
                end
                
            end

            
            %Plot projection length of core and periphery
            subplot(1,4,s);
            boxplot(vertcat(Core_Lengths,Periphery_Lengths), vertcat(ones(length(Core_Lengths),1),2.*ones(length(Periphery_Lengths),1)));
            
            c1=Core_Lengths;
            c2=Periphery_Lengths;
            
            if(strcmp(compare_distributions, 'ks'))
            
                [pValue, ks] = kstest2(c1, c2);
                
                CorePeriphery_Stats(s).stat_projection_length_core_periphery = ks; 
                CorePeriphery_Stats(s).pval_projection_length_core_periphery = pValue; 
            
            end
            
            if(strcmp(compare_distributions, 'statistical_energy'))
            
                [pValue, ME] = minentest(c1, c2, permutations);
                
                CorePeriphery_Stats(s).stat_projection_length_core_periphery = ME; 
                CorePeriphery_Stats(s).pval_rojection_length_core_periphery = pValue; 
   
            end
            
            
            title(num2str(pValue));
           

        end
        
        suptitle('Projection Length Core-Periphery');  
        
end

%Assign structures back to 2 basic output structures.
%This is not necessary to do, but for the sake of readability of the code,
%the struct names take the name of the what_to_do switch as prefix.

if(any(strcmp('existence',what_to_do)))
   
    Stats = Existence_Stats;
    Plotting = Existence_Plotting;
    
end

if(any(strcmp('core-periphery',what_to_do)))
   
    Stats = CorePeriphery_Stats;
    Plotting = CorePeriphery_Plotting;
    
end

return



