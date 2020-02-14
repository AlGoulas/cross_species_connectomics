
%what_to_do:    string indicating the following options:
%               'existence':
%               'core-periphery':
%               'connectome-cyto-geo-ap':
%
%compare_distributions: string indicating the following options:
%                       'ks': 2D Kolmogorov-Smirnov test
%                       'statistical-energy': 2D statistical energy test

function [Stats, Plotting]=CrossSpecies_NetworkMetrics(Species, indexes, weighted, permutations, what_to_do, compare_distributions)

%Meter index is used to store info durign the iteration of each species
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

    %Collect the data for the current species for the struct Species
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
        
        %Convert to binary if weights=0
        if(weighted==0)
            C=double(C > 0);
        end

        use=find(~isnan(Dist) & ~isnan(Delta) & (Dist > 0));
        conn=C(use);
        conn=double(conn>0);
        dist=Dist(use);
        delta=Delta(use);
        
        %Run a logistic regression to predict present and absent
        %conenctions.
        %The model has distance between areas and cortical type differences
        %as predictors.
        [~, ~, stats] = glmfit(Rescale01(horzcat(dist,abs(delta))),conn,'binomial','link','logit');
        
        %Save the returned stat struct for further analsysi or inspection
        Existence_Stats(meter).stat_logistic_regr=stats;
        
        if(strcmp('ks', compare_distributions))
            
            [~, pValue_delta, KSstatistic_delta] = kstest2(abs(delta(conn==1)), abs(delta(conn==0)));
            [~, pValue_dist, KSstatistic_dist] = kstest2(abs(dist(conn==1)), abs(dist(conn==0)));
            
            %Keep results in a struct
            Existence_Stats(meter).pval_delta=pValue_delta;
            Existence_Stats(meter).pval_dist=pValue_dist;
            Existence_Stats(meter).stat_delta=KSstatistic_delta;
            Existence_Stats(meter).stat_dist=KSstatistic_dist;
            
        end    
            
        if(strcmp('statistical-energy', compare_distributions))
        
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

        %calculate the core. This is the union of all cliques
        %of maximum size in the connectome
        fprintf('\n Calculating core and null size...\n');

        [~, InCore, ~, CliqueSize, ~]=CorePeriphery(C);

        %Calculate the pvalue through permutations controlling for
        %degree-distribution (this is needed because with an heterogenous 
        %degree-distribution network WILL exhibit a core.

        MaxCliqueSizeNull=zeros(1,permutations);
        
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
        
        if(weighted==1)
        
            if(~isempty(find((C~=0) & (C~=1))))
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 1, directed, 1);
            else
                [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, 0, directed, 1);
            end
        
        else
            
            [~, AllMetrics_AreaWise]=GetMetrics(C, Dist, weighted, directed, 1);
            
        end
        
        AllMetrics_AreaWise=AllMetrics_AreaWise(:,3:end);%Use only centrality
        
        %Normalize for homogenous visualization across species
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


    end


    %Compute network metrics and relate to cyto, geo and AP.
    if(any(strcmp('connectome-cyto-geo-ap', what_to_do)))
        
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
        

        %AllMetrics_AreaWise=AllMetrics_AreaWise(:,1:end);
            
        AllMetrics_AreaWise_Norm=Rescale01(AllMetrics_AreaWise);

        %Summarize the features in 2D by performing a PCA.
        %In this way we have a network topology 2D coordinate system as 
        %a summary for the connectomes
        [coeff, score, ~, ~, explained, ~] = pca(AllMetrics_AreaWise_Norm);

        %Store the results of the pca for each species
        Connectome_CrossLevels_Stats(meter).pca_coeff=coeff;
        Connectome_CrossLevels_Stats(meter).pca_score=score;
        Connectome_CrossLevels_Stats(meter).pca_explained=explained;
        Connectome_CrossLevels_Stats(meter).net_metrics=AllMetrics_AreaWise;
        Connectome_CrossLevels_Stats(meter).net_metrics_names=labels_net_metrics;
        
        %We only take into account PC1 for plotting and correlations
        %The returned Stats structure contains all info for further
        %analysis of other PCs.
        pc=1;      
        
        %Correlate 1st component with AP, ND, GC.
        use=find((~isnan(AP)) & (~isnan(ND)) & (~isnan(GC)));  
        
        Connectome_CrossLevels_Plotting(meter).pc = score(use,pc);
        Connectome_CrossLevels_Plotting(meter).nd = ND(use);
        Connectome_CrossLevels_Plotting(meter).ap = AP(use);
        Connectome_CrossLevels_Plotting(meter).all_names = Names;
        Connectome_CrossLevels_Plotting(meter).reduced_names = Names(use);
        
        %Save r and pvals of partial and univariate correlations
        [r, p]=partialcorr(ND(use),score(use,pc),AP(use),'type','Spearman');
        
        Connectome_CrossLevels_Stats(meter).rho_pc_cyto_partial=r;
        Connectome_CrossLevels_Stats(meter).p_pc_cyto_partial=p;
        
        [r, p]=corr(ND(use),score(use,pc),'type','Spearman');
        
        Connectome_CrossLevels_Stats(meter).rho_pc_cyto=r;
        Connectome_CrossLevels_Stats(meter).pval_pc_cyto=p;
        
        
        [r, p]=partialcorr(AP(use),score(use,pc),ND(use),'type','Spearman');
        
        Connectome_CrossLevels_Stats(meter).rho_pc_ap_partial=r;
        Connectome_CrossLevels_Stats(meter).p_pc_ap_partial=p;
        
        [r, p]=corr(AP(use),score(use,pc),'type','Spearman');
        
        Connectome_CrossLevels_Stats(meter).rho_pc_ap=r;
        Connectome_CrossLevels_Stats(meter).pval_pc_ap=p;
        
    end

    meter=meter+1;
    
end

%Here the total output across the indicated species will be rendered 
%in figures
if(any(strcmp('existence', what_to_do)))


    for s=1:length(indexes)
    
   
        figure('Position',[0 0 500 500]);

        conn= Existence_Plotting(s).conn;
        delta=Existence_Plotting(s).delta;
        dist=Existence_Plotting(s).dist;


        scatterhist(Rescale01(dist),Rescale01(delta),'Group',conn,'Direction','out','Kernel','on');
    
    
    end
    

end



%Block of analysis for the core-periphery centric approach
if(any(strcmp('core-periphery',what_to_do)))

        figure('Position',[0 0 500 500]);
       
        %index for subplots
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
        
        for s=1:length(indexes)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            Cyto=Rescale01(CorePeriphery_Plotting(s).Cyto);

            %Plot cyto
            subplot(1,length(indexes),s);
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
           

        end
        
        suptitle('Cyto');
        
        
        %Plot the Anterior-Posterior coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
        
        for s=1:length(indexes)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            AP=Rescale01(CorePeriphery_Plotting(s).AP);

            %Plot AP
            subplot(1,length(indexes),s);
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
           
            
 

        end
        
        suptitle('AP');
        
        
        
        figure('Position',[0 0 500 200]);
        
        %Plot the geometric centrality for the core and periphery 
        %areas
        for s=1:length(indexes)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            GC=CorePeriphery_Plotting(s).GC;

            %Plot
            subplot(1,length(indexes),s);
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
           
           

        end
        
        suptitle('GC');
        
        %Plot the Dorsal-Ventral coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(indexes)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            DV=Rescale01(CorePeriphery_Plotting(s).DV);

            %Plot cyto
            subplot(1,length(indexes),s);
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
            
            


        end
        
        suptitle('DV');
        
        
        
        %Plot the Lateral-Medial coordinates for the core and periphery 
        %areas
        figure('Position',[0 0 500 200]);
       
        
        for s=1:length(indexes)
            
            InCore=CorePeriphery_Plotting(s).InCore;
            LM=Rescale01(CorePeriphery_Plotting(s).LM);

            %Plot cyto
            subplot(1,length(indexes),s);
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
            
            

        end
        
        suptitle('LM');
       
        %Plot the wiring cost for the core and periphery areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(indexes)
            
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
            subplot(1,length(indexes),s);
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

        end
        
        suptitle('WiringCost Core-Periphery');
        
        
        %Plot the projection lengths for the core and periphery areas
        figure('Position',[0 0 500 200]);
       
        for s=1:length(indexes)
            
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
            subplot(1,length(indexes),s);
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
            
            
            
           

        end
        
        suptitle('Projection Length Core-Periphery');  
        
end



if(any(strcmp('connectome-cyto-geo-ap',what_to_do)))
    
    for s=1:length(indexes)

        %Get the to-be-plotted data
        
        %The only info for the plotting that is in the _Stats struct
        %is for the biplot (coeff and scores) and the labels of the 
        %net metrics
        coeff = Connectome_CrossLevels_Stats(s).pca_coeff;
        score = Connectome_CrossLevels_Stats(s).pca_score;
        labels_net_metrics = Connectome_CrossLevels_Stats(s).net_metrics_names;
        
        names = Connectome_CrossLevels_Plotting(s).all_names;
        names_reduced = Connectome_CrossLevels_Plotting(s).reduced_names
        pc = Connectome_CrossLevels_Plotting(s).pc;
        nd = Connectome_CrossLevels_Plotting(s).nd;
        ap = Connectome_CrossLevels_Plotting(s).ap;
        
        figure;
        biplot(coeff(:,1:2), 'scores', score(:,1:2), 'varlabels', labels_net_metrics);
         
        figure;
        scatter(score(:,1),score(:,2));
        text(score(:,1),score(:,2),names);
        title('PCA space of connectome');

        figure;
        scatter(pc,nd);
        text(pc,nd,names_reduced);
        lsline;  
        
        r = Connectome_CrossLevels_Stats(s).rho_pc_cyto;
        p = Connectome_CrossLevels_Stats(s).pval_pc_cyto;
        
        title(strcat(Species(i).name,' ND rho=',num2str(r),'pval=',num2str(p))); 
        
        figure;
        scatter(pc,ap);
        text(pc,ap,names_reduced);
        lsline;
        
        r = Connectome_CrossLevels_Stats(s).rho_pc_ap;
        p = Connectome_CrossLevels_Stats(s).pval_pc_ap;
        
        title(strcat(Species(i).name,' AP rho=',num2str(r),'pval=',num2str(p)));
   
    end
        
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


if(any(strcmp('connectome-cyto-geo-ap',what_to_do)))
   
    Stats = Connectome_CrossLevels_Stats;
    Plotting = Connectome_CrossLevels_Plotting;
    
end


return



