



function [All_Orig, All_Null, AllDelta, ConnOrig, DistOrig, DeltaOrig]=CorticalTypeScrambler(Species, indexes, percentages_range, what_to_do, iterations, allow_stretching)

AllDelta=[];ConnOrig=[]; DistOrig=[]; DeltaOrig=[];
flag=0;

for i=indexes


    fprintf('\nIterating %s...\n',Species(i).name);    


    C=Species(i).C;

    Dist=Species(i).Dist;
    ND=Species(i).CorticalType;
    %AP=Species(i).AP_Coords;
    %LM=Species(i).LM_Coords;
    %DV=Species(i).DV_Coords;
    %Names=Species(i).Names;
    DeltaOrig=Species(i).Delta;
    
    %C(isnan(C))=0;
    Dist(isnan(Dist))=0;

    %Calculate the geometric centrality. To adhere to the intuitive meaning of 
    %centrality, divide by max and reverse the order so that the higher value
    %correspond to areas that are central.
    %GC=sum(Dist,2)./(size(Dist,2)-1);
    %GC=GC./max(GC);
    %GC=1-GC;
    
    min_ND=min(ND);
    max_ND=max(ND);

    Orig_ND=ND;%Save the original measures since we have to "renew" it in
               %every cycle.
       
    %Perform the basic presence-absence analysis
    if(any(strcmp('existence',what_to_do)))
        
        
        for iter=1:iterations
        
       
            for perc=1:length(percentages_range)

                current_perc=percentages_range(perc);

                ind=1:length(ND);
                ind=ind(randperm(length(ind)));

                ind=ind(1:round(current_perc*length(find(~isnan(ND))))); 

                for l=1:length(ind)

                        if(~isnan(ND(ind(l))))

                            %If we hit the minimum, only option is to increase the
                            %index
                            if(ND(ind(l))==min_ND | ND(ind(l))==max_ND)

                                if(ND(ind(l))==min_ND)
                                   ND(ind(l))=ND(ind(l))+1; 
                                end

                                if((ND(ind(l))==max_ND) & (allow_stretching==1))
                                   ND(ind(l))=ND(ind(l))+1; 
                                end

                                if((ND(ind(l))==max_ND) & (allow_stretching==0))
                                   ND(ind(l))=ND(ind(l))-1; 
                                end 

                            else

                                %Decide if adding or subtracting should be
                                %performed by flipping a coin.
                                if(rand > 0.5)
                                     ND(ind(l))=ND(ind(l))+1;
                                else
                                     ND(ind(l))=ND(ind(l))-1;
                                end

                            end


                        end

                end

    
            for k=1:length(ND)
               
                for m=1:length(ND)
                    
                    Delta(k,m)=ND(k)-ND(m);
                    
                end
                
            end

            use=find(~isnan(Dist) & ~isnan(Delta) & (Dist > 0) & ~isnan(C));
            conn=C(use);
            conn=double(conn > 0);
            dist=Dist(use);
            delta=Delta(use);
            %ap=ap(use);
            
            if(flag==0)
                ConnOrig=conn;
                DistOrig=dist;
                DeltaOrig=DeltaOrig(use);
                flag=1;
            end
            
            AllDelta=horzcat(AllDelta, delta);

            [b, ~, ~] = glmfit(Rescale01(horzcat(dist,abs(delta))),conn,'binomial','link','logit');

            All_Orig(iter, perc)=b(3);
            
            %Compute stats of shuffled conn
            
            conn=conn(randperm(length(conn)));

            [b, ~, ~] = glmfit(Rescale01(horzcat(dist,abs(delta))),conn,'binomial','link','logit');

            All_Null(iter, perc)=b(3);

            ND=Orig_ND;

            end
    
        
        end
        
        
    end
    
    
    
    
    
    if(any(strcmp('core-periphery',what_to_do)))

        %calculate the core
        fprintf('\n Calculating core and null size...\n');

        [~, InCore, ~, ~, ~]=CorePeriphery(C);
        
    
    
        for iter=1:iterations


            for perc=1:length(percentages_range)

                current_perc=percentages_range(perc);


                ind=1:length(ND);
                ind=ind(randperm(length(ind)));

                ind=ind(1:round(current_perc*length(find(~isnan(ND)))));


                for l=1:length(ind)

                    if(~isnan(ND(ind(l))))

                        %If we hit the minimum, only option is to increase the
                        %index
                        if(ND(ind(l))==min_ND | ND(ind(l))==max_ND)

                            if(ND(ind(l))==min_ND)
                               ND(ind(l))=ND(ind(l))+1; 
                            end

                            if((ND(ind(l))==max_ND) & (allow_stretching==1))
                               ND(ind(l))=ND(ind(l))+1; 
                            end

                            if((ND(ind(l))==max_ND) & (allow_stretching==0))
                               ND(ind(l))=ND(ind(l))-1; 
                            end 

                        else

                            %Decide if adding or subtracting should be
                            %performed by flipping a coin.
                            if(rand > 0.5)
                                 ND(ind(l))=ND(ind(l))+1;
                            else
                                 ND(ind(l))=ND(ind(l))-1;
                            end

                        end


                    end

                end

                %Calculate the difference of the core and periphery.

                c1=ND(InCore==1);
                c1=c1(~isnan(c1));

                c2=ND(InCore==0);
                c2=c2(~isnan(c2));

                [~, ~, KSstatistic] = kstest2(c1, c2);
                All_Orig(iter, perc)=KSstatistic;

                %[p, phi_nm, phi_nm_boot] = minentest(c1, c2, 10);
                %All_Orig(iter, perc)=phi_nm;

                %Calculate also the SK stat based on permuted values.
                all=vertcat(c1,c2);
                all=all(randperm(length(all)));

                c1=all(1:length(c1));
                c2=all((length(c1)+1):end);

                [~, ~, KSstatistic] = kstest2(c1, c2);
                 All_Null(iter, perc)=KSstatistic;

                %[p, phi_nm, phi_nm_boot] = minentest(c1, c2, 10);
                %All_Null(iter, perc)=phi_nm;

                %Restore original ND
                ND=Orig_ND;

            end


        end
    
    
    end
    
       
end

return