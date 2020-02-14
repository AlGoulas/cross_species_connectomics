


function  [all_metric_across_iterations, combos] = Bootstrapping(grouping_index, variables, iterations)

%Possible combinations given the number of the variables
kk=1;

for i=1:size(variables,2)
    
   c = combnk(1:size(variables,2),i);
   
   for n=1:size(c,1)
      
       combos{kk}=c(n,:);
       kk=kk+1;
       
   end
    
end


%For all combinations of the variables, calculate the metric based on
%the grouping index in a bootstraping way.

%Get the 2 group indexes. The entries in grouping_index must be 1 or 0. 
gr1=find(grouping_index==1);
gr2=find(grouping_index==0);

for c=1:length(combos)
    
    vars_to_use=combos{c};
    
    gr1_variables=variables(gr1, vars_to_use);
    gr2_variables=variables(gr2, vars_to_use);
    
    %Perform the resampling as many times as the variable iterations.
    
    metric_across_iterations=zeros(iterations,1);
    
    for iter=1:iterations
        
        %Construct the bootstraped samples for gr1 and gr2.
        boot1=floor(rand(1,length(gr1_variables)).*length(gr1_variables))+1;
        boot2=floor(rand(1,length(gr2_variables)).*length(gr2_variables))+1;
        
        boot_sample1=gr1_variables(boot1,:);
        boot_sample2=gr2_variables(boot2,:);
        
        %Statistical energy
        %[~,D] = minentest(boot_sample1, boot_sample2, 1);%No need to run permutations since we need the actual value of the metric across bootstrapping iterations
        
        
        
        %2D KS test
        
        if(size(boot_sample1,2) > 1)
            
            [H, pValue, D] = kstest_2s_2d(boot_sample1,boot_sample2);
        
        else
            
            [H, pValue, D] = kstest2(boot_sample1,boot_sample2);
             
        end
        
%         try
%             [~,D] = kstest2d(boot_sample1,boot_sample2);
%         catch
%             f=1;
%         end
            
%         try
%             [H, pValue, D] = kstest2(boot_sample1,boot_sample2);
%         catch
%            f=1; 
%         end
        
        
        metric_across_iterations(iter,1)=D;
        
    end
    
    
    all_metric_across_iterations(:,c)=metric_across_iterations;
    
end




return