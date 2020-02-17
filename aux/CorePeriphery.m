



function [MC_Full, InCore, FreqInCore, CliqueSize, pval]=CorePeriphery(C)

%Take into account only the reciprocal conenctions - only such connections
%can be part of cliques by definition

C=double(C~=0);

S=C+C';
S(S==1)=0;
S(S==2)=1;

MC=ELSclique(S);
MC_Full=full(MC);

CliqueSize=sum(MC_Full,1);

%Get maximum clique size

max_size=max(CliqueSize);
ind=find(CliqueSize==max_size);

InCore=(MC_Full(:,ind));
InCore=sum(InCore,2);
FreqInCore=InCore;
InCore=double(InCore~=0);

%Compute analytically the probability of finding such core in such a network
%after Ercsey-Ravasz et al. 2013

%Number of nodes
%n=size(C,1);

%Network density
[kden,N,K]=density_dir(C);

%Core density
[kden_core,N_core,K_core]=density_dir(C(InCore==1,InCore==1));

%Possible rewirings in core
possible_rewirings=round((N_core^2-N_core)-K_core);

max_core=N_core^2-N_core;

kden_power_K_core=kden^K_core;

one_minus_kden_power_possible_rewirings=(1-kden)^possible_rewirings;

%Factorials of integers might lead to Infs so replace them by the highest
%integer resulting in the real maxvalue based on Matlab double
%specifications

the_sky_reached=0;

factorial_N=factorial(N);
factorial_N_core=factorial(N_core);
factorial_max_core=factorial(max_core);
factorial_possible_rewirings=factorial(possible_rewirings);
factorial_N_N_core=factorial(N-N_core);
factorial_max_core_possible_rewirings=factorial(max_core-possible_rewirings);

if(isinf(factorial(N)))
    factorial_N=realmax;
    the_sky_reached=1;
end

if(isinf(factorial(N_core)))
    factorial_N_core=realmax;
    the_sky_reached=1;
end

if(isinf(factorial(max_core)))
    factorial_max_core=realmax;
    the_sky_reached=1;
end

if(isinf(factorial(possible_rewirings)))
    factorial_possible_rewirings=realmax;
    the_sky_reached=1;
end

if(isinf(factorial(N-N_core)))
    factorial_N_N_core=realmax;
    the_sky_reached=1;
end

if(isinf(factorial(max_core-possible_rewirings)))
    factorial_max_core_possible_rewirings=realmax;
    the_sky_reached=1;
end
    
%p_value
if(the_sky_reached==1)
    
    pval=(factorial_N/(factorial_N_core*factorial_N_N_core))...
        *(factorial_max_core/(factorial_possible_rewirings*factorial_max_core_possible_rewirings))...
        *(kden_power_K_core)...
        *(one_minus_kden_power_possible_rewirings);

    
else
    
    pval=(factorial(N)/(factorial(N_core)*factorial(N-N_core)))...
        *(factorial(max_core)/(factorial(possible_rewirings)*factorial(max_core-possible_rewirings)))...
        *(kden_power_K_core)...
        *(one_minus_kden_power_possible_rewirings);

end

return;

