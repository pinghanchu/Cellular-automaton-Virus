%20210112 Pinghan Chu
% The sum of virus_prod and ifn_prod is a constant.
% Tune the sum so that <ifn_prod> = 1,3,6
% Only 10% of cells can generate ifn_prod
%20210120 Pinghan Chu
% 10% IFN active with the same value
% Protected cells have a lifespan 100
clear all
close all

P = py.sys.path

addpath 'model6';
date = '20211206';

num_model_runs=100;
num_steps=365*5;
num_model_states=20;
eventnumb=10000;
grid_size=[100,100];
array_size=[grid_size,num_model_states];
%vars = importdata('input.xlsx');
%parms=vars.data;

distribution=["negbinomial"];
virus_prod=[2,1.8];
virus_diff=[0.5];
ifn_prob=[10];%ifn cell percentage
virus_reduction_factor=[10];%percent

%ifn_prod=[0,5];
%res_ratio=[1,0.1];
%virus_type_ratio=[50,10];

%ifn_prod=[1,2,3,4,6,7,8,9,10];
%res_ratio=[1,0.1];
%virus_type_ratio=[50,10];

%ifn_prod=[5];
%res_ratio=[0.5,0.05,0.01];
%virus_type_ratio=[50,10];%
ifn_prod=[5];
res_ratio=[0.1];
virus_type_ratio=[0,5,20,30,100];
% 
%res_ratio=[1];
%virus_type_ratio=[30,100];
%0,5,20


lifespan_mean = 10;%2 day=(5 steps)
lifespan_sigma = 2;
virus_prod_delay=5;
ifn_prod_delay=virus_prod_delay;
protected_lifespan=25;
dead_lifespan=25;
prob_infect=0.2;

parms(1)=virus_prod(1);
parms(2)=virus_prod(2);
parms(3)=virus_diff(1);
parms(4)=virus_prod_delay; % 5
parms(5)=ifn_prod(1);
parms(6)=virus_diff(1)*5;
parms(7)=ifn_prod_delay; % 5
parms(8)=protected_lifespan;
parms(9)=dead_lifespan;
parms(10)=prob_infect;
parms(11)=virus_reduction_factor(1);

scale=1;
a_res=[0,1,2];

xvirus_proda = virus_prod(1);
xvirus_prodb = virus_prod(2);
parms(1)=xvirus_proda;
parms(2)=xvirus_prodb;
disp([xvirus_proda,xvirus_prodb])
for ivirus_diff = 1:length(virus_diff)
    xvirus_diff = virus_diff(ivirus_diff);
    parms(3)=xvirus_diff;
    parms(6)=5*parms(3);%ifn_diff = 5*virus_diff   
    for iifn_prod = 1:length(ifn_prod)
        xifn_prod = ifn_prod(iifn_prod);
        parms(5)=xifn_prod;
        for iifn_prob = 1:length(ifn_prob)
            xifn_prob = ifn_prob(iifn_prob)*0.01;
            a_ifn_prob=[0,1];%ifn_prods is not the same; only a centain amount of ifn has effects.
            p_ifn_prob=[1-xifn_prob, xifn_prob];
            for idistribution = 1:length(distribution)
                xdistribution = distribution(idistribution);
                for ivirus_reduction_factor = 1:length(virus_reduction_factor)
                    parms(11)=virus_reduction_factor(ivirus_reduction_factor);
                    for ires_ratio = 1:length(res_ratio)
                        xres_ratio=res_ratio(ires_ratio)/100;
                        for ivirus_type_ratio = 1:length(virus_type_ratio)
                            xvirus_type_ratio=virus_type_ratio(ivirus_type_ratio)/100;
                            p_res=[1-xres_ratio,xres_ratio*(1-xvirus_type_ratio),xres_ratio*xvirus_type_ratio];
                            filename = [num2str(xdistribution),'_virusproda',num2str(xvirus_proda),'_virusprodb',num2str(xvirus_prodb),'_virusdiff',num2str(xvirus_diff),'_ifnprod',num2str(xifn_prod),'_ifnprob',num2str(ifn_prob(iifn_prob)),'_resratio',num2str(res_ratio(ires_ratio)),'_virustyperatio',num2str(virus_type_ratio(ivirus_type_ratio)),'_dim',num2str(grid_size(1)),'_sim',num2str(num_model_runs),'_step',num2str(num_steps)];
                            disp(filename);
                            disp(p_res);
                            %initialize output parameters
                            infected_count=zeros(num_model_runs,num_steps);
                            protected_count=zeros(num_model_runs,num_steps);
                            dead_count=zeros(num_model_runs,num_steps);
                            reservoir_count=zeros(num_model_runs,num_steps);
                            virus_count=zeros(num_model_runs,num_steps);
                            ifn_count=zeros(num_model_runs,num_steps);
                            virusproda_cell=zeros(num_model_runs,num_steps);
                            virusprodb_cell=zeros(num_model_runs,num_steps);
                            ifn_cell=zeros(num_model_runs,num_steps);
                            virus_contact=zeros(num_model_runs,num_steps);
                            ifn_contact=zeros(num_model_runs,num_steps);
                            virus_type1 = zeros(num_model_runs,num_steps);
                            virus_type2 = zeros(num_model_runs,num_steps);

                            for j=1:num_model_runs
                                if mod(j,20)==0
                                    disp([num2str(j),' simulations complete']);
                                end
                                %intital input matres
                                grid=initialize_grid(array_size);
                                %create ifn_prod matrix
                                %py.py_random.choice(a_ifn_prob,p_ifn_prob,eventnumb);         
                                %A_ifn_prod = dlmread('choice.csv');  
                                A_ifn_prod = rand(eventnumb,1)<xifn_prob;
                                A_ifn_prod = A_ifn_prod*(xifn_prod/xifn_prob);
                                M_ifn_prod=reshape(A_ifn_prod,array_size(1),array_size(2));
                                grid(:,:,11)=M_ifn_prod;

                                %create virus_prod matrix
                                switch xdistribution
                                    case "delta"
                                        scale = xvirus_proda;
                                        py.py_random.delta(scale,eventnumb);
                                    case "geometric"
                                        scale=1./(xvirus_proda+1.);
                                        py.py_random.geometric(scale,eventnumb);                        
                                    case "poisson"
                                        scale=xvirus_proda;
                                        py.py_random.poisson(scale,eventnumb);                    
                                    case "negbinomial"
                                        r=0.5;
                                        p=1.-xvirus_proda/(xvirus_proda+r);
                                        perl('negbinomial.pl',compose("%9.2f",r),compose("%9.7f",p),compose("%9.0f",eventnumb))
                                    
                                    otherwise
                                        disp('no this distribution')
                                end
                                A_virus_proda = dlmread('virus_prod.csv');
                                M_virus_proda = reshape(A_virus_proda,array_size(1),array_size(2));
                                grid(:,:,7)=M_virus_proda;
                                grid(:,:,9)=M_virus_proda;
                                %create virus_prod matrix
%                                 
                                switch xdistribution
                                    case "delta"
                                        scale = xvirus_prodb;
                                        py.py_random.delta(scale,eventnumb);
                                    case "geometric"
                                        scale=1./(xvirus_prodb+1.);
                                        py.py_random.geometric(scale,eventnumb);                        
                                    case "poisson"
                                        scale=xvirus_prodb;
                                        py.py_random.poisson(scale,eventnumb);                    
                                    case "negbinomial"
                                        r=0.5;
                                        p=1.-xvirus_prodb/(xvirus_prodb+r);
                                        %py.py_random.negbinomial(r,p,eventnumb);
                                        
                                        perl('negbinomial.pl',compose("%9.2f",r),compose("%9.7f",p),compose("%9.0f",eventnumb))
                                    otherwise
                                        disp('no this distribution')
                                end
                                A_virus_prodb = dlmread('virus_prod.csv');
                                M_virus_prodb = reshape(A_virus_prodb,array_size(1),array_size(2));
                                grid(:,:,8)=M_virus_prodb;
                                grid(:,:,10)=M_virus_prodb;
                                
                                %create infected cell lifespan
                                %py.py_random.normal(lifespan_mean,lifespan_sigma,eventnumb);                             
                                %A_infected_lifespan = dlmread('guess.csv');
                                p1rand = rand(eventnumb,1);
                                p2rand = rand(eventnumb,1);
                                    
                                A_infected_lifespan = sqrt(-2*log(p1rand)).*cos(2*pi*p2rand)*lifespan_sigma+lifespan_mean;

                                A_infected_lifespan = round(A_infected_lifespan);
                                for k=1:array_size(1)*array_size(2)
                                    if A_infected_lifespan(k)<1
                                        A_infected_lifespan(k) = 0;
                                    end
                                end
                                M_infected_lifespan = reshape(A_infected_lifespan,array_size(1),array_size(2));
                                grid(:,:,16) = M_infected_lifespan;

                                %create reservoir cell lifespan
                                M_reservoir_lifespan = ceil(rand(array_size(1),array_size(2))*365*5);
                                grid(:,:,17)=M_reservoir_lifespan;

                                %create reservoir cell
                                %py.py_random.choice(a_res,p_res,eventnumb);         
                                %A_reservoir=dlmread('choice.csv');
                                
                                prand = rand(eventnumb,1);
                                
                                A_reservoir = zeros(eventnumb,1);
                                prand = rand(eventnumb,1);
                                for i = 1: eventnumb
                                    if prand(i)<p_res(1)
                                        A_reservoir(i) = a_res(1);
                                    elseif (prand(i)>= p_res(1)) && (prand(i)< p_res(1)+p_res(2))
                                        A_reservoir(i) = a_res(2);
                                    else 
                                        A_reservoir(i) = a_res(3);
                                    end
                                end
                                M_reservoir = reshape(A_reservoir,array_size(1),array_size(2)); 
                                grid(:,:,4) = M_reservoir;

                                for step=1:num_steps
                                    %disp([num2str(step),' step']);
                                    [grid,C,C2,C3,C4,C5] = NoVirusFIP_step(grid, parms);

                                    infected_count(j,step)=sum(sum(grid(:,:,1)));
                                    protected_count(j,step)=sum(sum(grid(:,:,2)));
                                    dead_count(j,step)=sum(sum(grid(:,:,3)));
                                    reservoir_count(j,step)=sum(sum(grid(:,:,4)));
                                    virus_count(j,step)=sum(sum(grid(:,:,5)));
                                    ifn_count(j,step)=sum(sum(grid(:,:,6)));
                                    virusproda_cell(j,step)=sum(sum(grid(:,:,7)));
                                    virusprodb_cell(j,step)=sum(sum(grid(:,:,8)));
                                    ifn_cell(j,step)=sum(sum(grid(:,:,11)));
                                    virus_contact(j,step)=sum(sum(grid(:,:,18)));
                                    ifn_contact(j,step)=sum(sum(grid(:,:,19)));
                                    type1= grid(:,:,20)==1;
                                    type2= grid(:,:,20)==2;
                                    virus_type1(j,step) = sum(sum(type1));
                                    virus_type2(j,step) = sum(sum(type2));    
                                end
                            end

                            %filename = [num2str(xdistribution),'_virusprod',num2str(xvirus_prod),'_virusdiff',num2str(xvirus_diff),'_ifnprod',num2str(xifn_prod),'_ifnprob',num2str(ifn_prob(iifn_prob)),'_dim',num2str(grid_size(1)),'_sim',num2str(num_model_runs),'_step',num2str(num_steps)];
                            save(['./data/results_',num2str(filename),'.mat'])
                        end
                    end
                end
            end
        end
    end
end

%%

