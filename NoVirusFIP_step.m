    %% Model without individual viruses step function
%20210108 Pinghan Chu
%Does not consider viruses, but rather any infection that is successful
%within a time step immediately infects the cell
%expontial decay
%negative binomial 
%poission distribution
%parameter mu = 1
%parameter distribution =0.1
%@param grid nxnxm grid representation of the simulation 
%       where nxn is the square matrix and m is the number of states
%@param arr array of variables from the variable text file
%@return oGrid the newly updated grid
%20210112 Pinghan Chu
% The sum of virus_prod and ifn_prod is a constant.
% Tune the sum so that <ifn_prod> = 1,3,6
% Only 10% of cells can generate ifn_prod


function [output_grid,C,C2,C3,C4,C5] = NoVirusFIP_step(grid, arr)
    %disp(size(grid));
    infected=grid(:,:,1); % infected cell
    protected=grid(:,:,2); % protected cell
    dead=grid(:,:,3); % dead cell
    reservoir=grid(:,:,4); % reservoir cell
    
    viruscell=grid(:,:,5); % cells generating virus
    ifncell=grid(:,:,6); % cells generating ifn
    
    virus_proda = grid(:,:,7); % number of virus generated
    virus_prodb = grid(:,:,8); % number of virus generated
    virus_prod0a = grid(:,:,9); % original number of virus generated
    virus_prod0b = grid(:,:,10); % orignal number of virus generated
    A_virus_prod0a = reshape(virus_prod0a,10000,1);
    A_virus_prod0b = reshape(virus_prod0b,10000,1);
    ifn_prod=grid(:,:,11); % number of ifn generated
    
    infected_timer=grid(:,:,12); % exposure timer
    reservoir_timer=grid(:,:,13); % reservoir timer
    protected_timer=grid(:,:,14);
    dead_timer=grid(:,:,15);
    
    infected_lifespan=grid(:,:,16); %infected cell life time
    reservoir_lifespan=grid(:,:,17); %reservoir cell life time
    
    virus_contact=grid(:,:,18);
    ifn_contact=grid(:,:,19);
    
    virus_type=grid(:,:,20);
        
    virus_proda_mean = arr(1); % virus production mean
    virus_prodb_mean = arr(2); % virus production mean
    virus_diff = arr(3); % virus diffusion constant
    virus_prod_delay = arr(4); % virus production delay
    ifn_prod_mean = arr(5); % interferon production mean
    ifn_diff = arr(6); % interferon diffusion constant
    ifn_prod_delay = arr(7); % interferon production delay
    protected_lifespan = arr(8);% protected cells become susceptible cells
    dead_lifespan = arr(9); % dead cells are replaced by new susceptible cells
    prob_infect = arr(10); % probability of infection
    virus_reduction_factor=[arr(11)/100,1-arr(11)/100];
    lambda = 1 / prob_infect;
    
    s = size(grid);
    n = s(1);
    
    infected1=infected; % infected cell
    
    protected1=protected; % protected cell
    dead1=dead; % dead cell
    reservoir1=reservoir; % reservoir cell
    
    viruscell1=viruscell; % cells generating virus
    ifncell1=ifncell; % cells generating ifn
    
    virus_proda1 = virus_proda; % number of virus a generated
    virus_prodb1 = virus_prodb; % number of virus b generated
    ifn_prod1 = ifn_prod; % number of ifn generated
    
    infected_timer1=infected_timer; % exposure timer
    reservoir_timer1=reservoir_timer; % reservoir timer
    protected_timer1=protected_timer;
    dead_timer1=dead_timer;
    %disp(reservoir1);
    
    virus_contact1=zeros(n,n); % cell exposured to virus
    ifn_contact1=zeros(n,n); % cell exposured to ifn
    virus_type1=virus_type; % virus type (a or b)
    
    for i=1:n
        for j=1:n    
            virus_proda1(i,j)=virus_proda(i,j);
            virus_prodb1(i,j)=virus_prodb(i,j);
            ifn_prod1(i,j)=ifn_prod(i,j);
            %   timer
            if (reservoir(i,j)>0)
                reservoir_timer1(i,j)= reservoir_timer(i,j) + 1;
            end
            infected_timer1(i,j) = infected_timer(i,j) + infected(i,j);
            protected_timer1(i,j) = protected_timer(i,j) + protected(i,j);
            dead_timer1(i,j) = dead_timer(i,j) + dead(i,j);
            
            if(reservoir_timer(i,j)>reservoir_lifespan(i,j))
                if(reservoir(i,j)>0)
                    %disp([i,j,reservoir(i,j),reservoir_timer(i,j),reservoir_lifespan(i,j)]);
                    infected1(i,j)=1;
                    reservoir1(i,j)=0;
                    virus_type1(i,j)=reservoir(i,j);
                    reservoir_timer1(i,j)=0;  
                    if(protected(i,j)>0)
                    	if(reservoir(i,j)==1)
                            virus_proda1(i,j) = round(virus_reduction_factor(1)*virus_prod0a(i,j));
                        elseif reservoir(i,j)==2
                            virus_prodb1(i,j) = round(virus_reduction_factor(2)*virus_prod0b(i,j));
                        end
                    end
                end
            end
            if infected_timer(i,j)>infected_lifespan(i,j)
                dead1(i,j)=1;
                infected1(i,j) = 0;
                protected1(i,j)=0;
                viruscell1(i,j)=0;
                ifncell1(i,j)=0;
                virus_type1(i,j)=0;
                infected_timer1(i,j)=0;
                protected_timer1(i,j)=0;
            end  
            if protected_timer(i,j)>protected_lifespan
                protected1(i,j)=0;
                protected_timer1(i,j)=0;
                infected1(i,j) = 0;
                virus_type1(i,j)=0;
            end
            if dead_timer(i,j)>dead_lifespan
                dead1(i,j)=0;
                dead_timer1(i,j)=0;
                infected1(i,j) = 0;
                protected1(i,j)=0;
                infected_timer1(i,j)=0;
                virus_type1(i,j)=0;
                xa = ceil(rand()*10000);
                xb = ceil(rand()*10000);
                virus_proda1(i,j)=A_virus_prod0a(xa);
                virus_prodb1(i,j)=A_virus_prod0b(xb);
            end
                     

            if (dead(i,j) > 0)
                infected1(i,j)=0;
                protected1(i,j)=0;
                viruscell1(i,j)=0;
                ifncell1(i,j)=0;
                virus_type1(i,j)=0;               
                infected_timer1(i,j)=0;
                protected_timer1(i,j)=0;
            else
                xvirus_prod=0;
                if(virus_type(i,j)==1)
                    xvirus_prod = virus_proda1(i,j);
                elseif (virus_type(i,j)==2)
                    xvirus_prod = virus_prodb1(i,j);
                end
                if (xvirus_prod>0 && infected_timer(i,j) > virus_prod_delay) 
                    %disp([i,j,virus_prod(i,j),infected_timer(i,j),infected(i,j),virus_prod_delay]);
                    viruscell1(i,j)=1;
                    dt = -log(rand(xvirus_prod,1))*lambda;
                    r = abs(sqrt(2*virus_diff*dt).*sqrt(-2*log(rand(xvirus_prod,1))).*cos(2*3.1415926*rand(xvirus_prod,1))+0);
                    theta = 2*pi*rand(xvirus_prod,1);
                    x = r.*sin(theta);
                    y = r.*cos(theta);
                    x = round(x);
                    y = round(y);
                    for k = 1:round(xvirus_prod)
                        a = i+x(k);
                        b = j+y(k);

                        if (~((x(k) == 0 && y(k) == 0) || a <= 0 || b <= 0 || a > n || b > n))
                            virus_contact1(a,b) = virus_contact1(a,b)+1;
                           	if (dead(a,b)==0 && infected(a,b)==0) % assuming the virus cannot infect reservoir cell
                                if (protected(a,b) == 0)
                                    infected1(a,b)=1;
                                    virus_type1(a,b)=virus_type(i,j);
                                    reservoir1(a,b)=0;
                                else  
                                    if (virus_type(i,j)==2) 
                                        infprob=rand();
                                        if(infprob>0.5)
                                            infected1(a,b)=1;
                                            virus_type1(a,b)=virus_type(i,j);
                                            protected1(a,b)=0;
                                            protected_timer1(a,b)=0;
                                        end
                                    end
                                end
                            end    
                        end
                    end
                end
                
                 %if steps since exposure long enough and if not dead
                if (ifn_prod(i,j) > 0 && infected_timer(i,j) > ifn_prod_delay)
                    ifncell1(i,j) = 1;
                    
                    %math by Michael Lavigne
                    %dt = exprnd(lambda,ifn_prod,1);
                    %r = normrnd(0, sqrt(2*ifn_diff.*dt));
                    dt = -log(rand(ifn_prod(i,j),1))*lambda;

                    % r = abs(sqrt(2*ifn_diff*dt)*sqrt(-2*log(rand([1,1]))*cos(2*3.1415926*rand([1,1]))+0);
                    r= abs(sqrt(2*ifn_diff*dt).*sqrt(-2*log(rand(ifn_prod(i,j),1)).*cos(2*3.1415926*rand(ifn_prod(i,j),1))));
                    theta = 2*3.1415926*rand(ifn_prod(i,j),1);
                    %disp(size(dt));
                    %change r and theta to x and y
                    x = r.*sin(theta);
                    y = r.*cos(theta);
                    
                    %change to valid coordinates
                    x = round(x);
                    y = round(y);
                    
                    %disp([r,theta,x,y]);
                    %%%%%%%%%%%%%%%%%%%
                    for k = 1:ifn_prod(i,j)
                        a = i+x(k);
                        b = j+y(k);
                        %disp([k,i,j,a,b,x(k),y(k)]);
                        %boundary checking
                        %%%%%%%%%%%%%%%%%%%
                        if (~((x(k) == 0 && y(k) == 0) || a <= 0 || b <= 0 || a > n || b > n))
                            %%%%%% MICHAEL EDIT
                            
                            ifn_contact1(a,b) = ifn_contact1(a,b)+1;
                            
                            if (dead(a,b) == 0 && protected(a,b)==0) % if neither exposed nor infected,
                                if (infected(a,b) == 0)
                                    protected1(a,b) = 1;                     % becomes protected
                                    %reservoir1(a,b) = 0;
                                elseif (infected(a,b) > 0)
                                    if(virus_type(a,b)==1)
                                        virus_proda1(a,b) = round(virus_reduction_factor(1)*virus_prod0a(a,b));
                                    elseif virus_type(a,b)==2
                                        virus_prodb1(a,b) = round(virus_reduction_factor(2)*virus_prod0b(a,b));
                                    end
                                end
                            end
                        %%%%%%%%%%%%%%%%%%%
                        end
                    %%%%%%%%%%%%%%%%%%%
                    end
                %%%%%%%%%%%%%%%%%%%
                end   
            end
        end
    end

    output_grid(:,:,1)=infected1;
    output_grid(:,:,2)=protected1;
    output_grid(:,:,3)=dead1;
    output_grid(:,:,4)=reservoir1;
    output_grid(:,:,5)=viruscell1;
    output_grid(:,:,6)=ifncell1;
    output_grid(:,:,7)=virus_proda1;
    output_grid(:,:,8)=virus_prodb1;
    output_grid(:,:,9)=virus_prod0a;
    output_grid(:,:,10)=virus_prod0b;
    output_grid(:,:,11)=ifn_prod1;
    output_grid(:,:,12)=infected_timer1;
    output_grid(:,:,13)=reservoir_timer1;
    output_grid(:,:,14)=protected_timer1;
    output_grid(:,:,15)=dead_timer1;
    output_grid(:,:,16)=infected_lifespan;
    output_grid(:,:,17)=reservoir_lifespan;
    output_grid(:,:,18)=virus_contact1;
    output_grid(:,:,19)=ifn_contact1;
    output_grid(:,:,20)=virus_type1;

    C=zeros(n,n,3);
    C2=zeros(n,n,3);
    C3=zeros(n,n,3);
    C4=zeros(n,n,3);
    C5=zeros(n,n,3);
    % true colors, 3 layers (red, green, blue)
    C(:,:,1)= infected1;                                          %set red
    C(:,:,2)= infected1 - (infected_timer1>virus_prod_delay);      %set green
    C(:,:,3)= protected1;                                         %set blue
    % if C(1)=1 and C(2)=1 -> yellow; it means exposure_timer>
    % virus_prod_delay is wrong. 
    C2(:,:,1) = virus_contact1; %red
    C2(:,:,2) = ifn_contact1; %green
    C2(:,:,3) = 0;
    C3(:,:,1) = ifncell1;
    C3(:,:,2) = ifn_contact1;
    C3(:,:,3) = 0;
    C4(:,:,1) = viruscell1;
    C4(:,:,2) = virus_contact1;
    C4(:,:,3) = 0;
    C5(:,:,1) = reservoir1;
    C5(:,:,2) = virus_type1==1;
    C5(:,:,3) = virus_type1==2;
    targets=(C(:,:,1)+ C(:,:,2)+ C(:,:,3)==0); 
 
    for i=1:3
        C(:,:,i)=C(:,:,i)+targets;
    end
 
    for i=1:n
        for j=1:n
            if dead(i,j)==1
                C(i,j,:)=0.5; % grey if a cell is dead, all colors are set 0.5 (0 will be black, 1 will be white)
            end
        end
    end   
end
