

% ************************************************** Input Parameters *****************************************

clc
clear all


I_NP = 50;                                           % Population Size
I_D = 24;                                            % Number of Variables
F_weight = 0.80;                                     % DE-stepsize F_weight ex [0, 2]
F_CR = 0.20;                                         % crossover probabililty constant ex [0, 1]

MinParValue = [0];                                   % Lower bound for Generator Capacity
MaxParValue = [66];                                  % Upper bound for Generator Capacity

FVr_minbound = [MinParValue];                        % Lower bound for Generator Capacity
FVr_maxbound = [MaxParValue];                        % Upper bound for Generator Capacity


PD = 636.42;                                         % Maximum Load Demand Capacity

I_itermax = 5000;                                   % Maximum Iteration

I_strategy = 4;                                      % I_strategy     1 --> DE/rand/1:


if (I_NP < 5)
    I_NP=5;
    fprintf(1,' I_NP increased to minimal value 5\n');
end


if ((F_CR < 0) | (F_CR > 1))
    F_CR=0.5;
    fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
end


if (I_itermax <= 0)
    I_itermax = 200;
    fprintf(1,'I_itermax should be > 0; set to default value 200\n');
end

% ******************************************** initialization of population ****************************************

FM_pop = zeros(I_NP,I_D);                      % Initialize FM_pop to gain speed


d=0;
i = 1;
while ( i ~= (I_NP + 1))
    
    j=1;
    while (j ~= I_D)
        d(i,j) = (FVr_minbound(1,1) + (FVr_maxbound(1,1) - FVr_minbound(1,1)) * rand);
        j=j+1;
    end
    
    x = sum((d(i,:))')';
    y = PD - x;
    
    if (y <= FVr_maxbound(1,1) & y >= FVr_minbound(1,1))
        d(i,j)=y;
        FM_pop(i,:) = d(i,:);
        i=i+1
    end
end
FM_pop
% pause

% ********************************************** End of the update ****************************************************
FM_popold     = zeros(size(FM_pop));             % toggle population
FVr_bestmem   = zeros(1,I_D);                    % best population member ever
FVr_bestmemit = zeros(1,I_D);                    % best population member in iteration
I_nfeval      = 0;                               % number of function evaluations
% ********************************************** calculation of Cost Function *****************************************

inelasticload
ineL = transpose(inelasticloade(:,1));
Tr = transpose(inelasticloade(:,2));


for j=1:I_NP
    
    F1=0;
    for i=1:I_D
        F1=F1+(Tr(i)*(FM_pop(j,i)+ineL(i)));
    end
    
    Dejed(j,1) = F1;

end
% pause

% ******************************************* Evalution of Best Cost population Set*********************
[S_bestval I_best_index] = min(Dejed);
S_bestvalit   = S_bestval;                    % best value of current iteration
FVr_bestmemit = FM_pop(I_best_index,:);       % best member of current iteration
FVr_bestmem = FVr_bestmemit;                  % best member ever

% ******************************************* Crossover Operation **************************************

FM_pm1   = zeros(I_NP,I_D);                   % initialize population matrix 1
FM_pm2   = zeros(I_NP,I_D);                   % initialize population matrix 2
FM_pm3   = zeros(I_NP,I_D);                   % initialize population matrix 3
FM_pm4   = zeros(I_NP,I_D);                   % initialize population matrix 4
FM_pm5   = zeros(I_NP,I_D);                   % initialize population matrix 5
FM_bm    = zeros(I_NP,I_D);                   % initialize FVr_bestmember  matrix
FM_ui    = zeros(I_NP,I_D);                   % intermediate population of perturbed vectors
FM_mui   = zeros(I_NP,I_D);                   % mask for intermediate population
FM_mpo   = zeros(I_NP,I_D);                   % mask for old population
FVr_rot  = (0:1:I_NP-1);                      % rotating index array (size I_NP)
FVr_rotd = (0:1:I_D-1);                       % rotating index array (size I_D)
FVr_rt   = zeros(I_NP);                       % another rotating index array
FVr_rtd  = zeros(I_D);                        % rotating index array for exponential crossover
FVr_a1   = zeros(I_NP);                       % index array
FVr_a2   = zeros(I_NP);                       % index array
FVr_a3   = zeros(I_NP);                       % index array
FVr_a4   = zeros(I_NP);                       % index array
FVr_a5   = zeros(I_NP);                       % index array
FVr_ind  = zeros(4);

FM_meanv = ones(I_NP,I_D);

disp(sprintf('Iterations\t\tS_bestval\t\t\t\t\t\t\t\tFVr_bestmemit\t\t\t\t\tLossF'))

tic

I_iter = 1;
while (I_iter ~= I_itermax+1)
    
    FM_popold = FM_pop;                              % save the old population
    S_struct.FM_pop = FM_pop;
    S_struct.FVr_bestmem = FVr_bestmem;
    
    y3 = 0;
    i60 = 1;
    while (i60~=I_NP+1)
        
        FVr_ind = randperm(4);
        
        FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
        FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);     % rotate indices by ind(1) positions
        FVr_a2  = FVr_a1(FVr_rt+1);                 % rotate vector locations
        FVr_rt  = rem(FVr_rot+FVr_ind(2),I_NP);
        FVr_a3  = FVr_a2(FVr_rt+1);
        FVr_rt  = rem(FVr_rot+FVr_ind(3),I_NP);
        FVr_a4  = FVr_a3(FVr_rt+1);
        FVr_rt  = rem(FVr_rot+FVr_ind(4),I_NP);
        FVr_a5  = FVr_a4(FVr_rt+1);
        
        FM_pm1 = FM_popold(FVr_a1,:);               % shuffled population 1
        FM_pm2 = FM_popold(FVr_a2,:);               % shuffled population 2
        FM_pm3 = FM_popold(FVr_a3,:);               % shuffled population 3
        FM_pm4 = FM_popold(FVr_a4,:);               % shuffled population 4
        FM_pm5 = FM_popold(FVr_a5,:);               % shuffled population 5
        
        for k=1:I_NP                                % population filled with the best member
            FM_bm(k,:) = FVr_bestmemit;             % of the last iteration
        end
        
        FM_mui = rand(I_NP,I_D) < F_CR;             % all random numbers < F_CR are 1, 0 otherwise
        
        FM_mpo = FM_mui < 0.5;                      % inverse mask to FM_mui
        
        
        Swarmxx(i60,:) = FM_ui(i60,:);
        
        if (I_strategy == 1)                        % DE/rand/1
            Swarmxx(i60,:) = FM_pm3(i60,:) + F_weight*(FM_pm1(i60,:) - FM_pm2(i60,:));             % differential variation
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);      % crossover
            FM_origin(i60,:) = FM_pm3(i60,:);
            
        elseif (I_strategy == 2)                    % DE/local-to-best/1
            Swarmxx(i60,:) = FM_popold(i60,:) + F_weight*(FM_bm(i60,:)-FM_popold(i60,:)) + F_weight*(FM_pm1(i60,:) - FM_pm2(i60,:));
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);
            FM_origin(i60,:) = FM_popold(i60,:);
            
        elseif (I_strategy == 3)                    % DE/best/1 with jitter
            Swarmxx(i60,:) = FM_bm(i60,:) + (FM_pm1(i60,:) - FM_pm2(i60,:)).*((1-0.999)*rand(1,I_D)+F_weight);
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);
            FM_origin(i60,:) = FM_bm(i60,:);
            
        elseif (I_strategy == 4)                    % DE/rand/1 with per-vector-dither
            f1 = ((1-F_weight)*rand(I_NP,1)+F_weight);
            
            for k=1:I_D
                FM_pm5(:,k)=f1;
            end
            
            Swarmxx(i60,:) = FM_pm3(i60,:) + (FM_pm1(i60,:) - FM_pm2(i60,:)).*FM_pm5(i60,:);      % differential variation
            FM_origin(i60,:) = FM_pm3(i60,:);
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);     % crossover
            
        elseif (I_strategy == 5)                    % DE/rand/1 with per-vector-dither
            f1 = ((1-F_weight)*rand+F_weight);
            Swarmxx(i60,:) = FM_pm3(i60,:) + (FM_pm1(i60,:) - FM_pm2(i60,:))*f1;                  % differential variation
            FM_origin(i60,:) = FM_pm3(i60,:);
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);     % crossover
            
        else                                        % either-or-algorithm
            if (rand < 0.5)                        % Pmu = 0.5
                Swarmxx(i60,:) = FM_pm3(i60,:) + F_weight*(FM_pm1(i60,:) - FM_pm2(i60,:));% differential variation
                FM_origin(i60,:) = FM_pm3(i60,:);
            else                                    % use F-K-Rule: K = 0.5(F+1)
                
                Swarmxx(i60,:) = FM_pm3(i60,:) + 0.5*(F_weight+1.0)*(FM_pm1(i60,:) + FM_pm2(i60,:) - 2*FM_pm3(i60,:));
                FM_origin(i60,:) = FM_pm3(i60,:);
            end
            
            Swarmxx(i60,:) = FM_popold(i60,:).*FM_mpo(i60,:) + Swarmxx(i60,:).*FM_mui(i60,:);      % crossover
        end
        
        % **************************************** boundary constraints via bounce back **************************************
        % ***************************** Select which vectors are allowed to enter the new population *************************
        
        % ************************************************* update this portion **********************************************
        
        j60=1;
        while ( j60 ~= I_D + 1)
            if (Swarmxx(i60,j60)>= FVr_maxbound(1,1))
                Swarmxx(i60,j60) = FVr_maxbound(1,1);
                
            elseif (Swarmxx(i60,j60)<= FVr_minbound(1,1))
                Swarmxx(i60,j60) = FVr_minbound(1,1);
            end
            j60=j60+1;
        end
        
        x3k = sum((Swarmxx(i60,:))')' - Swarmxx(i60,I_D);
        y3(i60,1) = PD - x3k;
        
        if ((y3(i60,1)<=FVr_maxbound(1,1)) & (y3(i60,1)>=FVr_minbound(1,1)))
            Swarmxx(i60,I_D)=y3(i60,1);
            i60 = i60 + 1;
        end
    end
    
    % ********************************************** End of update process *****************************************
    FM_ui = Swarmxx;
    pop11 = FM_ui;
    % ******************************************** Cost Function Calculation ***************************************
    
    for j=1:I_NP
        
        F2=0;
        for i=1:I_D
            F2=F2+(Tr(i)*(pop11(j,i)+ineL(i)));
        end
        
        S_tempval(j,1) = F2;
        
    end
    
    % ************************************************** End of update process  ***********************************************
    
    for k=1:I_NP
        if ( S_tempval(k,1) <= Dejed(k,1))
            FM_pop(k,:) = pop11(k,:);                    % replace old vector with new one (for new iteration)
            Dejed(k,1)  = S_tempval(k,1);                % save value in "cost array"
            
            % ****************************************** we update S_bestval only in case of success to save time ******************************************
            
            if ( S_tempval(k,1) <= S_bestval )
                S_bestval = S_tempval(k,1);              % new best value
                FVr_bestmem = pop11(k,:);                % new best parameter vector ever
            end
        end
    end
    
    FVr_bestmemit = FVr_bestmem;                         % freeze the best member of this iteration for the coming
    
    % iteration. this is needed for some of the strategies.
    % xd(e10,:) = FVr_bestmemit;
    % fxmin(e10,:) = S_bestval;
    
    % ************************************************** Output section *******************************************************
    XX(I_iter,:) = S_bestval;
    YY(I_iter,:) = I_iter;
    
    disp(sprintf('%7d\t\t\t%.12f\t\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f\t\t%.10f', I_iter, S_bestval, FVr_bestmemit))
    I_iter = I_iter + 1
    
end

toc

% ****************************************************** Output Plot ************************************************
plot(YY(:,1),XX(:,1))
% title([' The Model of Demand Side versus number of iteration ', num2str(S_bestval)]);
xlabel('Iteration');
ylabel('Objective Value');
legend('DE')
