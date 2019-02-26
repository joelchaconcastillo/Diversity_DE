function mvmo(fhd,iii,jjj,kkk,args)
    global proc
    global ps %disc
    global parameter table
    global changed_best
    global mod_Method
    persistent amax
    
name_dat_file=[pwd '\Parameter_MVMO-SH_' num2str(iii) '_D' num2str(ps.D) '.txt'] ;
%  the parameter files must be located in the current folder
%  for using a separate folder, change the command line above

AAAA2 = exist(name_dat_file,'file')  ;

if AAAA2~=0
    fidtunp=fopen(name_dat_file,'r');
    param_set=zeros(17,1);
    flag=0; 
    while flag<17
        flag=flag+1;
        tline=fgets(fidtunp);
        param_set(flag)=str2num(tline(17:31));
    end
    fclose(fidtunp);
    flagfile=1;
end         
%
    proc.finish=0;
    proc.i_eval=0;
	proc.last_improvement=1;

% only the following parameters are required from the user    
   switch iii   %  iii represents the function number
%  the parameters below are tuned for dimension D10
%  these parametrs will be used for all dimensions unless separate parameter files are provided
   case {1}
    parameter.n_par=10;                   % Number of particles  
    parameter.n_tosave=5 ;                % Archive size
    parameter.fs_factor_start=1;          % Initial fs-factor 
    parameter.fs_factor_end=20;           % Final fs-factor
    parameter.local_prob= 100;            % Local search probability factor
    min_eval_LS =round(0.25*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.8;         % Initial portion of good particles    
    parameter.ratio_gute_min=0.5 ;        % Final portion of good particles
    parameter.n_random_ini =5   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=5  ;          % Final number for variables selected for mutation
    local_i_max= 5     ;                  % max number of local search runs allowed
    alpha=1d0;                            % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {2}
    parameter.n_par=10;                   % Number of particles  
    parameter.n_tosave=5 ;                % Archive size
    parameter.fs_factor_start=1;          % Initial fs-factor 
    parameter.fs_factor_end=20;           % Final fs-factor
    parameter.local_prob= 100;            % Local search probability factor
    min_eval_LS =round(0.25*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.8;         % Initial portion of good particles    
    parameter.ratio_gute_min=0.5 ;        % Final portion of good particles
    parameter.n_random_ini =5   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=5  ;          % Final number for variables selected for mutation
    local_i_max= 5     ;                  % max number of local search runs allowed
    alpha=1d0;                            % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {3}
    parameter.n_par=103;                  % Number of particles  
    parameter.n_tosave=7 ;                % Archive size
    parameter.fs_factor_start=0.93186;    % Initial fs-factor 
    parameter.fs_factor_end=25.731;       % Final fs-factor
    parameter.local_prob= 100;            % Local search probability factor
    min_eval_LS =round(0.44521*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.91534    ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.90578    ; % Final portion of good particles
    parameter.n_random_ini =3   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=2  ;          % Final number for variables selected for mutation
    local_i_max= 5     ;                  % max number of local search runs allowed
    alpha=0.43719E+01;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {4}
    parameter.n_par=0.21230E+04;          % Number of particles  
    parameter.n_tosave=3 ;                % Archive size
    parameter.fs_factor_start=0.37449E+01 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.23913E+02;       % Final fs-factor
    parameter.local_prob=0.57590E+02;     % Local search probability factor
    min_eval_LS =round(0.36139E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.80189E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.86668E+00  ; % Final portion of good particles
    parameter.n_random_ini =8   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=6  ;          % Final number for variables selected for mutation
    local_i_max= 0.44000E+02     ;        % max number of local search runs allowed
    alpha=0.66869E+01;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {5}
    parameter.n_par=0.10900E+03;          % Number of particles  
    parameter.n_tosave=4 ;                % Archive size
    parameter.fs_factor_start=0.60309E+00 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.38264E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.36542E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.82783E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.35316E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.20000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=0.20000E+01  ;          % Final number for variables selected for mutation
    local_i_max= 0.63000E+02     ;        % max number of local search runs allowed
    alpha=0.51665E+01;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {6}
    parameter.n_par=0.10000E+03;          % Number of particles  
    parameter.n_tosave= 0.50000E+01 ;                % Archive size
    parameter.fs_factor_start=0.98734E+00 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.46405E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.35199E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.82733E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.45667E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.10000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.10000E+01  ;          % Final number for variables selected for mutation
    local_i_max= 0.10500E+03     ;        % max number of local search runs allowed
    alpha=0.85784E+01;                    % alpha-factor
     derivative='forward';                % Method for num. derivative calculation
%   derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {7}
    parameter.n_par=0.74000E+02;          % Number of particles  
    parameter.n_tosave= 0.50000E+01 ;                % Archive size
    parameter.fs_factor_start=0.71118E+00 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.37401E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.26380E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.84782E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.32648E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.10000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=0.20000E+01  ;          % Final number for variables selected for mutation
    local_i_max=  0.10000E+01     ;        % max number of local search runs allowed
    alpha=0.26386E+00;                    % alpha-factor
     derivative='forward';                % Method for num. derivative calculation
%   derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {8}
     parameter.n_par=0.11000E+03 ;          % Number of particles  
    parameter.n_tosave=0.30000E+01 ;                % Archive size
    parameter.fs_factor_start=0.90485E+00 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.34906E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.35614E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.87272E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.53172E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.30000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=0.10000E+01  ;          % Final number for variables selected for mutation
    local_i_max= 0.21000E+02     ;        % max number of local search runs allowed
    alpha=0.68679E+01;                    % alpha-factor
     derivative='forward';                % Method for num. derivative calculation
%   derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {9}
    parameter.n_par=0.20000E+03 ;          % Number of particles  
    parameter.n_tosave=0.50000E+01 ;                % Archive size
    parameter.fs_factor_start=0.50000E+01 ;    % Initial fs-factor 
    parameter.fs_factor_end=0.33492E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.38645E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.89347E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.80000E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.50000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last=0.50000E+01  ;          % Final number for variables selected for mutation
    local_i_max=  0.21000E+02     ;        % max number of local search runs allowed
    alpha=0.10000E+01;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {10}
    parameter.n_par=0.11100E+03 ;          % Number of particles  
    parameter.n_tosave=0.40000E+01 ;                % Archive size
    parameter.fs_factor_start=0.69667E+00  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.28481E+02;       % Final fs-factor
    parameter.local_prob=0.92240E+02;     % Local search probability factor
    min_eval_LS =round(0.44863E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.84790E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.35176E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.10000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.10000E+01  ;          % Final number for variables selected for mutation
    local_i_max=   0.14700E+03     ;        % max number of local search runs allowed
    alpha=0.37790E+01 ;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {11}
    parameter.n_par=0.26400E+03 ;          % Number of particles  
    parameter.n_tosave=0.60000E+01 ;                % Archive size
    parameter.fs_factor_start=0.26952E+00  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.47129E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.45156E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.91207E+00 ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.99158E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.90000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.30000E+01   ;          % Final number for variables selected for mutation
    local_i_max= 0.10000E+02     ;        % max number of local search runs allowed
    alpha=0.24451E+00 ;                    % alpha-factor
     derivative='forward';                % Method for num. derivative calculation
%   derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {12}
    parameter.n_par=0.12200E+03 ;          % Number of particles  
    parameter.n_tosave=0.40000E+01 ;                % Archive size
    parameter.fs_factor_start=0.18286E+01  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.43774E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.38340E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.89281E+00 ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.85773E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.11000E+02   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.13000E+02  ;          % Final number for variables selected for mutation
    local_i_max= 0.12600E+03     ;        % max number of local search runs allowed
    alpha=0.46777E+01 ;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {13}
    parameter.n_par=0.13900E+03 ;          % Number of particles  
    parameter.n_tosave=0.30000E+01 ;                % Archive size
    parameter.fs_factor_start=0.91627E+00  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.40912E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.68484E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.90695E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.50818E+00  ; % Final portion of good particles
    parameter.n_random_ini =0.10000E+01   ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.10000E+01  ;          % Final number for variables selected for mutation
    local_i_max=0.42000E+02     ;        % max number of local search runs allowed
    alpha=0.57327E+01 ;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {14}
    parameter.n_par= 0.70400E+03;          % Number of particles  
    parameter.n_tosave=0.90000E+01 ;                % Archive size
    parameter.fs_factor_start=0.32208E+01  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.17945E+02;       % Final fs-factor
    parameter.local_prob=0.14878E+02;     % Local search probability factor
    min_eval_LS =round(0.19684E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.98484E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.66685E+00  ; % Final portion of good particles
    parameter.n_random_ini = 0.10000E+02    ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.70000E+01  ;          % Final number for variables selected for mutation
    local_i_max=0.55000E+02     ;        % max number of local search runs allowed
    alpha=0.40232E+01 ;                    % alpha-factor
%    derivative='forward';                % Method for num. derivative calculation
    derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
    case {15}
    parameter.n_par=0.15000E+03 ;          % Number of particles  
    parameter.n_tosave=0.50000E+01 ;                % Archive size
    parameter.fs_factor_start=0.34119E+01  ;    % Initial fs-factor 
    parameter.fs_factor_end=0.14037E+02;       % Final fs-factor
    parameter.local_prob=0.10000E+03;     % Local search probability factor
    min_eval_LS =round(0.33857E+00*proc.n_eval); % Begin of local search (number in the expression) 
    parameter.ratio_gute_max=0.89109E+00  ; % Initial portion of good particles    
    parameter.ratio_gute_min=0.80000E+00  ; % Final portion of good particles
    parameter.n_random_ini = 0.14000E+02    ;         % Initial number for variables selected for mutation 
    parameter.n_random_last= 0.50000E+01  ;          % Final number for variables selected for mutation
    local_i_max=0.16800E+03     ;        % max number of local search runs allowed
    alpha=0.13886E+02 ;                    % alpha-factor
     derivative='forward';                % Method for num. derivative calculation
%   derivative='central';                 % Method for num. derivative calculation
    optimmethod='interior-point';         % Optimization method used for local search
%    optimmethod='sqp';                   % Optimization method used for local search
   end %case

    % end of parameter input

    
     if AAAA2~=0 & flagfile~=0  %  MVMO tuning parameters are taken from the corresponding parameter file  
        parameter.n_par=param_set(1)   ;
        parameter.n_tosave=param_set(4) ;
        parameter.fs_factor_start=param_set(7) ;
        parameter.fs_factor_end=param_set(8) ;
        parameter.ratio_gute_max=param_set(2) ; 
        parameter.ratio_gute_min=param_set(3) ;
        parameter.n_random_ini =param_set(5) ; 
        parameter.n_random_last=param_set(6);
        parameter.local_prob=param_set(9);
        min_eval_LS =round(param_set(10)*proc.n_eval); 
        local_i_max=param_set(11);
        alpha=param_set(12);
        method_deriv=param_set(13);
        if method_deriv==1
            derivative='forward';
        else
            derivative='central';
        end
        optimmethod='interior-point'; %  always used in the matlab version
    end

    parameter.scaling=(ps.x_max-ps.x_min);
    %%----------------- Create initial random population ------------------   
    xx=zeros(parameter.n_par,ps.D);
    x_norm=xx;
    for iijj=1:parameter.n_par %Start initialization: x_normalized
          for jjkk=1:ps.D
              xx(iijj,jjkk)=ps.x_min (jjkk) + rand(1)*(ps.x_max(jjkk)-ps.x_min(jjkk));
          end
        x_norm(iijj,:)=(xx(iijj,:)-ps.x_min)./parameter.scaling;
    end % End initialization
    
    x_normalized=x_norm;
    
    %% ------------------ Initialize control parameters --------------------
    n_eval = proc.n_eval;
    n_par = parameter.n_par; 
    D=ps.D ; 
    n_to_save = parameter.n_tosave; 
    fs_factor_start = parameter.fs_factor_start;
    fs_factor_end = parameter.fs_factor_end;
    Shape_dyn = ones(n_par,D);%
    Shape_dyn(:,:) =1.d0;%
    IX=ones(n_par); 
    shape = zeros(n_par,D);
    izm = zeros(1,n_par);
    for i=1:n_par
     IX(i)=i;
    end
%    local_search0_percentage=parameter.local_prob;
    ratio_gute_min=parameter.ratio_gute_min;
    ratio_gute_max=parameter.ratio_gute_max;
    n_randomly_ini = parameter.n_random_ini;
    n_randomly_last = parameter.n_random_last;
    local_i=0;    % local seurch runs counter
    local_search_noch_nicht=false;
    local_i_selected=0;
    mod_Method=3;
    l_vari=n_to_save  ; % can be changed
      if  l_vari<2 || l_vari > n_to_save 
        l_vari=n_to_save ;
      end  

    
    %% --------------------- Data structure for the table ---------------------
    table.bests = NaN*zeros(parameter.n_tosave,ps.D,n_par);
    table.fitness = Inf*ones(parameter.n_tosave,1,n_par);
    table.objective = Inf*ones(parameter.n_tosave,1,n_par);
    table.feasibility = zeros(parameter.n_tosave,1,n_par);
    
    %% ----------------------------- Mapping ----------------------------------
    local_search=zeros(n_par);
    meann = x_normalized;
    meann(:,:)= 0.5d0  ;
    meann_app = meann;
   %% ------------------------ Variable selection ----------------------------
    izz = zeros(1,n_par);
    considered = true(ps.D);
    probab=ones(n_par,ps.D);
    values_noneq = zeros(n_to_save,1);
    
   %% ---------------------- Checking correctness ------------------------   
   
    if (n_randomly_last<1)
        n_randomly_last=1;
    end
    
    if (n_randomly_last>ps.D)
        n_randomly_last=ps.D;
    end
    
    if (n_randomly_ini>ps.D)
        n_randomly_ini=ps.D;
    end  
    if (n_eval<=0)
        n_eval=100000d0;
    end

    if (fs_factor_start<=0.0)
        fs_factor_start=1d0;
    end
    
    if (fs_factor_end<=0)
        fs_factor_end=1d0;
    end

    fs_factor0=fs_factor_start;
    
    yes_n_randomly=true;
    if (n_randomly_ini==n_randomly_last)
        n_randomly_n=n_randomly_last;
        yes_n_randomly=false;
    end
    
    if (n_to_save<=1)
        n_to_save=2d0;
    end
    
    yes_fs_factor=true;
    if (fs_factor_start==fs_factor_end)
        yes_fs_factor=false;
    end

%% ----------------------------- Counters ---------------------------------

    no_in = zeros(1,n_par);
    no_inin = zeros(1,n_par);
    amax=0.d0;
    
    local_search0=parameter.local_prob/100; % Probability of local search (percentage / number of optimization variables)
    goodbad=ones(n_par);
    firsttime=true;
    

    delta_nrandomly=n_randomly_ini-n_randomly_last ; 
    A=zeros(n_par,1); 
%    local_search_called=0;
   
  while 1       
        %Evaluating the particles.....
            ff=proc.i_eval/n_eval   ;
            ff2=ff*ff;
            vvqq=10.d0^(-(5.3d0+ff*5.3d0));
            %Determining the relative number of particles belonging to the group of
            %good particles
            border_gute0=ratio_gute_max-ff*(ratio_gute_max-ratio_gute_min);
            border_gute=(n_par*border_gute0) ;
            border_gute=round(normalfunc(border_gute , border_gute*0.15d0 ))   ;

             if border_gute < 3 && n_par > 3
                 border_gute=3;
             end
            if border_gute > n_par
                border_gute=n_par ;
             end                     
            %Selecting the subset of variables to be mutated
            if yes_n_randomly
                n_randomly_X=round(n_randomly_ini-ff*delta_nrandomly) ;
                n_randomly_n=round(n_randomly_last+rand(1)*(n_randomly_X-n_randomly_last)) ;
             end
            %Quadratic variation of fs_factor0
            if yes_fs_factor
                fs_factor0=fs_factor_start+ff2*(fs_factor_end-fs_factor_start)  ;
            end
            
     if local_search0 > 0.d0 && proc.i_eval > min_eval_LS && ~local_search_noch_nicht
         for iuz=1:border_gute
          if  rand <= local_search0 && local_i_selected < local_i_max   
           if  rand(1) <= local_search0  
             local_search(IX(iuz))=true;
             local_i_selected=local_i_selected+1;
           if  local_i_selected >= local_i_max 
              local_search_noch_nicht=true;
              break 
           end
          end
          end
         end   
     end  
      
     ipx=0 ;
     while ipx  < n_par 
         ipx=ipx+1 ;
         ipp=IX(ipx);

             
% -------------MAPPING----------------
%   only for not local search runs
%      considered(1:nDim)=.true.
       if no_inin(ipp) >= 1    
        if ~local_search(ipp) 
         considered(1:D)=false;
         VariableSelect1(n_randomly_n); % Call random variable selection strategy
            if no_inin(ipp) < l_vari 
                considered(1:D)=true;
            end
         for ivar=1:D
             if considered(ivar) 
       %                 x_normalized(ipp,ivar) = rand(1)  ;    
                  
                    if rand(1) > ff2    
                      x_normalized(ipp,ivar) = rand(1)  ;    
                    else
                      x_normalized(ipp,ivar) = -1.  ; 
                       while (x_normalized(ipp,ivar) > 1.d0 || x_normalized(ipp,ivar) < 0.d0) 
                            x_normalized(ipp,ivar)=normalfunc(0.5d0 , (1.d0-ff2)*(1.d0-ff2)+1.d-3 )  ;     
                       end  
                    end  
                    if shape(ipp,ivar) > 0.d0 
                      sss1 = shape(ipp,ivar);
                      if (sss1 > Shape_dyn(ipp,ivar))
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)*1.1d0  ; 
                      else
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)/1.1d0  ; 
                      end
                        if  Shape_dyn(ipp,ivar) > sss1
                            grosser=Shape_dyn(ipp,ivar);
                            kleiner=sss1;
                        else
                            kleiner=Shape_dyn(ipp,ivar);
                            grosser=sss1;
                        end   
 %                       if   meann_app(ipp,ivar) > 0.5d0    
                        if   rand(1) > 0.5d0    
                                 sss1=   grosser;
                                 sss2=   kleiner;
                        else
                                 sss1=   kleiner ;
                                 sss2=   grosser ;
                        end  
                         fs_factor=fs_factor0*(1.0d0 + 4.d0*(rand(1)-0.2)) ; 
                         sss1=sss1*fs_factor;   
                         sss2=sss2*fs_factor;  
                      x_normalized(ipp,ivar)=h_function(meann_app(ipp,ivar),sss1,sss2,x_normalized(ipp,ivar)); %Mapping function
                    end
             end
         end
       end  
       end
% -------------END of MAPPING----------------
%Denormalize to the original [min, max] range 
              x_normalized(ipp,:) = ps.x_min+parameter.scaling.* x_normalized(ipp,:);
             if  local_search(ipp)   
                [msgstr, msgid] = lastwarn ;
                TFrcond = strcmp('MATLAB:nearlySingularMatrix',msgid); % Only informative from 'fmincon' function 
                if TFrcond~=0
                    rcond_value0=str2num(msgstr(81:end-1));
                end
               
 %proc.i_eval
                [ffx,oox,~,x_normalized(ipp,:),FEVALS] = LocalSearchMVMOSH(fhd,iii,jjj,kkk,args,x_normalized(ipp,:),proc.n_eval-proc.i_eval); %Local search
                local_i=local_i+1  ; 
                if  local_i >= local_i_max 
                  local_search0=0;
                  local_search(:)=0 ;
                  mod_Method=2;

                end
                 xt.fitness=ffx ;  % Constraint handling outsourced so far
                 xt.objective=oox;
                 if xt.fitness==xt.objective
                    xt.feasibility= true; 
                 else
                    xt.feasibility= false;
                 end

                [~, msgid1] = lastwarn;
                TFrcond1 = strcmp('MATLAB:nearlySingularMatrix',msgid1);
             else
              %   ff
                [ffx,oox,~,x_normalized(ipp,:)]=feval(fhd,iii,jjj,kkk,args,x_normalized(ipp,:)); %Problem evaluation
                xt.fitness=ffx ; % Constraint handling outsourced so far   
                xt.objective=oox;
                if xt.fitness==xt.objective
                    xt.feasibility= true; %
                else
                    xt.feasibility= false;
                end
            end
            
            if proc.finish
                return;
            end
            
            x_normalized(ipp,:) = (x_normalized(ipp,:)-ps.x_min)./parameter.scaling;
            irandom=ipp;
            meann_app(ipp,1:ps.D) = meann(ipp,1:ps.D);
            if  border_gute > 2 
            while ipp==irandom 
             irandom=round(rand(1)*(border_gute-1))+ 1;
             irandom=IX(irandom) ;
            end
             meann_app(ipp,1:ps.D) =meann(irandom,1:ps.D) ;     
            end
             % Store the n-best solution corresponding to the corresponding particle's archive
%
           Fill_solution_archive(vvqq);
           
             if  no_inin(ipp) > l_vari        
            %Determining the proportion of good particles
                   if n_par > 1
                    if changed_best || firsttime
                     A(1:n_par)=table.fitness(1,:,1:n_par);
                    if  firsttime
                      amax=max(A);
                      firsttime=false  ;
                    end
                    for ia=1:n_par
                      if ~table.feasibility(1,:,ia)
                        A(ia)=A(ia)+amax;
                      end
                    end 
                    [~,IX] = sort(A);
                    end
                    goodbad(1:n_par)=0;
                    goodbad(IX(1:border_gute))=1;
                   end
                   if n_par < 4 
                    goodbad(ipp)=1   ;
                   end
                    
                    %Multi-parent strategy for bad particles  
                 if ~goodbad(ipp)    
                        [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;  
                        bbb=1.1d0+(rand(1)-0.5d0)*2.0d0  ; 
                        beta1 = alpha*3.d0*bbb*((1.d0+2.5*ff2)*rand(1) - (1.d0-ff2)*0.8d0 )  ;  
                        beta10=0.d0     ;
                        while beta1 ~= beta10
                          beta10=beta1 ;
                          for jx=1:ps.D
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                    if abs(ccc) >  1.d-8  
                                     x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;   
                                    else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)    ;
                                    end
                                 if table.bests(1,jx,bestp) <= 0.85 && table.bests(1,jx,bestp) >= 0.15 && rand(1) > 0.015 
                                   while x_normalized(ipp,jx) > 1.0d0 ||  x_normalized(ipp,jx)<0.0d0    
                                   [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;       
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                   if abs(ccc) >  1.d-8    
                                     bbb=1.1d0+(rand(1)-0.5d0)*2.0d0 ; 
                                     beta1 = alpha*3.d0*bbb*((1.d0+2.5*ff2)*rand(1) - (1.d0-ff2)*0.3d0 )  ;  
                                     x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;  
                                   else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)  ; 
                                   end
                                   end
                                 elseif table.bests(1,jx,bestp) > 0.85 || table.bests(1,jx,bestp) < 0.15 && rand(1) < 0.985 
                                   while x_normalized(ipp,jx) > 1.0d0 ||  x_normalized(ipp,jx) < 0.0d0    
                                    [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;   
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                    if abs(ccc) >  1.d-8     
                                     bbb=1.1d0+(rand(1)-0.5d0)*2.0d0 ;%*1.5d0;
                                     beta1 = alpha*3.d0*bbb*((1.d0+2.5*ff2)*rand(1) - (1.d0-ff2)*0.3d0 )  ;  
                                     x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;    
                                  else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)   ;
                                    end
                                   end
                                 else
                                   if x_normalized(ipp,jx) > 1.0d0    
                                     x_normalized(ipp,jx)=1.0d0;
                                   elseif x_normalized(ipp,jx) < 0.0d0  
                                     x_normalized(ipp,jx)=0.0d0 ;
                                   end
                                 end
                          end 
                        end
                 else
                       for jx=1:ps.D
                           x_normalized(ipp,jx)= 0.999*table.bests(1,jx,ipp)+0.001d0*table.bests(1,jx,IX(1)) ;   % x_normalized_best(ipp,:); %Local best-based parent assignment for good particles
                           if x_normalized(ipp,jx) > 1.d0 
                               x_normalized(ipp,jx)=1.d0;
                           end
                       end
                 end
            else
                     x_normalized(ipp,1:ps.D)=table.bests(1,1:ps.D,ipp); %Local best-based parent assignment during independent evaluations
            end
            
     %      end
        end %End n_par loop
    end %End while loop        

    
     %% ----------------------- Complementary functions ------------------------
     function [ffx,oox,ggx,xn_out,FEVALS] = LocalSearchMVMOSH(testcase,iii,jjj,kkk,args,xx_yy,FEsAllowed)
        global PPL GGL
        if FEsAllowed <= 0, return, end
        options=optimset('Display','off','algorithm',optimmethod,'UseParallel','never','MaxFunEvals',FEsAllowed,'FinDiffType',derivative) ;
        [Xsqp, FUN , ~ , output]=...
            fmincon(@(xx_yy)LSearch(xx_yy,testcase,iii,jjj,kkk,args),xx_yy,[],[],[],[],ps.x_min,ps.x_max,[],options);
        
        FEVALS=output.funcCount  ;
        for nvar=1:size(xx_yy,2)
            if isnan(Xsqp(1,nvar))
                Xsqp=xx_yy;
                break;
            end
        end
        
        xn_out = Xsqp;
        ffx=FUN;
        oox=PPL; 
        ggx=GGL;
    end

    function J=LSearch(xx_yy2,testcase,iii,jjj,kkk,args)
        global PPL GGL 
        [J,PPL,GGL,~] = feval(testcase,iii,jjj,kkk,args,xx_yy2);
        
    end
        
     
        function Fill_solution_archive(vvqq)
%
        no_in(ipp) = no_in(ipp)+1;
        changed = false;
        changed_best=false;

        if no_in(ipp) ==1 % the solution coming to the table for the first time large
            table.fitness(1:n_to_save,:,ipp) = 1.e200;
            table.feasibility(1:n_to_save,:,ipp) = 0;
            table.bests(1,:,ipp)   = x_normalized(ipp,:) ;
            table.fitness(1,:,ipp) = xt.fitness ;
            table.objective(1,:,ipp) = xt.objective;
            table.feasibility(1,:,ipp) = xt.feasibility    ;  %repmat(,n_to_save,1);
            no_inin(ipp)=no_inin(ipp)+1;
            changed_best=true;
            
        else % not for the first time and check for the update
           i_position =0;
               for ij=1:n_to_save 
                   if (xt.fitness < table.fitness(ij,:,ipp) && xt.feasibility == table.feasibility(ij)) || (table.feasibility(ij,:,ipp) <  xt.feasibility)                                 
                       i_position = ij;
                       changed =true;
                       if (ij<n_to_save)
                           no_inin(ipp) = no_inin(ipp)+1;  % how many times good solutions were found   
                       end
                       break;
                   end
               end
        end

        if changed   % if the new individual is better than any archived individual.
                     % Move the individuals and corresponding fitness values
                     % downward so that the individuals are sorted based on the
                     % fitness value in a descending order             
            nnnnnn = n_to_save;
            if i_position==1
              changed_best=true;
            end     
            if (no_inin(ipp) < n_to_save); nnnnnn = no_inin(ipp); end
            isdx=nnnnnn:-1:i_position+1;
            table.bests(isdx,:,ipp) = table.bests(isdx-1,:,ipp);
            table.fitness(isdx,:,ipp)= table.fitness(isdx-1,:,ipp);
            table.objective(isdx,:,ipp)= table.objective(isdx-1,:,ipp);
            table.feasibility(isdx,:,ipp)= table.feasibility(isdx-1,:,ipp);
 
            % save the new best
            table.bests(i_position,:,ipp) = x_normalized(ipp,:);
            table.fitness(i_position,:,ipp) = xt.fitness;
            table.objective(i_position,:,ipp) = xt.objective;
            table.feasibility(i_position,:,ipp) = xt.feasibility;

            % calculation of mean and variance
            if ((no_inin(ipp)>=n_to_save))
                for ivvar = 1:D
                    [meann(ipp,ivvar),shape(ipp,ivvar)] = mv_noneq(nnnnnn,table.bests(1:nnnnnn,ivvar,ipp),meann(ipp,ivvar),shape(ipp,ivvar),vvqq);
                end

            end
        end
        end
    
        function VariableSelect1(n_randomly_n)
          mode=4 ;
          n_var=ps.D;
          switch mode
            case 1
                for ii=1:n_randomly_n
                    isrepeat = false;
                    while ~isrepeat
                        inn=round(rand(1)*(n_var-1))+1;
                        if (~considered(inn))
                            isrepeat = true;
                        end
                    end
                    considered(inn)=true;
                end
            case 2
                in_randomly=0;
                isrepeat = false;
                izz(ipp)=round(rand(1)*(n_var-1))+1;  
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                        izz(ipp)=n_var;
                    end
                    considered(izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly_n)) 
                        isrepeat = true;
                    end
                end
            case 3
                in_randomly=0;
                izm(ipp)=izm(ipp)-1;
                izm(ipp)=round(rand(1)*(n_var-1))+1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                izz(ipp)=izm(ipp);
                isrepeat = false;
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                         izz(ipp)=n_var;
                    end
                    considered(izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly_n)) 
                        isrepeat = true;
                    end
                end   
            case 4
                izm(ipp)=izm(ipp)-1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                 considered(izm(ipp))=true;
                 if (n_randomly_n > 1)  
                    for ii=1:n_randomly_n-1
                        isrepeat = false;
                        while ~isrepeat
                            inn=round(rand(1)*(n_var-1))+1 ;
                            if (~considered(inn))
                                isrepeat = true;
                            end
                        end
                        considered(inn)=true;
                    end
                 end
             case 5
                  summep=sum(probab(ipp,:));
                  wahr=probab(ipp,:)/summep;
                  SS0=0.d0;
                  SS=zeros(1,n_var);
                  for imm=1:(n_var-1)
                      SS0=SS0+wahr(imm);
                      SS(imm+1)=SS0;
                  end
                  for ijlr=2:n_var
                      wahr(ijlr)=wahr(ijlr)+SS(ijlr);
                  end 
                  for ltt=1:n_randomly_n
                       isrepeat = false;
                       while  ~isrepeat
                        rnn=rand(1);
                        for irw=1:n_var
                          if considered(irw)
                           continue
                          end
                          if (irw==1)
                              unten=0.d0;
                          else
                              unten=wahr(irw-1);
                          end
                         if (rnn>=unten)&&(rnn<wahr(irw))
                              isrepeat = true;
                             considered(irw)=true;
                             break;
                          end
                        end
                       end
                  end
          end
        end
%
    function [vmean,vshape] = mv_noneq(nnnnnn,values,vmeann,vshape,vvqq)
           values_noneq(1:nnnnnn)=0.d0;
           iz =1; 
      %      vvqq
            values_noneq(iz)=values(1);
            for ii_jj=2:nnnnnn
                izz = iz;
                gleich = false;
                for kk_ii=1:izz
                    if  abs(values_noneq(kk_ii) - values(ii_jj)) <  vvqq ; 
                        gleich = true;
                        break;
                    end
                end
                if ~gleich;
                    iz = iz+1;
                    values_noneq(iz)=values(ii_jj);
                end
            end
            if (iz>1)
                   evmean=sum(values_noneq(1:iz))/iz;
                   vmean = 0.1d0*vmeann+0.9*evmean;
                   if vmean > 1.d0
                       vmean=1.d0;
                   end
                   values_noneq(1:iz)=values_noneq(1:iz)-vmean;
                   vv =sum(values_noneq(1:iz).^2)/iz;
                   if vv > 1.d-50 
                      vshape=-log(vv);
                   end
            else
              vmean=vmeann;
            end
        end 
    end

    %% Evacuated h-function
     function xnew = h_function(x_mean,s1,s2,x_curr)
%
    global mod_Method
%        
%    mod_Method=3  ; %  this is the best mapping function so far, therefore, default value
%    After local search finished, it will be changed to mod_Method=2
    %% Evacuated h-function
  %  function xnew = h_function(x_mean,s1,x_curr)
    if mod_Method==2   % this is also a new mapping function 
     if  x_curr < 0.5d0
%        %  %Forward
          if  (1.d0-x_mean) > 1.d-15     
           s11=s1/(1.d0-x_mean) ;
          else
           s11=s1/1.d-15;
          end
% 
         hmitte_f=x_mean*(1.d0-exp(-0.5d0*s11))  ; 
         hforw=x_mean*(1.d0-exp(-x_curr*s11)) ;
         hcorr=(x_mean-hmitte_f)*2d0*x_curr;
         xnew=hforw+hcorr ;
     else  
     %    %Backward
           if  x_mean > 1.d-15     
             s11=s2/x_mean ;
           else
             s11=s2/1.d-15;
           end
          
         hmitte_b=(1.0d0 - x_mean)* exp(-0.5d0*s11);
         hback=(1.0d0 - x_mean)* exp(-(1.d0-x_curr)*s11)+ x_mean;
         hcorr=hmitte_b*2d0*(1.d0-x_curr);
         xnew=hback-hcorr ;
     end
     
    elseif mod_Method==3
       if  x_curr < 0.5d0 
       %  %Forward
          if  (1.d0-x_mean) > 1.d-15     
            s11=s1/(1.d0-x_mean) ;
          else
            s11=s1/1.d-15;
          end
         hmitte_f=-x_mean/(0.5d0*s11+1.d0)+x_mean ;
         hforw=-x_mean/(x_curr*s11+1.d0)+x_mean;
         hcorr=(x_mean-hmitte_f)*2.d0*x_curr;
         xnew=hforw+hcorr ;
       else 
     %    %Backward
         if  x_mean > 1.d-15     
           s11=s2/x_mean ;
         else
           s11=s2/1.d-15;
         end
           
         hmitte_b=(1.d0-x_mean)/(0.5d0*s11+1.d0);  
         hback=(1.d0-x_mean)/((1.d0-x_curr)*s11+1.d0) + x_mean;
         hcorr=hmitte_b*2.d0*(1.d0-x_curr);
         xnew=hback-hcorr; 
       end  
    else  % this is the old mapping function
          H = x_mean .* (1.d0 - exp(-x_curr.*s1)) + ...
              (1.0d0 - x_mean) .* exp(-(1.d0-x_curr).*s2);              
          H0 = (1.d0 - x_mean) .* exp(-s2);
          H1 = x_mean .* exp(-s1);
          xnew = H + H1 .* x_curr + H0 .* (x_curr - 1.d0);
     end        
    end
    
    
    function [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)
                                              
%
    iup=round(5.d0*(1.d0-ff2))+1;
    if iup > border_gute
       iup=round((border_gute-1)*(1.d0-ff2))+1; 
    end
    ilow=1;
    bestp=round(rand(1)*(iup-ilow)+ ilow);
    worstp=-1 ;
    iup=15;
    if iup > n_par
        iup=n_par;
    end
    ilow=0;
    while (worstp <= bestp) ||  (worstp > n_par)
     worstp=round(rand(1)*(iup-ilow)+ ilow);   
     worstp=round(border_gute+(worstp-3)) ;    
    end  
    iup=worstp-1;
    ilow=bestp+1;
    onep1=round(rand(1)*(iup-ilow)+  ilow)  ;  
    onep1= IX(onep1) ;
    bestp= IX(bestp) ;
    worstp= IX(worstp) ;
   end
    
   function normalnum = normalfunc(mittel, streuung)
%
   normalnum =mittel+streuung*(sum(rand(1,12))-6.d0)/6.d0;
%
   end
%
