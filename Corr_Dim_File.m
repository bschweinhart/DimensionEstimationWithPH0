function [ output ] = Corr_Dim_File(str_file_name,num_N_values)
%Corr_Dim_File 
%
% This function loads the data from the file "str_file_name" and computes
% the correlation dimension using the function "CorrelationDimension4.m". 
% See the Readme file for proper input data formating. 
% The dimension is calculated with a various subsamplings of the total
% points. 
% 
% --NOTE: Below there is an option ("PLOT_FIGURES") to graph the data used
% calculate the correlation dimension. 
%
% --INPUT
% str_file_name -- The file name of the data we are looking at
% num_N_values  -- The number of different subsamplings to compute the
%                  dimension with. Each subsample has (i/num_N_values) many 
%                  points for 1<=i<=num_N_values.
% 
% --OUTPUT
%  A matrix of size num_N_values * 4
%       output(i,1) = number of points used for i^th dimension estimate
%       output(i,2) = mean of correlation dimension estimate over trials
%       output(i,3) = standard deviation of correlation dimension estimate over trials
%       output(i,4) = the mean  std of corr. dim. resulting from statistical test 
% 
% The program will also write its output to a .txt file with the same name 
% as the data file, but with '_correlationDimension' appended to the end. 
% This output may be turned off if desired. 
FILE_OUTPUT = true;

% PLOT_FIGURES is an option to plot a graph of "log r" vs. "log C(N,r)"
% The slope of this graph is used to estimate the correlation dim.
% This will produce 'num_N_values' figures. For each number of points subsampled,
% the results from each trial are plotted on the same graph.
% If there is a systematic bias in how the data is being fit, then 
% producing these plots is a good way to find it. 
PLOT_FIGURES = false;

tic
    output = zeros(num_N_values,4);

    fileID = fopen(str_file_name,'r');  
    

    
    
    disp('Getting Data Dimensions')
    
     A = fgetl(fileID); %% Gives first Trial Number 
     A = fgetl(fileID); %% We get the first line of data 
     dim = length(str2num(A ));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Reading Files')    
    frewind(fileID);

    Master_data ={};   

    trial =0;
    while ~feof(fileID)
        
        A = fgetl(fileID);
        trial = trial+1;
        
        disp(A)
        local_data = fscanf(fileID,'%f',[dim,Inf]);
        local_data = local_data';
        Master_data{trial} = local_data;
    end
toc

tic

    disp('Computing dimension')

    N_max   = length(Master_data{1}(:,1));
    num_trials = length(Master_data);
    N_list = ceil((1:num_N_values).* N_max/num_N_values );


    r_0 = 10^-4;
    for i = 1:num_N_values
        
        if PLOT_FIGURES 
            figure;
        end
        N_local = N_list(i);
        disp(N_local)
        corr_dim  = [];
        corr_dim_error =[];
        for j =1:num_trials
            [m,e,r_0] = CorrelationDimension4(Master_data{j}(1:N_local,:) ,r_0,PLOT_FIGURES);
            corr_dim(j)         = m;
            corr_dim_error(j)   = e;
        end
            
        output(i,1) = N_local;
        output(i,2) = mean(corr_dim);
        output(i,3) = std(corr_dim);
        output(i,4) = mean(corr_dim_error);
    end
    figure 
    errorbar(output(:,1),output(:,2),output(:,3))
    if FILE_OUTPUT 
        dlmwrite([str_file_name '_correlationDimension'],output);
    end

toc
end

