function [ corr_dim , slope_error ,r_0_new] = CorrelationDimension4( data ,r_0,PLOT_FIGURES)
%CorrelationDimension4 -- Computes an estimate of the correlation dimension
%   of the dataset "data". Uses "r_0" as the initial guess of the maximal
%   interpoint distance to consider. Also outputs the slope error from the
%   statistical test, as well as r_0_new, which is 5% larger than the
%   maximal interpoint distance used.  
% 
% PLOT_FIGURES is an option to plot a graph of "log r" vs. "log C(N,r)"
% The slope of this graph is used to estimate the correlation dim.

% The dimension of the data.
% disp(['dimension = ', num2str(length(data(1,:)))])

% Choose r_min so there are at least N^.75    distances |p_i - p_j| <r_min. 
% Choose r_max so there are no more than 50*N distances |p_i - p_j| <r_max.   
r_min = 10^-8;
r_max = 10^-1;
n_fit = 1000;

% We define R_List to be n_fit number of points logarithmically-equally
% spaced between r_min & r_max.
r_delta = (log(r_max)-log(r_min))/n_fit;
R_List = log(r_min) + r_delta * (0:n_fit) ;
R_List = exp(R_List);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(data(:,1));
model = KDTreeSearcher(data,'Distance','euclidean');

% % We want to choose r_0 so that C_RN = s*N.
ideal_multiple_s = 50.0;
appropriate_multiplier = 0;
while (appropriate_multiplier < 1)
    r0_Dist_length = 0;
    
    [~, D]= rangesearch(model,data,r_0);

    for j =1:length(D)
        r0_Dist_length = r0_Dist_length + length(D{j})-1;
    end
    appropriate_multiplier = r0_Dist_length/N;
    appropriate_multiplier = appropriate_multiplier /2; % Account for double counting.
    
    if (appropriate_multiplier<= ideal_multiple_s)
        r_0=r_0*1.5;
    end
    
    if (r_0 > 1e100)
        break
    end
end

% % % % 
%  Sort the dist list

Dist_list = [];
for j =1:length(D)
    if (length(D{j})>1)
        Dist_list = [Dist_list,D{j}(2:end)];
    end
end
Dist_list =sort(Dist_list);


% We update the r_0
if length(Dist_list) > ideal_multiple_s*N
    r_0 = Dist_list(floor(ideal_multiple_s*N));
else 
    r_0 = Dist_list(end);
end
% We slightly inflate r_0_new, so we don't risk underestimating next time,
% which is a bit expensive.
r_0_new = r_0*1.05; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_List = R_List(R_List<r_0);


counter = 1; 
% Computes C_NR
C_NR = 1:length(R_List);
for i = 1 : length(Dist_list)
    if Dist_list(i) > R_List(counter) 
        C_NR(counter) = i-1;
        counter = counter +1;
        if counter >length(R_List)
            break
        end
    end
end

if counter>length(R_List)
	counter=length(R_List);
end

% We make sure to go through all distances
if counter < length(R_List)
    R_List = R_List(1:counter);
    C_NR = C_NR(1:counter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% C_NR'
% We make sure that C_NR(r_min) > N^(.75)

while C_NR(1) < 2*N^.75
    C_NR(1)=[];
    R_List(1)=[];
    if length(C_NR) == 0 
        disp('This is a strange error')
        corr_dim =0;
        slope_error =-1;
        return
    end
end

% Here we counted all the distances twice.
C_NR = C_NR./(2* N*(N-1));


LC_NR = log(C_NR) ;
LR = log(R_List);

% Defines output
[P,S] = polyfit(LR,LC_NR,1);
corr_dim = P(1);
slope_error = sqrt((1/(length(LR)-2))*sum((LC_NR - (P(1).* LR +P(2))).^2) /  sum( (LR-mean(LR)).^2) );

if PLOT_FIGURES
    hold on
    scatter(LR,LC_NR)
    xlabel('log r')
    ylabel('log C(N,r)')
    title(['N = ',num2str(N)])
    hold off
end

end

