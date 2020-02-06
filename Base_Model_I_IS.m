% This script implements the day-to-day dynamic traffic assignment model
% developed in the reference below. It simulates the daily evolution of
% dynamic traffic flow on networks. The script simulate a scenario of link 
% disruption followed by a full recovery. The dynamic network loading 
% sub-routine is documented in 
%         Han, K., Eve, G., Friesz, T.L., 2019. Computing dynamic user 
%         equilibria on large-scale networks with software implementation. 
%         Networks and Spatial Economics, Volume 19, Issue 3, pp 869–902. 
%         Open-access URL: https://doi.org/10.1007/s11067-018-9433-y.
%
%
% Base Model I + IS assumes that the (departure time, route) choice pair 
% follow a multinomial Logit model, while the information sharing behavior
% is considered; see the reference below for more details
%
%
% INPUTS:
%         See individual parameters/variables defined below
%
%
% OUTPUTS: 
%  aggE - 
%         effective path delay aggregated (averaged) by the departure time
%         window. aggE is a 3-d matrix where the 1st dimension indicates
%         paths, the 2nd dimension indicates time window, and the 3rd 
%         dimension indicates day.
%
%  aggPath_flow - 
%         path flow within a departure window. aggPath_flow is a 3-d 
%         matrix with the same format as aggE.
% 
%  PC - 
%         perceived cost for each path (1st dimension) and time window 
%         (2nd dimension) on a given day (3rd dimension)
%
%
% DOCUMENTATION AND CITE AS:
%         Yu, Y., Han, K., Ochieng, W.Y., 2020. Day-to-Day Dynamic Traffic
%         Assignment with Imperfect Information, Bounded Rationality and 
%         Information Sharing. Transportation Research Part C, forthcoming
%


clear
clc

load 'Path_flow_data.mat'; % simulation time step (180s) and initial path departure rates
load 'OD_info.mat'; % origin-destination structure
load 'Network_planning_parameters'; % O-D demand and target arrival times (for departure time choices)


NumOD=size(OD_set,1);
time_horizon=[0, 5*3600]; % time horizon of the within-day dynamics, in seconds
T_A=T_A*3600;  % target arrival times (in second)
n_paths=size(pathDepartures,1); % number of paths
%% User-defined parameters
factor=1; % Total demand scaling factor
N=6; % memory days
lambda=0.7; % memory weight 
theta=0.004;  % Logit model dispersion parameter
Num_days=150; % total number of days for the DTD simulation
DT=900; % departure time window in seconds
TSPW=DT/dt; % Time Steps Per Window
NT=range(time_horizon)/DT; % number of departure time windows
nt=range(time_horizon)/dt; %number of time step in DNL

%% Initialize variables
pathDepartures=factor*pathDepartures; OD_demand=OD_demand*factor;
aggPath_flow=zeros(n_paths,NT,Num_days); % aggregated path departure rates
for i=1:NT
    aggPath_flow(:,i,1)=sum(pathDepartures(:,(i-1)*TSPW+1:i*TSPW,1),2)/TSPW;
end
aggE=zeros(n_paths, NT, Num_days); % aggregated travel costs for each path, departure window, and day
PC=zeros(n_paths,NT,Num_days);  % perceived travel cost for each path, departure window, and day
Flow_split=zeros(n_paths,NT,Num_days);

%% loop for T simulation days

for T=1:Num_days
    fprintf('Day no. %4.0f \n\n', T);    
    %% Dynamic Network Loading
    
    if T>50 && T<=100
        delay=DYNAMIC_NETWORK_LOADING(pathDepartures,nt,dt,'SiouxFalls6180_pp_68.mat');
    else         
        delay=DYNAMIC_NETWORK_LOADING(pathDepartures,nt,dt,'SiouxFalls6180_pp.mat');
    end
    
    %% Arrival penalty and travel cost
    time_grid=linspace(time_horizon(1),time_horizon(end),nt);
    gamma_early=0.8; gamma_late=1.8;      % coefficients of early and late arrival penalties
    AP=zeros(size(delay));             % initialize arrival penalty
    Arrival_Time=ones(n_paths,1)*time_grid+delay;
    for k=1:NumOD
        for i=1:length(ODpath_set{k,1})
            dummy=Arrival_Time(ODpath_set{k,1}(i),:)-T_A(k);
            dummy(dummy>0)=dummy(dummy>0)*gamma_late;
            dummy(dummy<=0)=-dummy(dummy<=0)*gamma_early;

            AP(ODpath_set{k,1}(i),:)=dummy;
        end
    end
    E=delay+AP;   % travel cost, 'E' stands for effective delay, meaning generalized cost

    for i=1:NT
        aggE(:,i,T)=sum(E(:,(i-1)*TSPW+1: i*TSPW),2)/TSPW;
    end
    
    %% Perceived cost and multinomial logit model
    ave_cost=zeros(1,NumOD);
    for i=1:NumOD
        index=ODpath_set{i};
        Flow_split(index,:,T)=aggPath_flow(index,:,T)/sum(sum(aggPath_flow(index,:,T)));
        ave_cost(i)=sum(sum(aggPath_flow(index,:,T).*aggE(index,:,T)*DT))/OD_demand(i);
    end
    alpha=2; % the weighting function for information sharing is g(x)=x^alpha
    if T>=N  % N: number of memory days, T: day index
        weight=Flow_split(:,:,T).^(alpha)*1;
        dummy=weight.*aggE(:,:,T);
        Sumweight=weight;
        for i=T-1:-1:T-N+1
            weight=Flow_split(:,:,i).^(alpha) *lambda^(T-i);
            dummy=dummy + weight.*aggE(:,:,i);
            Sumweight=Sumweight+weight;
        end
    else
        weight=Flow_split(:,:,T).^(alpha)*1;
        dummy=weight.*aggE(:,:,T);
        Sumweight=weight;
        for i=T-1:-1:1
            weight=Flow_split(:,:,i).^(alpha)*lambda^(T-1);
            dummy=dummy + weight.*aggE(:,:,i);
            Sumweight=Sumweight+weight;
        end
    end
    if T==1
        PC(:,:,T)=aggE(:,:,T);
    else
        PC(:,:,T)=1./Sumweight.*dummy;
    end

    
    for i=1:NumOD
        PC_alt=PC(ODpath_set{i},:,T); % matrix of perceived costs of all alternatives in OD pair i
        Den=sum(sum(exp(-theta*PC_alt)));
        for j=1:length(ODpath_set{i})
            for k=1:NT
                aggPath_flow(ODpath_set{i}(j),k,T+1)=OD_demand(i)*exp(-theta*PC_alt(j,k))/Den/DT;
            end
        end
    end
        
    %% Disaggregate path flows into smaller time steps for DNL on the next day
    for i=1:n_paths
        for j=1:NT
            pathDepartures(i , (j-1)*TSPW+1 : j*TSPW)=aggPath_flow(i,j,T+1);
        end
    end
end

