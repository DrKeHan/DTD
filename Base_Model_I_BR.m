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
% Base Model I + BR assumes that the (departure time, route) choice pair 
% follow a multinomial Logit model, while the boundedly rational behavior
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
delta=400; % bounded rationality 'indifference band' in second
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
 
    if T>=N  % N: number of memory days, T: day index
        dummy=aggE(:,:, T);
        for i=T-1:-1:T-N+1
            dummy=dummy + lambda^(T-i)*aggE(:,:,i);
        end
        PC(:,:,T)=(1-lambda)/(1-lambda^N)*dummy;
    else
        dummy=aggE(:,:, T);
        for i=T-1:-1:1
            dummy=dummy + lambda^(T-i)*aggE(:,:,i);
        end
        PC(:,:,T)=(1-lambda)/(1-lambda^T)*dummy;
    end
    

    for i=1:NumOD
        PC_alt=PC(ODpath_set{i},:,T); % matrix of perceived costs of all alternatives in OD pair i
        Den=sum(sum(exp(-theta*PC_alt))); 
        for j=1:length(ODpath_set{i})
            for k=1:NT
                FromOthers=0;
                for jj=1:length(ODpath_set{i})
                    for kk=1:NT
                        if jj~=j || kk~=k
                            FromOthers=FromOthers+...
                                aggPath_flow(ODpath_set{i}(jj),kk,T)*exp(-theta*PC_alt(j,k))/(Den-exp(-theta*PC_alt(jj,kk))+exp(-theta*(PC_alt(jj,kk)-delta))); % switched from other alternatives
                        end
                    end
                end
                aggPath_flow(ODpath_set{i}(j),k,T+1)=...
                    aggPath_flow(ODpath_set{i}(j),k,T)*exp(-theta*(PC_alt(j,k)-delta))/(Den-exp(-theta*PC_alt(j,k))+exp(-theta*(PC_alt(j,k)-delta)))... % Not switching
                    +FromOthers;
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

