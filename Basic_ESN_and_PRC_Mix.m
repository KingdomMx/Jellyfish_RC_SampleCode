function [] = Basic_ESN_and_PRC_Mix()
% Sample code to run a version of the Hybrid PRC and ESN used for the Jellyfish data analysis
% The Mux length sets the leaky-intergrator and currently set to define the future prediction length
% Created By: Max Austin
close all

%%%%% Time and Network Settings
MPLX_PRED_Len = 121;  %%% Multiplex (and Prediction) length
TWo_Pulse = 1000; %%%% Washout Time
T_Train_Frac = 1/1; %%% Fraction of data used for training: currently set to all
N = 100; %%% Number of Nodes
Spect= 0.35; %%%  Spectral Radius
Uscale = 1; %%% Baseline Input scalling, is rescaled by number of sensors and Mux len
Dt = 0.0167; %%%% Real time DT between data points


%%%% Load Your data: Read as Structs with fields for sensor and target time series
%%%% Surrogate Data for the example: Loads stimulated data from one specimen
load('SampleDataJF41.mat','TDJF41','BDJF41') 
BodyData = BDJF41;
TargetData = TDJF41;
PulseTargets = {'BODY_VX','BODY_VY','BODY_VZ'}; %%% Setting which data sets to estimate

%%% Selected Sensor Configurations and compile the PRC states and ESN input
SenseConfig =  {'Outer','R2toO2','Inner','Y2toO1'}; %%% Sensors used for the Data 
NumSensors = length(SenseConfig);
Uscale = Uscale/NumSensors; %%%% 
Ut = 2*flipud(MultiplexData(BodyData.Input,MPLX_PRED_Len))-1;
BodySenseData = zeros(MPLX_PRED_Len*NumSensors,size(Ut,2));
for LND = 1:NumSensors
    %%%% Normalize sensor data [-1,1]
    TMPSens = BodyData.(SenseConfig{LND});
    TMPSens = 2*(TMPSens-range(TMPSens)/2-min(TMPSens))/range(TMPSens);

    %%%% Multiplex Input and flip so it shift backwards in time T+n to T
    BodySenseData(MPLX_PRED_Len*(LND-1)+(1:MPLX_PRED_Len),:) = flipud(MultiplexData(TMPSens,MPLX_PRED_Len));
end

%%%%% Define Specific data lengths based on the size of the data and mux length
TTotal = length(Ut(1,:));
T = TTotal-TWo_Pulse;
T_Train = floor((T-MPLX_PRED_Len)*T_Train_Frac);

%%%%%%%%%%%%%%% Done with settings: Below runs the RC %%%%%%%%%%%%%%%%
fprintf('\n')
disp('Running Hybrid ESN Estimation')
fprintf('For Multiplex and Prediction Length: %i',MPLX_PRED_Len)
fprintf('\n')

%%%%% Set up the Multiplexed Body Data for PRC and ESN Input
 FlippedIND = fliplr(1:MPLX_PRED_Len);
CurSense = zeros(size(BodySenseData));
for LND = 1:NumSensors
    CurSense(MPLX_PRED_Len*(LND-1)+(FlippedIND(MPLX_PRED_Len):MPLX_PRED_Len),:) = BodySenseData(MPLX_PRED_Len*(LND-1)+(FlippedIND(MPLX_PRED_Len):MPLX_PRED_Len),:);
end

%%%%% Create the ESN Object and Run it
 ESN_Object = Gen_Hybrid_ESN(N,Spect,Uscale,'Activate','Tanh','InpNum',MPLX_PRED_Len*NumSensors);
InputScale =1/MPLX_PRED_Len; %%%%% Scaling the input based on Leaky Integrator size
ESN_States = ESN_Object.ESN(TTotal,CurSense*InputScale);


%%%% Check All Targets
for LND = 1:length(PulseTargets)
    %%%%% Set up Storage
    CurTarget = MultiplexData(TargetData.(PulseTargets{LND}),MPLX_PRED_Len);

    %%% Align First Term in Time for ESN and PRC with targets
    CurTarget_Align  = CurTarget(:,MPLX_PRED_Len:end)';
    ESN_States_Align  = ESN_States(1:(end-MPLX_PRED_Len+1),:);
    PRC_States_Align = CurSense(:,1:(end-MPLX_PRED_Len+1))';

    %%% Washout and Add Bias Term
    Target_Clean  = ESN_Object.Washout_Targets(TWo_Pulse,CurTarget_Align);
    ESN_States_Clean  = ESN_Object.Washout(TWo_Pulse,ESN_States_Align);
    PRC_States_Clean =  [PRC_States_Align((TWo_Pulse+1):end,:),ones([length(PRC_States_Align((TWo_Pulse+1):end,1)),1])];
    Hybrid_States_Clean = [ESN_States_Clean,PRC_States_Align((TWo_Pulse+1):end,:)];

    %%%%% Training
    W_Esn = ESN_Object.Train(T_Train,ESN_States_Clean,Target_Clean);
    W_PRC = ESN_Object.Train(T_Train,PRC_States_Clean,Target_Clean);
    W_Hybrid = ESN_Object.Train(T_Train,Hybrid_States_Clean,Target_Clean);

    %%%%% Compute Approximation
    Approx_Esn = ESN_States_Clean*W_Esn;
    Approx_PRC = PRC_States_Clean*W_PRC;
    Approx_Hyb = Hybrid_States_Clean*W_Hybrid;

    %%%%% Compute NRMSE
    R2_Esn = ESN_Object.R2_1D(Approx_Esn,Target_Clean);
    R2_PRC = ESN_Object.R2_1D( Approx_PRC,Target_Clean);
    R2_Hyb = ESN_Object.R2_1D( Approx_Hyb,Target_Clean);


    %%%%%% Saving Results data
    ESNResults.R2.(PulseTargets{LND}) = R2_Esn;
    PRCResults.R2.(PulseTargets{LND}) = R2_PRC;
    HybridResults.R2.(PulseTargets{LND}) = R2_Hyb ;

    MuxedTarget.(PulseTargets{LND}) = Target_Clean;
    ESNResults.Approx.(PulseTargets{LND}) = Approx_Esn;
    PRCResults.Approx.(PulseTargets{LND}) = Approx_PRC;
    HybridResults.Approx.(PulseTargets{LND}) =  Approx_Hyb;
end


%%%%% Plot the Results
figure
tiledlayout(1,length(PulseTargets),"TileSpacing","compact","Padding","compact")
for LND = 1:length(PulseTargets)
nexttile()
plot((0:MPLX_PRED_Len-1)*Dt,HybridResults.R2.(PulseTargets{LND}),'k','LineWidth',2)
title(PulseTargets{LND})
ylabel('Coefficient of Determination')
ylim([0,1])
xlabel('Future Prediction')
end
sgtitle('R^2 vs Prediction')


figure
tiledlayout(3,length(PulseTargets),"TileSpacing","compact","Padding","compact")
for LND = 1:length(PulseTargets)
nexttile(LND)
plot((0:size(Approx_Hyb,1)-1)*Dt,MuxedTarget.(PulseTargets{LND})(:,1),'k','LineWidth',2)
hold on
plot((0:size(Approx_Hyb,1)-1)*Dt,HybridResults.Approx.(PulseTargets{LND})(:,1),'Color',[255,20,147]/255,'LineWidth',2)
ylabel(PulseTargets{LND})
xlabel('Time')
xlim(Dt*([-1000,0]+size(Approx_Hyb,1)))
if LND == 2
    title(strcat('Predictions ',num2str(0),'s Into the future'))
end

nexttile(LND+3)
plot((0:size(Approx_Hyb,1)-1)*Dt,MuxedTarget.(PulseTargets{LND})(:,2),'k','LineWidth',2)
hold on
plot((0:size(Approx_Hyb,1)-1)*Dt,HybridResults.Approx.(PulseTargets{LND})(:,2),'Color',[255,20,147]/255,'LineWidth',2)
title(PulseTargets{LND})
ylabel(PulseTargets{LND})
xlabel('Time')
xlim(Dt*([-1000,0]+size(Approx_Hyb,1)))
if LND == 2
    title(strcat('Predictions ',num2str(1*Dt),'s Into the future'))
end

nexttile(LND+6)
plot((0:size(Approx_Hyb,1)-1)*Dt,MuxedTarget.(PulseTargets{LND})(:,end),'k','LineWidth',2)
hold on
plot((0:size(Approx_Hyb,1)-1)*Dt,HybridResults.Approx.(PulseTargets{LND})(:,end),'Color',[255,20,147]/255,'LineWidth',2)
title(PulseTargets{LND})
ylabel(PulseTargets{LND})
xlabel('Time')
xlim(Dt*([-1000,0]+size(Approx_Hyb,1)))
if LND == 2
    title(strcat('Predictions ',num2str((MPLX_PRED_Len-1)*Dt),'s Into the future'))
end
end
legend({'Target','Prediction'})
sgtitle('Prediction Over Time')


assignin("caller","MuxedTarget",MuxedTarget);
assignin("caller","ESNResults",ESNResults);
assignin("caller","PRCResults",PRCResults);
assignin("caller","HybridResults",HybridResults);

end

function [MultData,varargout] = MultiplexData(Data,CycleLength)
%%%%% Store data that is shifted by time and the length
MultData = zeros([CycleLength,length(Data)-CycleLength]);
    for ind = 1:CycleLength
        %%% Cut off end
        MultData(ind,:) =  Data(ind:end-(CycleLength-(ind-1)));
    end

    if nargout >= 2
        TMP = MultData;
        TMP(1,:) = 0; 
        varargout{1}= cumsum(TMP,1);
    end

end

%%%%%%%%%%%%%% GENERAL ESN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
function [ESN_OBJ] = Gen_Hybrid_ESN(N,Spect,Input_Scale,varargin)
%Generates an ESN object
% N is number of nodes, Spect ois spectral raduis, 
% Input_scale is the scale of the builtin random input (only used if you create the noise here)
% ActivateFn is the activation function selected
% Created By: Max Austin

%%%%%% Variable Input Arguments
if  any(strcmp(varargin,'RSeed')) %%% Random seed internal weights
    index = find(strcmp(varargin,'RSeed'));
    RFun_Seed = varargin{index+1} ;
else
    RFun_Seed = 1 ;
end

if  any(strcmp(varargin,'InpRSeed')) %%% Random seed input weights
    index = find(strcmp(varargin,'InpRSeed'));
    RFun_Seed2 = varargin{index+1} ;
else
    RFun_Seed2 = RFun_Seed+1 ;
end

if  any(strcmp(varargin,'InpNum')) %%% Random seed input weights
    index = find(strcmp(varargin,'InpNum'));
    N_Inputs = varargin{index+1} ;
else
      N_Inputs = 1;
end

if any(strcmp(varargin,'Activate')) %%% activation function
    index = find(strcmp(varargin,'Activate'));
    ActivateFn = varargin{index+1} ;
else
    ActivateFn = 'Tanh';
end

if any(strcmp(varargin,'Inbias')) %%% Input bias
    index = find(strcmp(varargin,'Inbias'));
    Offset = varargin{index+1} ;
else
    Offset = 0;
end

if any(strcmp(varargin,'Afun_Param'))  %%% Parameters for the activation function
    index = find(strcmp(varargin,'Afun_Param'));
    alpha = varargin{index+1} ;
else
   alpha = 1;
end

if any(strcmp(varargin,'Init'))  %%% Parameters for the inital condition
    index = find(strcmp(varargin,'Init'));
    ESN_OBJ.Init = varargin{index+1} ;
else
   ESN_OBJ.Init = ones([1,N]);
end

%%%%%%% Parameters related to external nodes %%%%%%%%%
if any(strcmp(varargin,'External_Size'))  %%% Parameter with number of external nodes
    index = find(strcmp(varargin,'External_Size'));
    ESN_OBJ.NEx = varargin{index+1};
    ESN_OBJ.SpectEx = Spect; 
else
   ESN_OBJ.NEx = 0;
   ESN_OBJ.SpectEx = 0;
end

if any(strcmp(varargin,'External_Spect'))  %%%  Parameter with external Spectral density
    index = find(strcmp(varargin,'External_Spect'));
    ESN_OBJ.SpectEx = varargin{index+1};
end

if  any(strcmp(varargin,'ExtRSeed')) %%% Random seed input weights
    index = find(strcmp(varargin,'ExtRSeed'));
    RFun_Seed3 = varargin{index+1} ;
else
    RFun_Seed3 = RFun_Seed2 + 1;
end


%%%%%%%%%   ESN Object  %%%%%%%%%%%

% Save number of Nodes
ESN_OBJ.N = N;
ESN_OBJ.AF = ActivateFn;

% Random value generation
ESN_OBJ.Rnd1D = @(Low,Upp,N) (Low+(Upp-Low)*rand([N,1]));
ESN_OBJ.Rnd2D = @(Low,Upp,N) (Low+(Upp-Low)*rand([N(1),N(2)]));

% Random input weights
rng(RFun_Seed2); % creates a reproducible random set
% ESN_OBJ.InputW = ESN_OBJ.Rnd1D(-Input_Scale,Input_Scale,N); % input weights
ESN_OBJ.InputW = ESN_OBJ.Rnd2D(-Input_Scale,Input_Scale,[N,N_Inputs]); % input weights %%%% For upt to N inputs
ESN_OBJ.InputBias = Offset;

% Set the spectral raduis of the internal weights
rng(RFun_Seed); % creates a reproducible random set
TMP_IW = ESN_OBJ.Rnd2D(-1,1,[N,N]).*ESN_OBJ.Rnd2D(0,1,[N,N]);  % Internal weights
Eval = eig(TMP_IW);
ESN_OBJ.Internal = Spect*TMP_IW /max(abs(Eval));  % Internal weights

% Set the spectral raduis of the Extenal Nodes
if ESN_OBJ.NEx >0 &&  ESN_OBJ.NEx ==N
    rng(RFun_Seed3);
    TMP_IW = ESN_OBJ.Rnd2D(-1,1,[N,ESN_OBJ.NEx]).*ESN_OBJ.Rnd2D(0,1,[ESN_OBJ.NEx,N]);  % Internal weights
    Eval = eig(TMP_IW);
    ESN_OBJ.External = ESN_OBJ.SpectEx*TMP_IW/max(abs(Eval));  % Internal weights
elseif ESN_OBJ.NEx >0
    rng(RFun_Seed3);
    ESN_OBJ.External = ESN_OBJ.SpectEx*ESN_OBJ.Rnd2D(0,1,[N,ESN_OBJ.NEx]);  % Internal weights
else
    ESN_OBJ.External = 0;  % Internal weights
end


ESN_OBJ.alpha  = alpha; 

%%% Run the ESN
if ESN_OBJ.NEx >0
    ESN_OBJ.ESN = @(Time,u,X_External)  (ForFN_ESN(Time,u,X_External,ESN_OBJ));
else
    ESN_OBJ.ESN = @(Time,u)  (ForFN_ESN(Time,u,0,ESN_OBJ));
end

%%% Washout the ESN and set up states for estimate
ESN_OBJ.Washout = @(TWo,X) ( [X((TWo+1):end,:),ones([length(X((TWo+1):end,1)),1])]);
ESN_OBJ.Washout_Targets = @(TWo,Y) (Y((TWo+1):end,:));


%%%% Linear Regression
ESN_OBJ.Train = @(T_Train,X,Y) (X(1:T_Train,:)\Y(1:T_Train,:));
ESN_OBJ.Train_Ridge = @(T_Train,X,Y,ESN_OBJ) ComputRidge(T_Train,X,Y,ESN_OBJ.N+ESN_OBJ.NEx);

ESN_OBJ.NRMSE = @(Y_hat,Y_Eval) (sqrt(mean((Y_hat-Y_Eval).^2))./std(Y_Eval));
ESN_OBJ.R2 = @(Estimate,Real) (1-((sum((Real-Estimate).^2))/(sum((Real - mean(Real)).^2))));
ESN_OBJ.R2_1D = @(Estimate,Real) (1-((sum((Estimate-Real).^2,1))./(sum((Real - mean(Real)).^2,1))));

end

function [Weights] = ComputRidge(T_Train,X,Y,N)
% Runs the ridge regression algorithm based on the 

% X_Train = X(1:T_Train,1:end-1); %%% training without Bias
X_Train = X(1:T_Train,:); %%% training without Bias
Y_Train = Y(1:T_Train,:);
lambda = eig(X_Train'*X_Train); %%% Eigenvalue
Sigma = svd(X_Train); %%% singular value decomposition


%%% Find beta
beta = 1e-13; %% initial conditions
Betas = zeros([1,N+1]);
for df = (N+1):-1:1 %% select df to be whole numbers 
    dF_dbeta = sum(lambda./((lambda+beta).^2));
    F = df-sum(lambda./(lambda+beta));
    Betas(df) = beta -  F/dF_dbeta;  %% Newtons methods 
    beta = Betas(df);
end


%%%  Compute the ridge regression
Weights = (X_Train'*X_Train+Betas*eye(N+1))\(X_Train'*Y_Train); % may need to remove bias though


end

function [X] = ForFN_ESN(Time,u,X_External,ESN_OBJ)
%Generate ESN STate For loop

X = ones([Time,ESN_OBJ.N]);

X(1,:) = ESN_OBJ.Init;

if size(u,1)<Time
    u = u';
end

%%%%% Adjusting size incase ther are no external Nodes
if all(size(X_External) == [1,1])
    X_External = zeros([Time,1]);
else
    X_External = X_External./(max(abs(X_External),[],1));
end

for inc = 2:Time
    Activate_Input = ((ESN_OBJ.Internal*X(inc-1,:)') + (ESN_OBJ.External*X_External(inc-1,:)') + (ESN_OBJ.InputW*u(inc-1,:)' + ESN_OBJ.InputBias))';
    
    %%%% Switches the activation function defaults to Tanh
    if strcmp(ESN_OBJ.AF,'Sigmoid')
        X(inc,:)= 1./(1+exp(-1*Activate_Input));
    elseif strcmp(ESN_OBJ.AF,'Atan')
        X(inc,:)= atan(Activate_Input);
    elseif strcmp(ESN_OBJ.AF,'Ln')
        X(inc,:)= log(1+Activate_Input);
   elseif strcmp(ESN_OBJ.AF,'ELU')
       if  Activate_Input < 0
           X(inc,:)= ESN_OBJ.alpha*exp(Activate_Input-1);
       else
           X(inc,:)= Activate_Input;
       end
    elseif strcmp(ESN_OBJ.AF,'Linear')
        X(inc,:)= Activate_Input;
    elseif strcmp(ESN_OBJ.AF,'Logistic')
        X(inc,:)= ESN_OBJ.alpha*Activate_Input.*(1-Activate_Input);
    else 
        X(inc,:)= tanh(Activate_Input);
    end
end

end


