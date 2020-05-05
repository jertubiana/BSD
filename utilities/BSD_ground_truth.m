function [N,C,Palg,Pphys,O]=BSD_ground_truth(Nground_truth, Fluorescence , O , varargin)
% Performs the Blind Sparse Deconvolution Algorithm on Fluorescence data.
% Solves the optimization problem N = argmin_{N' >=0} ( 0.5 * ||F - K * N' - b||^2 + lambda
% * |N'|_1), where K is a double exponential convolution kernel and lambda
% is a sparsity coefficient. More references in [....]
% Usage: [N,C,Palg,Pphys,O] = BSD(Fluorescence,O [, P] )
% Inputs:
% - Fluorescence: A (T X Nneurons) or (Nneurons X T) matrix of DeltaF/F signals.
% - O: A structure of mandatory inputs and optional algorithm options:
%       Mandatory fields:
%         - O.Time, the number of time frames.
%         - O.nNeurons, the number of neurons/voxels.
%         - O.dt, the time interval between two measurements (= 1/f, with f the acquisition frequency)
%         Optional fields:
%         - O.adaptive (def. 1). Whether or not to perform adaptive BSD, with several (Deconvolution/Parameter refinement) steps.
%         - O.iterations (def. 5 if adaptive, else 1). The maximum number of adaptive BSD iterations.
%         - O.early_stopping (def. 1). Whether or not to stop if cost function doesn't vary anymore.
%         - O.tol_iter (def 1e-4). The relative variation threshold to determine whether the cost function varies.
%         - O.superResolution (def 1). An integer >=1. Estimate the
%         spike train at larger than sampling resolution, with time interval O.dt/O.superResolution (frequency
%         f * superResolution)
%         - O.est_a (def. 1 if adaptive, else 0). Whether to refine or not a during adaptive BSD.
%         - O.est_b (def. 1 if adaptive, else 0). Whether to refine or not b during adaptive BSD.
%         - O.est_sigma (def. 0). Whether to refine or not sigma during adaptive BSD (note: if sigma is refined, the cost function may increase during the parameter refinement step).
%         - O.est_tauRise (def. 1 if adaptive, else 0). Whether to refine or not tauRise during adaptive BSD.       
%         - O.est_tauDecay (def. 1 if adaptive, else 0). Whether to refine or not tauDecay during adaptive BSD.
%         - O.est_lambda (def. 1 if adaptive, else 0). Whether to refine or not the sparsity coefficient lambda during adaptive BSD.
%         - O.est_threshold (def. 1 if adaptive, else 0). Whether to refine or not the threshold theta during adaptive BSD.
%         - O.z1 (def 99 percentile of the normal distribution). An adimensional number to evaluate the sparsity coefficient, see article. Higher means less false positive / more false negative
%         - O.z2 (def 99 percentile of the normal distribution). An adimensional number to evaluate the sparsity coefficient, see article. Higher means less false negative / more false positive
%         - O.z3 (def 99 percentile of the normal distribution). An adimensional number to evaluate the threshold coefficient, see article. Higher means less false positive / more false positive
%         - O.u (def 0.5). An adimensional number to evaluate the threshold coefficient, see article. Lower means less false negatives/ more false positive.
%         - O.nStepsKernelInference (def 7). The maximum delay of the autocorrelation function, for initial evaluation of the convolution kernel. Increase for less variance; may lead to larger decay time than true.
%         - O.slowOptimizer (def 0) slow but more accurate deconvolution optimizer, useful sometimes for SuperResolution.
%         - O.conditioner (default: 1e-8). Conditioner for the hessian, useful sometines when using superResolution
%         - O.thresholdBeforeKernelInference (default: 0). Whether or not to threshold the spikes before the kernel optimization. Useful when tauR is of same order as Delta t or if sigma is large.
%        
%
%
% - P (optional): A structure of already known generative parameters, that serve as initializers for adaptive BSD. 
%                 For instance, one can specify a single raise & decay time for all neurons with P.tauRise = 0.1, P.tauDecay = 0.5.
%                 By default the specified parameters serve are refined with further iterations. If one wants to stick with the initial estimation, use either O.adaptive =0 or O.est_*paramName* =0, e.g. O.est_tauRise = O.est_tauDecay = 0.
%                 P can be the parameters learnt from a previous applied BSD.
%     Fields: (nNeurons X1) or (1X1, implicitely copied among neurons).
%         - P.a: transient amplitude
%         - P.b: baseline
%         - P.sigma: noise level
%         - P.tauRise: raise time
%         - P.tauDecay: decay time
%         - P.lambda: sparsity prior
%         - P.threshold: threshold
%
% Outputs:
% - N: the inferred spike train, of size (superResolution X T) * nNeurons or transpose (same format as the input data). N is not normalized by N; same scale as the original data.
% - C: the inferred convolved spike train; C = K * N + b. C can be compared with the original fluorescence signal.
% - Palg: the list of all parameters  used by algorithm. Fields include: a,b,sigma,tauRise,tauDecay,lambda,threshold. Each field is a nNeurons x 1 vector.
% - Pphy: the list of all physical parameters. Fields: tauRise,tauDecay,SNR=sigma/a,lambda,threshold.
% - O: the options structure. Contains the algorithm parameters used, and the values of the cost function.
% 
% Workflow of the algorithm:
% - Complete all empty fields of O and P to default.
% - For each neuron/voxel:
%     - (1) Perform initial estimation for all the generative parameters of the model that are not passed in P (BSD_initialization).
%     - (2) Compute sparsity coefficient lambda / threshold theta (BSD_initialization).
%     - (3) Perform sparse deconvolution with given coefficients (BSD_deconvolution).
%     - (4) If O.adaptive =1, refine generative model parameters (BSD_parameter_estimation).
%     - If O.adaptive =1, repeat steps 3&4 until convergence of maximum number of iteration is reached.
    


if length(varargin) ==1;
    P = varargin{1};
else
    P = struct;
end;



%% Check Input Size.
siz=size(Fluorescence);

if siz~=[O.Time,O.nNeurons]; Fluorescence=Fluorescence'; end; % Make Fluorescence Time X nNeurons inside the algorithm.

siz2 = size(Nground_truth);
if siz2~=[O.Time,O.nNeurons]; Nground_truth = Nground_truth'; end;
Nground_truth = double(Nground_truth);
Nground_truth = min(Nground_truth ,1);

%% Initialize the algorithm options.

if ~isfield(O,'adaptive'); O.adaptive=1; end; % By default, estimate and optimize all parameters.

if O.adaptive == 1; % If adaptive is 1 and nothing else is specified, will estimate all the parameters. If you want all params but one, O.adaptive=1 and O.est_x=0.
    
    if ~isfield(O,'iterations'); O.iterations=5; end; % Default number of iterations is 5. 
    
    if ~isfield(O,'early_stopping');   O.early_stopping  = 0; end; % whether to stop early if cost function does not vary
    
    if ~isfield(O,'tol_iter'); O.tol_iter = 1e-4; end;   % The relative variation threshold to determine whether the cost function varies.     
    
    if ~isfield(O,'est_a');     O.est_a     = 1; end;    % whether to refine a (height of spikes) or not
    
    if ~isfield(O,'est_b');     O.est_b     = 1; end;    % whether to refine b (constant baseline) or not
    
    if ~isfield(O,'est_sigma');     O.est_sigma    = 0; end;    % whether to refine sigma (noise amplitude) or not
    
    if ~isfield(O,'est_tauRise');   O.est_tauRise   = 1; end;    % whether to refine tauRise or not
    
    if ~isfield(O,'est_tauDecay');   O.est_tauDecay  = 1; end;    % whether to refine tauDecay or not

    if ~isfield(O,'est_lambda');   O.est_lambda  = 1; end;    % whether to refine the sparsity prior or not
    
    if ~isfield(O,'est_threshold');   O.est_threshold  = 1; end;    % whether to refine the threshold parameter or not       
        
else 
    if ~isfield(O,'iterations'); O.iterations=1; end; % Default number of iterations is 1.
    
    if ~isfield(O,'early_stopping');   O.early_stopping  = 0; end;
    
    if ~isfield(O,'tol_iter'); O.tol_iter = 1e-4; end;    
    
    if ~isfield(O,'est_a');     O.est_a     = 0; end;    % whether to refine a (height of spikes) or not

    if ~isfield(O,'est_b');     O.est_b     = 0; end;    % whether to refine b (offset) or not

    if ~isfield(O,'est_sigma');     O.est_sigma    = 0; end;    % whether to refine sigma (noise amplitude) or not

    if ~isfield(O,'est_tauRise');   O.est_tauRise   = 0; end;    % whether to refine tauRise or not

    if ~isfield(O,'est_tauDecay');   O.est_tauDecay  = 0; end;    % whether to refine tauDecay or not

    if ~isfield(O,'est_lambda');   O.est_lambda= 0; end;    % whether to refine the sparsity parameter or not   

    if ~isfield(O,'est_threshold');   O.est_threshold  = 0; end;    % whether to refine the threshold parameter or not   
       
end

if ~isfield(O,'superResolution');
    O.superResolution =1;
end;

if ~isfield(O,'slowOptimizer'); O.slowOptimizer = 0; end; % Set to 1 for a slow but more accurate deconvolution optimizer, useful sometimes for SR.

if ~isfield(O,'thresholdBeforeKernelInference'); O.thresholdBeforeKernelInference = 0; end; % Set to 1 to threshold the spikes before kernel inference. Useful when tauR is small or noise is large.


%% Initialize the default values of some numerical parameters.

if ~isfield(O,'z1'); O.z1=norminv(1-1e-2,0,1); end; % This number controls the trade off Deformation-Sparsity in the estimation of spikes trains.

if ~isfield(O,'z2'); O.z2=norminv(1-1e-2,0,1); end; % To compute the binarized spikes and reestimate the convolutional kernel, we threshold by default at z2*sigma/sqrt(Ak).

if ~isfield(O,'z3'); O.z3=2; end; % To compute the threshold.

if ~isfield(O,'u'); O.u=0.5; end; % To compute the threshold.


if ~isfield(O,'nStepsKernelInference');
    O.nStepsKernelInference = 7; % Number of time steps used for kernel inference.
end;
    
if ~isfield(O,'conditioner');
    if O.superResolution>1
        O.conditioner = 1e-8;
    else
        O.conditioner = 0;
    end;
end;

if ~isfield(O,'tauRiseMin');
  O.tauRiseMin = 0;
end;
if ~isfield(O,'tauDecayMin');
  O.tauDecayMin = 0;
end;
if ~isfield(O,'tauRiseMax');
  O.tauRiseMax = inf;
end;
if ~isfield(O,'tauDecayMax');
  O.tauDecayMax = inf;
end;


O.sTime = O.Time * O.superResolution;


%% Initialize the algorithm variables and parameters.

a = zeros(1,O.nNeurons);
b = zeros(1,O.nNeurons);
sigma = zeros(1,O.nNeurons);
tauRise = zeros(1,O.nNeurons);
tauDecay = zeros(1,O.nNeurons);
lambda = zeros(O.nNeurons,O.superResolution);
mu = zeros(O.nNeurons,O.superResolution);
threshold = zeros(1,O.nNeurons);
gamma = zeros(1,O.nNeurons);
delta = zeros(1,O.nNeurons);
eta = zeros(1,O.nNeurons);

N = zeros(O.sTime, O.nNeurons);
C = zeros(O.sTime, O.nNeurons);

for k = 1: O.nNeurons;
% parfor k = 1: O.nNeurons;
    [a(k),b(k),sigma(k), tauRise(k), tauDecay(k), gamma(k), delta(k), eta(k), lambda(k,:),mu(k,:), threshold(k)] = BSD_initialization(Fluorescence(:,k),O,P,k);
end;

Palg=P;
%% preallocate memory for spike train estimation.

M   =spdiags([ones(O.sTime,1) ones(O.sTime,1) ones(O.sTime,1)], -2:0,O.sTime,O.sTime); % Initialize memory for the matrix transforming calcium into spikes

I   = speye(O.sTime);                      % Identity
tmp = zeros(O.sTime,1); tmp(O.superResolution:O.superResolution:end) =1; % Position of the measurements
I1 = spdiags([tmp],0,O.sTime,O.sTime);


d0  = 1:O.sTime+1:O.sTime^2;        % index of diagonal elements of Time x Time matrices
d1  = 2:O.sTime+1:O.sTime*(O.sTime-1); % index of subdiagonal 1 elements of Time x Time matrices
d2  = 3:O.sTime+1:O.sTime*(O.sTime-2); % index of subdiagonal 2 elements of Time x Time matrices


cost = zeros(O.nNeurons,O.iterations+1,2);  % initialize cost function
iterations = zeros(O.nNeurons,1);

%% Main File: Pseudo EM Steps.

for k = 1:O.nNeurons
% parfor k = 1:O.nNeurons
    conv = (O.iterations==0); % if 0 iterations, no loop
    q = 1;
    tmp_cost = Inf(O.iterations+1,2);
    
    tmpN = zeros(O.sTime,1);
    tmpC = zeros(O.sTime,1);
    
    while conv == 0
       q = q+1;
%        tmpBigLambda = repmat(lambda(k,:)',[O.Time,1]);
       tmpBigMu = repmat(mu(k,:)',[O.Time,1]);
       tmpBigLambda = 20 * (1-Nground_truth);

       
       [tmpN,~,tmp_cost(q,1)]  = BSD_deconvolution(Fluorescence(:,k),tmpN,b(k)...
            ,gamma(k),delta(k),eta(k),tmpBigLambda,tmpBigMu,O.conditioner, O.slowOptimizer,O.superResolution, I,I1,M,d0,d1,d2);       % update inferred spike train based on estimated parameters

       [tmpC,tmp_cost(q,2),b(k),a(k),sigma(k),lambda(k,:),mu(k,:),threshold(k),gamma(k),delta(k),eta(k),tauRise(k),tauDecay(k)] ...
             = BSD_parameter_estimation(tmpN,Fluorescence(:,k),O,tauRise(k),tauDecay(k),b(k),a(k),sigma(k),threshold(k),0);    % update parameters based on previous iteration.         
       if (O.early_stopping == 0) % If cost function does not vary, and early_stoppings =0.
           if (abs((tmp_cost(q,1)-tmp_cost(q-1,1))/tmp_cost(q,1))<O.tol_iter )
               conv =1;
           end;
       end;
       if (q> O.iterations)
           conv = 1;
       end;
    end;
    iterations(k) = q-1;
    cost(k,:,:) = tmp_cost;
    N(:,k) = tmpN;
    C(:,k) = tmpC;
end;


    
%% Exit.  



Palg.a=a;
Palg.b=b;
Palg.sigma=sigma;
Palg.gamma=gamma;
Palg.delta=delta;
Palg.eta=eta;
Palg.tauRise=tauRise;
Palg.tauDecay=tauDecay;
Palg.lambda=lambda;
Palg.mu = mu * O.conditioner;
Palg.threshold=threshold;


O.posts=cost;
O.iterations = iterations;                        % record of total # of iterations

Pphys=struct;
Pphys.SNR=Palg.a./Palg.sigma;
Pphys.tauRise=Palg.tauRise;
Pphys.tauDecay=Palg.tauDecay;
Pphys.lambda=Palg.lambda;
Pphys.threshold=Palg.threshold;

O  = orderfields(O);                   % order fields alphabetically o they are easier to read
Palg      = orderfields(Palg);
Pphys=orderfields(Pphys);

if siz~=[O.sTime,O.nNeurons]; N=N'; C=C'; end; % Return Spikes and calciumTraces of same size as initial fluorescence.


end