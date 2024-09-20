% Quorum Sensing based 1-D Biofilm Model for Bacterial Populations
% with the modified explicit Tau-Leap method in "Cao et al., Efficient step size selection for the tauleaping
% simulation method, 2006"
close all
clear 
clc


tic
%% Definitions
% S1 is the downregulation state
% S2 is the upregulation state
% A = Autoinducer AHL molecules
% B = Bacteria
% E = Extracellular Polymeric Substance (EPS)

%% Parameters (Frederick et al., A mathematical model of quorum sensing regulated EPS production in biofilm communities ,2011)
N_av = 6.02214076 * 10^23; %Avogadro constant (mol^-1)
M_max = 24*10^3; % Max bacterial cell density (g m^-3)

% Bacterial growth experimental data (Fekete et al., 2010)
C_B_exp_t = [6.54E+05, 7.27E+05, 1.11E+06, 1.24E+06, 2.16E+06, 6.23E+06, 1.84E+07, 2.67E+07, 7.46E+07, 1.04E+08, 1.37E+08]; % Exp data (baffled flask) from (Fekete et al., 2010) # per ml
C_B_exp_t = C_B_exp_t.*10^3; % convert the concentration to # l^-1

% Initial values of A, B and E
init_num = 20;

% Simulation
MC = 1000; % Number of Monte Carlo loops for Stochastic Tau Leap Simulation
K = 32; % Number of compartments
% Find the domain length by using the initial concentration value of experimental bacterial concentration data
% V = h^2*L*10^3 (l) and h = L/K => V = h^2*L*10^3 = (L/K)^2*L*10^3 = N_B(0) / C_B_exp_t(0). Hence L is derived as
% L = nthroot( (K^2*init_num)/(C_B_exp_t(1)*10^3), 3); % The domain length (m)
L = 70*10^-4; % The domain length (m) (Khassekhan et al., A Nonlinear Master Equation for a Degenerate Diffusion Model of Biofilm Growth, 2009)
h = L/K; % The compartment length (m)
V = h^2*L*10^3; %Volume of the domain (l)
N_r = 10; % Number of reaction types
delta_t = 0.1; % step size (h)
t_s = 10; % Total simulation time (h)
t_vec = 0:delta_t:t_s; % Time vector (h)
N_reac = 4*K + 6*(K-1); % Number of reactions
N_spec = 3*K; % Number of species
n_c = 10; % Critical reaction threshold (Cao et al., 2006)
epsilon = 0.03; % Bounding parameter (Cao et al., 2006)

% State S1
% h_unit = 10^-4/128;
% r_alpha = 920 * 10^-6 * M_max * N_av * h^3/24; % Production rate of A at S1 (h^-1)
mu_alpha = 2.3*10^-19*N_av; % Specific production rate of A (mol per cell h^-1) (Fekete et al., 2010)
r_1 = 0.84/24; % EPS production rate at state S1 (h^-1)

% State S2
% r_beta = ((9200+920) * 10^-6 * M_max * N_av * h^3/24); % Production rate of A at S2 (h^-1)
mu_beta = (2.3*10^-18 + 2.3*10^-19)*N_av; % Specific production rate of A (mol per cell h^-1) (Fekete et al., 2010)
r_2 = 8.4/24; % EPS production rate at state S2 (h^-1)

% All states
r_sigma = 0.1109/24; % Degradation rate of A (h^-1)
% r_sigma = 0.005545; % Degradation rate of A (h^-1) (Fekete et al., 2010)
mu_g = 1/24; % Specific Growth rate of B (h^-1) (Frederick et al., 2011)
% D_A = ((0.26+0.52)/2)/24; % Diffusion coeff. of A - average value of D_A in water and D_A in fully developed biofilm (m^2 h^-1)
D_B = 6.67*10^-9/24; % Diffusion coeff. of B and E - average value, i.e., motility coeff. (m^2 h^-1)
D_A = 0; % Exclude autoinducer diffusion
% D_A = D_B;
% D_A = 10^-12/24; D_B = D_A; % (m^2 h^-1) (Khassekhan et al.,2009) 
d_A = D_A/h^2; %reaction rate for the diffusion of A (h^-1)
d_B = D_B/h^2; %reaction rate for the diffusion of B and E (h^-1)

% Scaled reaction rate constants
f_s = 10^-6; % Scaling factor
mu_alpha = f_s*mu_alpha;
mu_beta = f_s*mu_beta;
% d_A = f_s*d_A;

% gamma_QS = 6*10^-6*N_av*h^2*L; % !!! Detection threshold for quorum sensing (number of particles) (Fekete et al..2010)- multiplied by h^2*L and N_av. The original value is 70*10^-9 n mol L^-1
gamma_QS = f_s*5*10^-9*N_av; % !!! Quorum sensing detection threshold (number of particles l^-1) (Fekete et al..2010)- multiplied by h^3 and N_av. The original value is 70*10^-9 n mol L^-1



% Generation of Stoichiometric change matrix
% Columns represent species: A_1 E_1 B_1 A_2 E_2 B_2 ... A_K E_K B_K
% Rows represent reactions R_1 (K reactions) R_2 (K reactions) ... R_10 (K-1 reactions) 
nu = zeros(N_reac, N_spec); % Stoichiometric change matrix (state)
for i = 0:K-1
    nu(i+1, i*3 + 1) = 1; % Reaction 1 - Production of A
    nu(i+1+K, i*3 + 1) = -1; % Reaction 2 - Degradation of A
    nu(i+1+2*K, i*3 + 2) = 1; % Reaction 3 - Production of E
    nu(i+1+3*K, i*3 + 3) = 1; % Reaction 4 - Production of B
end

for i = 0:K-2
    nu(i+1+4*K, i*3 + 1) = -1; nu(i+1+4*K, i*3 + 4) = 1; % Reaction 5 - Forward Diffusion of A
    nu(i+1+4*K + K-1, i*3 + 2) = -1; nu(i+1+4*K + K-1, i*3 + 5) = 1; % Reaction 6 - Forward Diffusion of E
    nu(i+1+4*K + 2*(K-1), i*3 + 3) = -1; nu(i+1+4*K + 2*(K-1), i*3 + 6) = 1; % Reaction 7 - Forward Diffusion of B
    nu(i+1+4*K + 3*(K-1), i*3 + 1) = 1; nu(i+1+4*K + 3*(K-1), i*3 + 4) = -1; % Reaction 8 - Backward Diffusion of A
    nu(i+1+4*K + 4*(K-1), i*3 + 2) = 1; nu(i+1+4*K + 4*(K-1), i*3 + 5) = -1; % Reaction 9 - Backward Diffusion of E
    nu(i+1+4*K + 5*(K-1), i*3 + 3) = 1; nu(i+1+4*K + 5*(K-1), i*3 + 6) = -1; % Reaction 10 - Backward Diffusion of B
end



%% Modified Tau-Leap Algorithm (Cao et al., 2006)
C_A_t = zeros(length(t_vec),MC); %Initialize MC variables
C_B_t = zeros(length(t_vec),MC);
C_E_t = zeros(length(t_vec),MC);
C_A_x = zeros(MC,K);
C_B_x = zeros(MC,K);
C_E_x = zeros(MC,K);
state = ones(MC,length(t_vec)); % State matrix (up and downregulation)
t_detection_MC = zeros(MC,1);

% parpool('local',5) %starts the parallel pool 
% parfor i_mc = 1:MC
for i_mc = 1:MC
    %Initialize variables
    N_B = zeros(length(t_vec), K); % Number of B
    N_E = zeros(length(t_vec), K); % Number of E
    N_A = zeros(length(t_vec), K); % Number of A
    B = zeros(1, K); % Number of B - temp
    E = zeros(1, K); % Number of E - temp
    A = zeros(1, K); % Number of A - temp
    a = zeros(1, N_reac); % Propensity functions matrix
    X = zeros(1, N_spec); % State vector holding the number of molecules for each species (A_1 E_1 B_1 ... A_K E_K B_K)
    
    % Initial conditions
    init_x = K/2; %initial compartment to grow
    N_A(1,init_x) = init_num/2; % Initial number of A
    N_A(1,init_x+1) = init_num/2;
    A = N_A(1,:);
    X(1:3:end) = A;
    C_A = sum(A)./V;
    N_E(1,init_x) = init_num/2; % Initial number of E
    N_E(1,init_x+1) = init_num/2;
    E = N_E(1,:);
    X(2:3:end) = E; 
    N_B(1,init_x) = init_num/2; % Initial number of B
    N_B(1,init_x+1) = init_num/2; 
    B = N_B(1,:);
    X(3:3:end) = B; 
    
    r_g = mu_g.*B; % Growth rate of B (h^-1)
    r_alpha = mu_alpha.*B; % Production rate of A at S1 (h^-1)  (Fekete et al., 2010) ?????? Multiplication with B?
    r_beta = mu_beta.*B; % Production rate of A at S2 (h^-1)  (Fekete et al., 2010) ?????
    
    i = 1;
    t = 0;
    while t <= t_s
        clearvars ind_set_a rows cols rowsn J_ncr temp_cr J_cr mu_hat var_hat P_c
        % Decision of the states and updates of propensities depending on state
%         if sum(A) >= gamma_QS
        if C_A >= gamma_QS
            state_temp = 2; % upregulation state - S2
            a(1:K) = r_beta.*B; % Reaction 1 - Production of A
            a(2*K+1:3*K) = r_2*B; % Reaction 3 - Production of E
        else
            state_temp = 1; % downregulation state - S1
            a(1:K) = r_alpha.*B; % Reaction 1 - Production of A
            a(2*K+1:3*K) = r_1*B; % Reaction 3 - Production of E
        end
        
        % Propensity Functions independent of the state
        a(K+1: 2*K) = r_sigma*A; % Reaction 2 - Degradation of A
        a(3*K+1: 4*K) = r_g.*B; % Reaction 4 - Production of B
        a(4*K+1: 4*K+(K-1)) = d_A*A(1:K-1); % Reaction 5 - Forward Diffusion of A
        a(4*K+(K-1)+1: 4*K+2*(K-1)) = d_B*E(1:K-1);  % Reaction 6 - Forward Diffusion of E
        a(4*K+2*(K-1)+1: 4*K+3*(K-1)) = d_B*B(1:K-1); % Reaction 7 - Forward Diffusion of B
        a(4*K+3*(K-1)+1: 4*K+4*(K-1)) = d_A*A(2:K); % Reaction 8 - Backward Diffusion of A
        a(4*K+4*(K-1)+1: 4*K+5*(K-1)) = d_B*E(2:K); % Reaction 9 - Backward Diffusion of E
        a(4*K+5*(K-1)+1: 4*K+6*(K-1)) = d_B*B(2:K); % Reaction 10 - Backward Diffusion of B
    
        % Step 1 - Determine the noncritical reactions (not likely to exhaust all molecules in one step)
        Lj = zeros(1, length(a)); % maximum number of times (Lj) that Rj can fire before exhausting one of its reactants
        ind_set_a = find(a > 0); % find the indices of a > 0
        [rows, cols] = find( nu(ind_set_a,:) < 0 ); % find the indices of nu < 0  where a > 0
        rowsn(:,1) = ind_set_a(rows);
        for i_L = 1:length(rowsn)
%             Lj(1,rowsn(i_L)) = min(X(1, cols(i_L))./abs( nu( rowsn(i_L), cols(i_L) ) ) ); %Eq. 10 in (Cao et al., 2006)
            Lj(1,rowsn(i_L)) = min(floor(X(1, cols(i_L))./abs( nu( rowsn(i_L), cols(i_L) ) ) ) ); %Eq. 10 in (Cao et al., 2006)
        end
    
        J_ncr = find(Lj >= n_c); % indices of noncritical reactions
        temp_cr = find(Lj(1,ind_set_a) < n_c); % indices of critical reactions
        J_cr = ind_set_a(temp_cr);
        % Step 2
        if isempty(J_ncr)
            tau_p = Inf;
        else
            mu_hat = a(J_ncr)*nu(J_ncr,:); % Eq. 32a in (Cao et al., 2006)
            var_hat = a(J_ncr)*nu(J_ncr,:).^2; % Eq. 32b in (Cao et al., 2006)
            tau_p = min( min( max(epsilon.*X, 1)./abs(mu_hat), max(epsilon.*X, 1).^2 ./abs(var_hat) ) ); % Eq. 33 in (Cao et al., 2006)
        end
    
        % Step 3
        a0 = sum(a,'all');
        control_tau_p = 0; %while loop with control_tau_p is created to  go to step 3 from step 6 if there is a negative element in X
        while control_tau_p < 1
            if tau_p < (10/a0)
                % EXECUTE 100 TIMES SINGLE REACTION SSA (first reaction method) and return to step 1
                for i_s = 1:100
                    s_1 = rand(1, N_reac); % Generate random numbers from U(0,1)
                    tau_temp = (1./a).*log(1./s_1); %calculate the putative times of the reactions
                    [tau, I] = min(tau_temp,[],'all');  % Determine which reaction will occur
                    X = X + nu(I,:);
                    t = t + tau; % at which step the reaction occurs
                    
                    A = X(1:3:end);
                    E = X(2:3:end);
                    B = X(3:3:end);
                    % Decision of the states and updates of propensities depending on state
                    if sum(A) >= gamma_QS
                        state_temp = 2; % upregulation state - S2
                        a(1:K) = r_beta*B; % Reaction 1 - Production of A
                        a(2*K+1:3*K) = r_2*B; % Reaction 3 - Production of E
                    else
                        state_temp = 1; % downregulation state - S1
                        a(1:K) = r_alpha*B; % Reaction 1 - Production of A
                        a(2*K+1:3*K) = r_1*B; % Reaction 3 - Production of E
                    end
                    % Propensity Functions independent of the state
                    a(K+1: 2*K) = r_sigma*A; % Reaction 2 - Degradation of A
                    a(3*K+1: 4*K) = r_g.*B; % Reaction 4 - Production of B
                    a(4*K+1: 4*K+(K-1)) = d_A*A(1:K-1); % Reaction 5 - Forward Diffusion of A
                    a(4*K+(K-1)+1: 4*K+2*(K-1)) = d_B*E(1:K-1);  % Reaction 6 - Forward Diffusion of E
                    a(4*K+2*(K-1)+1: 4*K+3*(K-1)) = d_B*B(1:K-1); % Reaction 7 - Forward Diffusion of B
                    a(4*K+3*(K-1)+1: 4*K+4*(K-1)) = d_A*A(2:K); % Reaction 8 - Backward Diffusion of A
                    a(4*K+4*(K-1)+1: 4*K+5*(K-1)) = d_B*E(2:K); % Reaction 9 - Backward Diffusion of E
                    a(4*K+5*(K-1)+1: 4*K+6*(K-1)) = d_B*B(2:K); % Reaction 10 - Backward Diffusion of B
                end
                break; % go to step 1
            end
            
            % Step 4
            a0_cr = sum(a(J_cr)); % sum of all propensities of critical reactions
            tau_pp = exprnd(1/a0_cr);
        
            % Step 5
            k = zeros(1,size(nu,1)); % initialize k - # of firings of each reaction for each tau
            if tau_p < tau_pp
                tau = tau_p;            
                k(J_cr) = 0; % no critical reaction will take place 
                k(J_ncr) = poissrnd(a(J_ncr).*tau);
            else
                tau = tau_pp;
                P_c = a(J_cr)./a0_cr;
                jc = randsample(length(J_cr), 1, true, P_c); % Sample 1 point from J_cr and probabilities P
                k(J_cr) = 0; k(J_cr(jc)) = 1; % one critical reaction takes place one time
                k(J_ncr) = poissrnd(a(J_ncr).*tau);
            end
        
            % Step 6
            temp_X = X + (k*nu);
            if any(temp_X < 0)
                tau_p = tau_p/2;
                continue; % go to step 3
            else
                X = temp_X;
                t = t + tau;
                control_tau_p = 1;
            end    
      
        end
    
        % Sample the number of each species at each time step
        if t >= i*delta_t
            i = i + 1;
            N_A(i,:) = X(1:3:end);
            N_E(i,:) = X(2:3:end);
            N_B(i,:) = X(3:3:end);
            state(i_mc,i) = state_temp;            
        end
    
        A = X(1:3:end);
        E = X(2:3:end);
        B = X(3:3:end);

        C_A = sum(A)./V;
        
    end
    temp_det = find(state(i_mc,:) == 2, 1, 'first');
    if ~isempty(temp_det)
        t_detection_MC(i_mc) = t_vec( temp_det );
    else
        temp_det = NaN;
    end

    N_A_sum_t = sum(N_A,2); % sum of A wrt time
    N_B_sum_t = sum(N_B,2); % sum of B wrt time
    N_E_sum_t = sum(N_E,2); % sum of E wrt time
    
    C_A_t(:,i_mc) = N_A_sum_t./V; % concentration of A (# l^-1)
    C_B_t(:,i_mc) = N_B_sum_t./V; % concentration of B (# l^-1)
    C_E_t(:,i_mc) = N_E_sum_t./V; % concentration of E (# l^-1)

    C_A_x(i_mc,:) = N_A(end,:)./h^3; % concentration of A in each compartment (# l^-1)
    C_B_x(i_mc,:) = N_B(end,:)./h^3; % concentration of B in each compartment (# l^-1)
    C_E_x(i_mc,:) = N_E(end,:)./h^3; % concentration of E in each compartment (# l^-1)
%     toc
end

C_A_t = C_A_t*(1/f_s);
% Mean and standard deviations of the Monte Carlo results
% Time
C_A_t_mean = mean(C_A_t,2); C_A_t_std = std(C_A_t,0,2); 
C_E_t_mean = mean(C_E_t,2); C_E_t_std = std(C_E_t,0,2);
C_B_t_mean = mean(C_B_t,2); C_B_t_std = std(C_B_t,0,2);
% Space
C_A_x_mean = mean(C_A_x); C_A_x_std = std(C_A_t); 
C_E_x_mean = mean(C_E_x); C_E_x_std = std(C_E_t); 
C_B_x_mean = mean(C_B_x); C_B_x_std = std(C_B_t); 

%% Deterministic Model
% Initial values of A, B and E
init_B = 20;
init_A = 20;
init_E = 20;

% Model
model = sbiomodel('Biofilm_QS');

% Add compartments
for i_c = 1:K
    comp(1,i_c) = addcompartment(model, ['C_',num2str(i_c)]);
    comp(1,i_c).Value = 1;
    speciesB = addspecies (comp(1,i_c), 'B'); % Add Species - B
    speciesA = addspecies (comp(1,i_c), 'A'); % Add Species - A
    speciesE = addspecies (comp(1,i_c), 'E'); % Add Species - E
end

% Specify Initial Amounts of Each Species
comp(1,K/2).Species(1).InitialAmount = init_B/2; % B
comp(1,K/2+1).Species(1).InitialAmount = init_B/2; % B
comp(1,K/2).Species(2).InitialAmount = init_A/2; % A
comp(1,K/2+1).Species(2).InitialAmount = init_A/2; % A
comp(1,K/2).Species(3).InitialAmount = init_E/2; % E
comp(1,K/2+1).Species(3).InitialAmount = init_E/2; % E

% Reactions taking place in the same compartment
for j_c = 1:K
    % Reaction 1 - Production of A
    reac(1, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_', num2str(j_c), '.A + C_',num2str(j_c),'.B']); 
    kl(1, j_c) = addkineticlaw(reac(1,j_c), 'MassAction'); % Kinetic Law for reaction 1
    r_alpha(j_c) = (1/f_s)*mu_alpha.*comp(1,j_c).Species(1).Value;
    r_beta(j_c) = (1/f_s)*mu_beta.*comp(1,j_c).Species(1).Value;
    p(1, j_c) = addparameter(kl(1,j_c), ['c_',num2str(j_c)],  'Value', r_alpha(j_c)); % Rate constants 
    kl(1, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 2 - Degradation of A
    reac(2, j_c) = addreaction(model, ['C_',num2str(j_c),'.A -> null']); 
    kl(2, j_c) = addkineticlaw(reac(2, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(2, j_c) = addparameter(kl(2, j_c), ['c_',num2str(j_c)],  'Value', r_sigma); % Rate constants 
    kl(2, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 3 - Production of E
    reac(3, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_',num2str(j_c),'.E + C_',num2str(j_c),'.B']); 
    kl(3, j_c) = addkineticlaw(reac(3, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(3, j_c) = addparameter(kl(3, j_c), ['c_',num2str(j_c)],  'Value', r_1); % Rate constants 
    kl(3, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 4 - Production of B
    reac(4, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_',num2str(j_c),'.B + C_',num2str(j_c),'.B']); 
    kl(4, j_c) = addkineticlaw(reac(4, j_c), 'MassAction'); % Kinetic Law for reaction 1
    r_g = mu_g.*comp(1,j_c).Species(1).Value;
    p(4, j_c) = addparameter(kl(4, j_c), ['c_',num2str(j_c)],  'Value', r_g); % Rate constants 
    kl(4, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law
end

% Reactions taking place between subsequent compartments - Diffusion
for j_c = 1:K-1
    % Reaction 5 - Forward Diffusion of A
    reac(5, j_c) = addreaction(model, ['C_',num2str(j_c),'.A -> C_',num2str(j_c+1),'.A']); 
    kl(5, j_c) = addkineticlaw(reac(5, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(5, j_c) = addparameter(kl(5, j_c), ['c_',num2str(j_c)],  'Value', d_A); % Rate constants 
    kl(5, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 6 - Forward Diffusion of E
    reac(6, j_c) = addreaction(model, ['C_',num2str(j_c),'.E -> C_',num2str(j_c+1),'.E']); 
    kl(6, j_c) = addkineticlaw(reac(6, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(6, j_c) = addparameter(kl(6, j_c), ['c_',num2str(j_c)],  'Value', d_B); % Rate constants 
    kl(6, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 7 - Forward Diffusion of B
    reac(7, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_',num2str(j_c+1),'.B']); 
    kl(7, j_c) = addkineticlaw(reac(7, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(7, j_c) = addparameter(kl(7, j_c), ['c_',num2str(j_c)],  'Value', d_B); % Rate constants 
    kl(7, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 8 - Backward Diffusion of A
    reac(8, j_c) = addreaction(model, ['C_',num2str(j_c+1),'.A -> C_',num2str(j_c),'.A']); 
    kl(8, j_c) = addkineticlaw(reac(8, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(8, j_c) = addparameter(kl(8, j_c), ['c_',num2str(j_c)],  'Value', d_A); % Rate constants 
    kl(8, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law

    % Reaction 9 - Backward Diffusion of E
    reac(9, j_c) = addreaction(model, ['C_',num2str(j_c+1),'.E -> C_',num2str(j_c),'.E']); 
    kl(9, j_c) = addkineticlaw(reac(9, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(9, j_c) = addparameter(kl(9, j_c), ['c_',num2str(j_c)],  'Value', d_B); % Rate constants 
    kl(9, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law
    
    % Reaction 10 - Backward Diffusion of B
    reac(10, j_c) = addreaction(model, ['C_',num2str(j_c+1),'.B -> C_',num2str(j_c),'.B']); 
    kl(10, j_c) = addkineticlaw(reac(10, j_c), 'MassAction'); % Kinetic Law for reaction 1
    p(10, j_c) = addparameter(kl(10, j_c), ['c_',num2str(j_c)],  'Value', d_B); % Rate constants 
    kl(10, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law
end

%Get the Active Configuration Set for the Model.
cs = getconfigset(model,'active');
cs.SolverType = 'ode15s'; % 
cs.TimeUnits = 'hour';
solver = cs.SolverOptions;
cs.CompileOptions.DimensionalAnalysis = false;

for t_temp = 1:0.1:t_s
    cs.StopTime = t_temp;
    [t_det_f, N_det_f, names] = sbiosimulate(model);
    N_A_det_f = N_det_f(:,2:3:end); % # of A - deterministic
    N_A_det_t_f = sum(N_A_det_f,2); % # of A - deterministic - time evolution
    C_A = N_A_det_t_f(end)/V; % concentration of A (# l^-1)
    if C_A >= (1/f_s)*gamma_QS
        init2_B = N_det_f(end,1:3:end); % # of B - deterministic
        init2_A = N_det_f(end,2:3:end); % # of A - deterministic
        init2_E = N_det_f(end,3:3:end); % # of E - deterministic

        % Specify Initial Amounts of Each Species
        for ii = 1:K
            comp(1,ii).Species(1).InitialAmount = init2_B(ii); % B
            comp(1,ii).Species(2).InitialAmount = init2_A(ii); % A
            comp(1,ii).Species(3).InitialAmount = init2_E(ii); % E
        end
        % Update the state-dependent rates
        for j_c = 1:K
            % Reaction 1 - Production of A
            reac(1, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_', num2str(j_c), '.A + C_',num2str(j_c),'.B']); 
            kl(1, j_c) = addkineticlaw(reac(1,j_c), 'MassAction'); % Kinetic Law for reaction 1
            p(1, j_c) = addparameter(kl(1,j_c), ['c_',num2str(j_c)],  'Value', r_beta(j_c)); % Rate constants 
            kl(1, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law        
      
            % Reaction 3 - Production of E
            reac(3, j_c) = addreaction(model, ['C_',num2str(j_c),'.B -> C_',num2str(j_c),'.E + C_',num2str(j_c),'.B']); 
            kl(3, j_c) = addkineticlaw(reac(3, j_c), 'MassAction'); % Kinetic Law for reaction 1
            p(3, j_c) = addparameter(kl(3, j_c), ['c_',num2str(j_c)],  'Value', r_2); % Rate constants 
            kl(3, j_c).ParameterVariableNames = {['c_',num2str(j_c)]}; % Set the Kinetic Law Constants for Each Kinetic Law        
        end
        cs.StopTime = t_s - t_temp;
        [t_det, N_det, names] = sbiosimulate(model);
        break;
    end
end
t_det_c = [t_det_f(1:end-1); t_temp+t_det];

N_B_det = [N_det_f(1:end-1,1:3:end); N_det(:,1:3:end)]; % # of B - deterministic
N_A_det = [N_det_f(1:end-1,2:3:end); N_det(:,2:3:end)]; % # of A - deterministic
N_E_det = [N_det_f(1:end-1,3:3:end); N_det(:,3:3:end)]; % # of E - deterministic

N_B_det_t = sum(N_B_det,2); % # of B - deterministic - time evolution
N_A_det_t = sum(N_A_det,2); % # of A - deterministic - time evolution
N_E_det_t = sum(N_E_det,2); % # of E - deterministic - time evolution

C_A_det_t = N_A_det_t./V; % concentration of A (# l^-1)
C_B_det_t = N_B_det_t./V; % concentration of B (# l^-1)
C_E_det_t = N_E_det_t./V; % concentration of E (# l^-1)

N_B_det_x = N_B_det(end,:); % # of B - deterministic - spatial evolution
N_A_det_x = N_A_det(end,:); % # of A - deterministic - spatial evolution
N_E_det_x = N_E_det(end,:); % # of E - deterministic - spatial evolution



%% Results
t_vec = t_vec';
% Time Evoulution of the Biofilm formation in the whole domain
% All stochastic
% close all;
% h_f = figure;
% plot(t_vec, C_A_t_mean,'b-', t_vec, C_B_t_mean, 'r--', t_vec, C_E_t_mean,'k-.','LineWidth',1.25);
% xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');  xlim([0 10.25]);
% hold on; yyaxis right
% scatter(t_vec,state,10,'filled'); ylim([0 3]); ylabel('States'); 
% grid on; legend('Autoinducer (stochastic-mean)', 'Bacteria (stochastic-mean)', 'EPS (stochastic-mean)', 'States');

% Validation - Comparison with in vitro experimental data
t_exp = 0:t_s;
C_A_exp_t = [0.0, 33.4, 57.8, 11.6, 8.8, 6.0, 12.0, 19.3, 55.8, 72.8, 124.2]; % Exp autoinducer data (baffled flask) from (Fekete et al., 2010) (n mol l^-1) 
C_A_exp_t = C_A_exp_t.*10^-9.*N_av; % convert the concentration to # l^-1
C_B_exp_t = [6.54E+05, 7.27E+05, 1.11E+06, 1.24E+06, 2.16E+06, 6.23E+06, 1.84E+07, 2.67E+07, 7.46E+07, 1.04E+08, 1.37E+08]; % Exp bacteria data (baffled flask) from (Fekete et al., 2010) # per ml
C_B_exp_t = C_B_exp_t.*10^3; % convert the concentration to # l^-1

% % Autoinducer comparison - Errorbar (mean+std)
% h_f = figure;
% % h_plot = plot(t_exp, C_A_exp_t,'r*', t_det_c, C_A_det_t,'b-','LineWidth',1.25);   % concentration
% h_plot = plot(t_exp, C_A_exp_t,'r*');   % concentration
% hold on; 
% errorbar(t_vec, C_A_t_mean, 1*C_A_t_std, 1*C_A_t_std,'Color', [0.4940 0.1840 0.5560], 'Linestyle', '-','Marker', 's','MarkerSize',3,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1);
% xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
% yyaxis right
% scatter(t_vec,state,10,'filled','m'); ylim([0 3]); 
% ylabel('States'); set(gca,'ycolor','m');
% legend('Autoinducer (in vitro)', 'Autoinducer (stochastic)', 'States');
% grid on; xlim([-0.1 10.25]);

% Autoinducer comparison - Box plot
h_f = figure;
t_exp_b = 1:11;
t_vec_b = 1:delta_t:t_s+1;
t_det_c_b = t_det_c+1;
h1 = boxplot(C_A_t(1:10:end,:)',t_exp,'Whisker', Inf, 'Color', [0.4940 0.1840 0.5560]);
hold on;
h_plot = plot(t_exp_b, C_A_exp_t,'r*','LineWidth',1.25);   % concentration
h_plot2 = plot(t_det_c_b, C_A_det_t,'Color', [0.6350 0.0780 0.1840],'LineStyle', '-.', 'LineWidth',1.25);   % deterministic concentration
xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
yyaxis right
h_s = scatter(t_vec_b,state,10,'filled','m','LineWidth',1.25);
ylabel('States'); set(gca,'ycolor','m'); ylim([0 3]);
grid on;
legend([h1(1,1), h_plot(1,1), h_plot2(1,1), h_s], 'Autoinducer (stochastic)', 'Autoinducer (in vitro)', 'Autoinducer (deterministic)', 'States');


% % Bacteria comparison - Errorbar (mean+std)
% h_f = figure;
% % h_plot = semilogy(t_exp, C_B_exp_t,'r*', t_det_c, C_B_det_t,'b-','LineWidth',1.25); % concentration
% h_plot = semilogy(t_exp, C_B_exp_t,'r*','LineWidth',1.25); % concentration
% hold on;
% errorbar(t_vec, C_B_t_mean, 1*C_B_t_std, 1*C_B_t_std, 'b-','Marker', 's','MarkerSize',3,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1);
% legend('Bacteria (in vitro)', 'Bacteria (stochastic)');
% xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
% grid on; xlim([-0.1 10.25]); ylim([10^7 10^12]);

% Bacteria comparison - Box plot
h_f = figure;
h1 = boxplot(C_B_t(1:10:end,:)',t_exp,'Whisker', Inf, 'Color', 'b');
hold on;
h_plot = plot(t_exp_b, C_B_exp_t,'r*','LineWidth',1.25); % concentration
xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
grid on; ylim([10^7 10^12]);
legend([h1(1,1), h_plot(1,1)], 'Bacteria (stochastic)', 'Bacteria (in vitro)');
ax = gca;
ax.YAxis.Scale ="log";

% % EPS - Errorbar (mean+std)
% h_f = figure;
% errorbar(t_vec, C_E_t_mean, 1*C_E_t_std, 1*C_E_t_std, 'Color', [0.9290 0.6940 0.1250], 'Marker', 's','MarkerSize',3,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1); % Concentration
% % h_plot = plot(t_det_c, C_E_det_t,'b-', 'LineWidth',1.25);
% hold on;  
% xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
% yyaxis right
% scatter(t_vec,state,10,'filled','m'); ylim([0 3]); 
% ylabel('States'); set(gca,'ycolor','m');
% grid on; legend('EPS (Stochastic)', 'States');
% xlim([-0.1 10.25]);
% % plot(t_det_c, C_A_det_t,'b-', t_det_c, C_B_det_t,'r-', t_det_c, C_E_det_t,'k-','LineWidth', 1.25); %Deterministic plot concentration

% EPS Box plot
h_f = figure;
h1 = boxplot(C_E_t(1:10:end,:)',t_exp,'Whisker', Inf, 'Color', [0.8500 0.3250 0.0980]);
hold on;
xlabel('Time (h)'); ylabel('Concentration (number of molecules/l)');
yyaxis right
h_s = scatter(t_vec_b,state,10,'filled','m'); ylim([0 3]); 
ylabel('States'); set(gca,'ycolor','m');
grid on; legend([h1(1,1), h_s],'EPS (Stochastic)', 'States');
xlim([0.5 11.5]);
grid on;


% % Detection time
% h_f = figure;
% pd = fitdist(t_detection_MC,'Normal')
% histfit(t_detection_MC,10);
% title('Histogram of the detection time');
% xlabel('Time (h)'); ylabel('Number of appearance');  xlim([0 10.25]);
% grid on; 
% legend('Histogram', ['$\mathcal{N}$(', num2str(pd.mu),', ', num2str(pd.sigma^2),')'],'Interpreter','latex');

% Spatial Evolution of the Biofilm
x_plot = h:h:h*K;
h_f = figure;
C_bio_x_mean = C_B_x_mean + C_E_x_mean;
C_bio_x_mean_n = C_bio_x_mean./max(C_bio_x_mean);
C_E_x_mean_n = C_E_x_mean./ max(C_bio_x_mean);
C_B_x_mean_n = C_B_x_mean./ max(C_bio_x_mean);
plot(x_plot, C_B_x_mean_n, 'b--','LineWidth',1.25); hold on;
plot(x_plot, C_E_x_mean_n, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.','LineWidth',1.25);
plot(x_plot, C_bio_x_mean_n,  'Color', [0.1 0.6 0.3], 'LineStyle', '-', 'LineWidth',1.25);
legend('Bacteria', 'EPS' , 'Biofilm (Bacteria+EPS)');
xlabel('x-axis (m)'); ylabel('Concentration ratio');
grid on; ylim([0 1.1]);

% Concentration ratios of the biofilm
h_f = figure;
C_bio_t_mean = C_E_t_mean + C_B_t_mean;
C_bio_t_mean_n = C_bio_t_mean./ max(C_bio_t_mean);
C_E_t_mean_n = C_E_t_mean./ max(C_bio_t_mean);
C_B_t_mean_n = C_B_t_mean./ max(C_bio_t_mean);
plot(t_vec, C_B_t_mean_n, 'b--', 'LineWidth',1.25); hold on;
plot(t_vec, C_E_t_mean_n, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth',1.25);
plot(t_vec, C_bio_t_mean_n, 'Color', [0.1 0.6 0.3], 'LineStyle', '-', 'LineWidth',1.25);
xlabel('Time (h)'); ylabel('Concentration ratio');
yyaxis right
scatter(t_vec,state,10,'filled','m'); 
ylabel('States'); set(gca,'ycolor','m');
grid on; legend('Bacteria (Stochastic-mean)','EPS (Stochastic-mean)', 'Biofilm (Stochastic-mean)','States');
xlim([-0.1 10.25]); ylim([0 3]);




set(h_f,'Units','Inches');
pos = get(h_f,'Position');
set(h_f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(h_f,sprintf('Plot_Biofilm_t.pdf'),'-dpdf','-r0') %save as pdf 
% savefig(h_f,sprintf('SIR_TLW_plot_dinf_%.1f_gamma_%d_MC_%d.fig',d_inf,gamma,MC)); %save the figure file
% save('I_alpha_t.mat','I','t')
toc

save(['data_MC_', num2str(MC), '.mat']); %save all workspace variables
