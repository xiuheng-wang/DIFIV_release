function [SRI_hs, P3Psi_est] = hosvdvaropt(MSI,HSI,P1,P2,Pm,ranksZ,ranksPsi,initOpt)
% =============================================================
% Computes the image fusion problem with spectral variability 
% with a coupled tensor factorization model, using a block 
% coordinate descent approach.
%
% non-intuitive paramters
% initOpt:  (integer) initialization option:
%           1 - Use the SVD of HSI and MSI matricizations,
%               combined with a coarse/smooth estimation of 
%               the variability Psi, to initialize the factrs. 
%           2 - Initialize using the results from hosvdvar.m,
%               i.e., the albegraic algorithm (this entails
%               stricter constraints on the ranks!)
%           3 - Initialize using an uniform distribution
%
% Author: Ricardo Borsoi
%
% Code related to the paper: 
%   Coupled tensor decomposition for hyperspectral and multispectral image fusion with inter-image variability
%   R.A. Borsoi, C. Prévost, K. Usevich, D. Brie, J.C.M. Bermudez, C. Richard
%   IEEE Journal of Selected Topics in Signal Processing 15 (3), 702-717, 2021.
% =============================================================


if nargin < 8
    initOpt = 1;
end

% Adjust the flag to the variability compensator, interpolation or
% pseudoinverse
if initOpt == 1
    varCompensatorFlag = 2;
elseif initOpt == 4
    initOpt = 1;
    varCompensatorFlag = 3;
end


[N1,N2,L]   = size(HSI);
[M1,M2,L_l] = size(MSI);


if initOpt == 1
    % ----> compute the initialization with the compensated variability model! <----
    % compensate variability approximately using simple normalization
    MSI_compensated = variabilityCompensator(MSI,HSI,P1,P2,Pm,varCompensatorFlag);
    P3Psi_est = MSI - MSI_compensated;

    % approximate initialization of the spatial factors with compensated variability model
    [B1z0, ~, ~] = svds(tens2mat(MSI_compensated,1,[]), ranksZ(1));
    [B2z0, ~, ~] = svds(tens2mat(MSI_compensated,2,[]), ranksZ(2));
    % Initialize spectral factor from HSI
    [B3z0, ~, ~] = svds(tens2mat(HSI,3,[]), ranksZ(3));

    % approximate initialization of Dz with HOSVD
    Dz0 = estimateCoreTensorScott(HSI,MSI_compensated,ranksZ,P1,P2,Pm,B1z0,B2z0,B3z0);

    % initialize PSI factors from the compensated residuals
    [B1psi0,   ~, ~] = svds(tens2mat(P3Psi_est,1,[]), ranksPsi(1));
    [B2psi0,   ~, ~] = svds(tens2mat(P3Psi_est,2,[]), ranksPsi(2));
    [P3B3psi0, ~, ~] = svds(tens2mat(P3Psi_est,3,[]), ranksPsi(3));
    Dpsi0 = tmprod(P3Psi_est,{B1psi0', B2psi0', P3B3psi0'},[1,2,3]);

elseif initOpt == 2
    % ----> compute the initialization with the albegraic algorithm! <----
    SRtemp = hosvdvar(MSI,HSI,ranksZ,ranksPsi,P1,P2,Pm);
    P3Psi_est = MSI - tmprod(SRtemp,Pm,3);
    [B1z0, ~, ~] = svds(tens2mat(SRtemp,1,[]), ranksZ(1));
    [B2z0, ~, ~] = svds(tens2mat(SRtemp,2,[]), ranksZ(2));
    [B3z0, ~, ~] = svds(tens2mat(SRtemp,3,[]), ranksZ(3));

    AA = kron(B3z0'*B3z0, kron(B2z0'*B2z0, B1z0'*B1z0));
    bb = tmprod(SRtemp,{B1z0', B2z0', B3z0'},[1,2,3]);
    Dz0 = reshape(AA\ bb(:), ranksZ);

    % initialize PSI factors from the compensated residuals
    [B1psi0,   ~, ~] = svds(tens2mat(P3Psi_est,1,[]), ranksPsi(1));
    [B2psi0,   ~, ~] = svds(tens2mat(P3Psi_est,2,[]), ranksPsi(2));
    [P3B3psi0, ~, ~] = svds(tens2mat(P3Psi_est,3,[]), ranksPsi(3));
    Dpsi0 = tmprod(P3Psi_est,{B1psi0', B2psi0', P3B3psi0'},[1,2,3]);
    clear SRtemp;


elseif initOpt == 3
    B1z0     = rand(M1,ranksZ(1));
    B2z0     = rand(M2,ranksZ(2));
    B3z0     = rand(L,ranksZ(3));
    Dz0      = rand(ranksZ);
    B1psi0   = rand(M1,ranksPsi(1));
    B2psi0   = rand(M2,ranksPsi(2));
    P3B3psi0 = rand(L_l,ranksPsi(3));
    Dpsi0    = rand(ranksPsi);
end





fprintf('\n\n')
niterm = 250;
rec_err = zeros(niterm,1);
for i=1:niterm
    % ---------------
    err_HSI = HSI - tmprod(Dz0,{P1*B1z0, P2*B2z0, B3z0},[1,2,3]);
    err_MSI = MSI - tmprod(Dz0,{B1z0, B2z0, Pm*B3z0},[1,2,3]) - tmprod(Dpsi0,{B1psi0, B2psi0, P3B3psi0},[1,2,3]);
    
    rec_err(i) = norm(err_HSI(:))^2+norm(err_MSI(:))^2;
    if i > 1
        if abs(rec_err(i) - rec_err(i-1))/abs(rec_err(i)) < 0.001
            fprintf('\n Finished at iteration %d... reconstr. error: %f \n', i-1, rec_err(i))
            break;
        end
    end

    fprintf('\n Iteration %d  of  %d... reconstr. error: %f', i, niterm, rec_err(i))
    
    % optimize wrt Bz1
    % optimize wrt Bz2
    % optimize wrt Bz3
    % optimize wrt Dz

    % optimize wrt Bpsi1
    % optimize wrt Bpsi2
    % optimize wrt Bpsi3
    % optimize wrt Dpsi

    [Dz0]=update_Dz(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);
    


    [B1z0,Dz0]=update_B1z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);
    [B2z0,Dz0]=update_B2z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);
    [B3z0,Dz0]=update_B3z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);
    % [B3z0,Dz0]=update_B3z_woMSI(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);



    % estimate variability
    [B1psi0,B2psi0,P3B3psi0,Dpsi0]=update_psis(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi);


end

SRI_hs = lmlragen({B1z0,B2z0,B3z0},Dz0);
P3Psi_est = lmlragen({B1psi0,B2psi0,P3B3psi0},Dpsi0);
end
% end of main function
% ..............................................................................





%*******************************************************************************  
function [B1z,Dz]=update_B1z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi)
% First remove the influence of Psi
% MSI_til = MSI - tmprod(tmprod(tmprod(Dpsi0, B1psi0,1), B2psi0,2), P3B3psi0,3);
MSI_til = MSI - tmprod(Dpsi0,{B1psi0, B2psi0, P3B3psi0},[1,2,3]);
Dz0_mat = tens2mat(Dz0,    1,[]);
MSI_1   = tens2mat(MSI_til,1,[]);
HSI_1   = tens2mat(HSI,    1,[]);

% compute the matrices related to the matricization
XX_hs = kron(B3z0, P2*B2z0) * Dz0_mat';
XX_ms = kron(Pm*B3z0, B2z0) * Dz0_mat';

% Solve generalized Sylvester equation:
% ((XX_ms'*XX_ms)\(XX_hs'*XX_hs)) * A' * (P1'*P1) + A' = 
%          ((XX_ms'*XX_ms)\(XX_ms')) * MSI_3 +
%          ((XX_ms'*XX_ms)\XX_hs') * HSI_3 * P1
% X = dlyap(A,B,C) solves the Sylvester equation A*X*B - X + C = 0
B1z = dlyap(((XX_ms'*XX_ms)\(XX_hs'*XX_hs)),...
           - (P1'*P1)',...
           + ((XX_ms'*XX_ms)\(XX_ms')) * MSI_1' + ((XX_ms'*XX_ms)\XX_hs') * HSI_1' * P1)';

% % attempted speedup
% XX_msT_XX_ms = Dz0_mat * kron(B3z0'*(Pm'*Pm)*B3z0, B2z0'*B2z0) * Dz0_mat';
% XX_hsT_XX_hs = Dz0_mat * kron(B3z0'*B3z0, B2z0'*(P2'*P2)*B2z0) * Dz0_mat';
% B1z = dlyap((XX_msT_XX_ms \ XX_hsT_XX_hs),...
%            - (P1'*P1)',...
%            + XX_msT_XX_ms \ (XX_ms' * MSI_1' + XX_hs' * HSI_1' * P1) )';

% make the columns of B orthonormal to avoid convergence issues
[B1z, RB1z] = qr(B1z,0);
if nargout==2  % Incorporate the RB1z matrices in Dz
    Dz = mat2tens(RB1z*tens2mat(Dz0,1,[]),ranksZ,1,[]);
end
end



%*******************************************************************************  
function [B2z,Dz]=update_B2z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi)
% First remove the influence of Psi
% MSI_til = MSI - tmprod(tmprod(tmprod(Dpsi0, B1psi0,1), B2psi0,2), P3B3psi0,3);
MSI_til = MSI - tmprod(Dpsi0,{B1psi0, B2psi0, P3B3psi0},[1,2,3]);
Dz0_mat = tens2mat(Dz0,    2,[]);
MSI_2   = tens2mat(MSI_til,2,[]);
HSI_2   = tens2mat(HSI,    2,[]);

% compute the matrices related to the matricization
XX_hs = kron(B3z0, P1*B1z0) * Dz0_mat';
XX_ms = kron(Pm*B3z0, B1z0) * Dz0_mat';

% Solve generalized Sylvester equation:
% ((XX_ms'*XX_ms)\(XX_hs'*XX_hs)) * A' * (P1'*P1) + A' = 
%          ((XX_ms'*XX_ms)\(XX_ms')) * MSI_3 +
%          ((XX_ms'*XX_ms)\XX_hs') * HSI_3 * P1
% X = dlyap(A,B,C) solves the Sylvester equation A*X*B - X + C = 0
B2z = dlyap(((XX_ms'*XX_ms)\(XX_hs'*XX_hs)),...
           - (P2'*P2)',...
           + ((XX_ms'*XX_ms)\(XX_ms')) * MSI_2' + ((XX_ms'*XX_ms)\XX_hs') * HSI_2' * P2)';

% make the columns of B orthonormal to avoid convergence issues
[B2z RB2z] = qr(B2z,0);
if nargout==2  % Incorporate the RB1z matrices in Dz
    Dz = mat2tens(RB2z*tens2mat(Dz0,2,[]),ranksZ,2,[]);
end
end




%*******************************************************************************  
function [B3z,Dz]=update_B3z(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi)
% First remove the influence of Psi
% MSI_til = MSI - tmprod(tmprod(tmprod(Dpsi0, B1psi0,1), B2psi0,2), P3B3psi0,3);
MSI_til = MSI - tmprod(Dpsi0,{B1psi0, B2psi0, P3B3psi0},[1,2,3]);
Dz0_mat = tens2mat(Dz0,    3,[]);
MSI_3   = tens2mat(MSI_til,3,[]);
HSI_3   = tens2mat(HSI,    3,[]);

% compute the matrices related to the matricization
XX_hs = kron(P2*B2z0, P1*B1z0) * Dz0_mat';
XX_ms = kron(B2z0, B1z0) * Dz0_mat';

% Solve generalized Sylvester equation:
% ((XX_ms'*XX_ms)\(XX_hs'*XX_hs)) * A' * (P1'*P1) + A' = 
%          ((XX_ms'*XX_ms)\(XX_ms')) * MSI_3 +
%          ((XX_ms'*XX_ms)\XX_hs') * HSI_3 * P1
% X = dlyap(A,B,C) solves the Sylvester equation A*X*B - X + C = 0
B3z = dlyap((XX_hs'*XX_hs)\(XX_ms'*XX_ms),...
           - (Pm'*Pm),...
           + ((XX_hs'*XX_hs)\(XX_ms')) * MSI_3' * Pm + ((XX_hs'*XX_hs)\XX_hs') * HSI_3')';

% make the columns of B orthonormal to avoid convergence issues
[B3z RB3z] = qr(B3z,0);
if nargout==2  % Incorporate the RB1z matrices in Dz
    Dz = mat2tens(RB3z*tens2mat(Dz0,3,[]),ranksZ,3,[]);
end
end









%*******************************************************************************  
function [Dz]=update_Dz(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi)
% First remove the influence of Psi
% MSI_til = MSI - tmprod(tmprod(tmprod(Dpsi0, B1psi0,1), B2psi0,2), P3B3psi0,3);
MSI_til = MSI - tmprod(Dpsi0,{B1psi0, B2psi0, P3B3psi0},[1,2,3]);

if (norm(eye(ranksZ(1))-B1z0'*B1z0,'fro')>1e-8) || (norm(eye(ranksZ(2))-B2z0'*B2z0,'fro')>1e-8) || (norm(eye(ranksZ(3))-B3z0'*B3z0,'fro')>1e-8) 
    [B1z0, ~, ~] = svds(tens2mat(B1z0,1,[]), ranksZ(1));
    [B2z0, ~, ~] = svds(tens2mat(B2z0,2,[]), ranksZ(2));
    [B3z0, ~, ~] = svds(tens2mat(B3z0,3,[]), ranksZ(3));
end

lam = 1;
aalpha = 0;

% Adapted from Clémence's code:
if (ranksZ(1)>size(HSI,1)||ranksZ(2)>size(HSI,2)) && ranksZ(3)>size(MSI_til,3)
    % as a last resort, try a least squares solution:
    AAA = (eye(ranksZ(3)));
    BBB = (B2z0'*(P2'*P2)*B2z0);
    CCC = (B1z0'*(P1'*P1)*B1z0);
    DDD = (B3z0'*(Pm'*Pm)*B3z0);
    EEE = (eye(ranksZ(1)*ranksZ(2)));
    AA = kron(AAA,kron(BBB, CCC)) + lam * kron(DDD, EEE);
    bb = tmprod(HSI,{B1z0'*P1', B2z0'*P2', B3z0'},[1,2,3]) + lam * tmprod(MSI_til,{B1z0', B2z0', B3z0'*Pm'},[1,2,3]);
    Dz = reshape(AA\ bb(:), ranksZ);
    Dz = (Dz);

elseif ((ranksZ(1)>size(HSI,1)||ranksZ(2)>size(HSI,2)) && ranksZ(3)<=size(MSI_til,3)) 
    A = B1z0'*(P1'*P1)*B1z0;
    B = kron(eye(ranksZ(3)), B2z0'*(P2'*P2)*B2z0);
    C = eye(ranksZ(1));
    D = kron(lam*(B3z0'*(Pm'*Pm)*B3z0), eye(ranksZ(2)));
    b_old = tmprod(HSI,{(P1*B1z0)', (P2*B2z0)', B3z0'},[1,2,3]) + lam * tmprod(MSI_til,{B1z0', B2z0', (Pm*B3z0)'},[1,2,3]);
    E = reshape(b_old, ranksZ(1), ranksZ(2)*ranksZ(3));
    if ranksZ(3) > size(MSI_til, 3)
        Dz = reshape(bartelsStewart(A,B,C,D,E),ranksZ);
    else
        Dz = reshape(bartelsStewart(C,D,A,B,E),ranksZ);
    end

else
    A = kron(B2z0'*(P2'*P2)*B2z0, B1z0'*(P1'*P1)*B1z0);
    B = lam*(B3z0'*(Pm'*Pm)*B3z0);
    if aalpha ~= 0
        A = A + aalpha*norm(A,2)^2*eye(size(A,2));
        B = B + aalpha*norm(A,2)^2*eye(size(B,2));
    end
    b_old = tmprod(HSI,{(P1*B1z0)', (P2*B2z0)', B3z0'},[1,2,3]) + lam * tmprod(MSI_til,{B1z0', B2z0', (Pm*B3z0)'},[1,2,3]);
    C = reshape(b_old, ranksZ(1)*ranksZ(2), ranksZ(3));
    %S = reshape(bartelsStewart(A,[],[],B,C), ranksZ);
    Dz = reshape(sylvester(A,B,C), ranksZ);
end
end







% =========================================================================
% =========================================================================
% =========================================================================


function [B1psi,B2psi,P3B3psi,Dpsi]=update_psis(HSI,MSI,B1z0,B2z0,B3z0,Dz0,B1psi0,B2psi0,P3B3psi0,Dpsi0,P1,P2,Pm,ranksZ,ranksPsi)
% First remove the influence of Z
P3Psi_est = MSI - tmprod(Dz0,{B1z0, B2z0, Pm*B3z0},[1,2,3]);

% initialize PSI factors from the compensated residuals
[B1psi,   ~, ~] = svds(tens2mat(P3Psi_est,1,[]), ranksPsi(1));
[B2psi,   ~, ~] = svds(tens2mat(P3Psi_est,2,[]), ranksPsi(2));
[P3B3psi, ~, ~] = svds(tens2mat(P3Psi_est,3,[]), ranksPsi(3));
Dpsi = tmprod(P3Psi_est,{B1psi', B2psi', P3B3psi'},[1,2,3]);
end

% =========================================================================
% =========================================================================
% =========================================================================

function [S]= estimateCoreTensorScott(HSI,MSI,R,P1,P2,Pm,U,V,W)
P1 = sparse(P1); P2 = sparse(P2); Pm = sparse(Pm);
U_tilde = P1*U; V_tilde = P2*V; W_tilde = Pm*W;
opts.lambda = 1;
opts.alpha  = 0;
if (R(1)>size(HSI,1)||R(2)>size(HSI,2)) && R(3)>size(MSI,3)
    fprintf('Out of the identifiability region ! \n')
%     SRI_hat = NaN; info = 'Non-identifiable';
    % approximate initialization of Dz with HOSVD
    lam = 1;
    AAA = single(eye(R(3)));
    BBB = single(V'*(P2'*P2)*V);
    CCC = single(U'*(P1'*P1)*U);
    DDD = single(W'*(Pm'*Pm)*W);
    EEE = single(eye(R(1)*R(2)));
    AA = kron(AAA,kron(BBB, CCC)) + lam * kron(DDD, EEE);
    bb = tmprod(HSI,{U'*P1', V'*P2', W'},[1,2,3]) ...
         + lam * tmprod(MSI_compensated,{U', V', W'*Pm'},[1,2,3]);
    S = reshape(AA\ bb(:), R);
    S = double(S);

elseif ((R(1)>size(HSI,1)||R(2)>size(HSI,2)) && R(3)<=size(MSI,3))
    A = U_tilde'*U_tilde;
    B = kron(eye(R(3)), V_tilde'*V_tilde);
    C = eye(R(1));
    D = kron(opts.lambda*(W_tilde'*W_tilde), eye(R(2)));
    b_old = tmprod(HSI,{U_tilde', V_tilde', W'},[1,2,3]) + opts.lambda * tmprod(MSI,{U', V', W_tilde'},[1,2,3]);
    E = reshape(b_old, R(1),R(2)*R(3));
    if R(3) > size(MSI, 3)
      S = reshape(bartelsStewart(A,B,C,D,E),R);
    else
     S = reshape(bartelsStewart(C,D,A,B,E),R);
    end 
%     SRI_hat = lmlragen({U,V,W},S);
%     info.factors = {U,V,W}; info.core = {S}; info.rank = {R};
else
    A = kron(V_tilde'*V_tilde, U_tilde'*U_tilde);
    A = A+ opts.alpha*norm(A,2)^2*eye(size(A,2));
    B = opts.lambda*(W_tilde'*W_tilde);
    B = B+ opts.alpha*norm(A,2)^2*eye(size(B,2));
    b_old = tmprod(HSI,{U_tilde', V_tilde', W'},[1,2,3]) + opts.lambda * tmprod(MSI,{U', V', W_tilde'},[1,2,3]);
    C = reshape(b_old, R(1)*R(2), R(3));
    %S = reshape(bartelsStewart(A,[],[],B,C), R);
    S = reshape(sylvester(A,B,C), R); 
%     SRI_hat = lmlragen({U,V,W},S);
%     info.factors = {U,V,W}; info.core = {S}; info.rank = {R};
end
end








% =========================================================================
function [Ym] = variabilityCompensator(Ym,Yh,P1,P2,Pm,flag_option)
[M1, M2, L_l] = size(Ym);
[N1, N2, L]   = size(Yh);

if ~exist('flag_option','var') % if flag is not selected, use spatial interpolation
    flag_option = 2;
end

if flag_option == 2 % additive scaling factors interpolation
    % Downgrade the images
    MSI_til = tmprod(tmprod(Ym,P1,1), P2, 2);
    HSI_til = tmprod(Yh,Pm,3);

    Psi_til = MSI_til - HSI_til;
    Psi_est = zeros(M1,M2,L_l);
    for ll=1:L_l
        Psi_est(:,:,ll) = imresize(Psi_til(:,:,ll), [M1 M2], 'bicubic');
    end
    % Compensate the spectral variability in Ym
    Ym = Ym - Psi_est;

elseif flag_option == 3 % plain pseudoinverse
    % Downgrade the images
    MSI_til = tmprod(tmprod(Ym,P1,1), P2, 2);
    HSI_til = tmprod(Yh,Pm,3);

    Psi_til = MSI_til - HSI_til;
    Psi_est = tmprod(tmprod(Psi_til,pinv(P1),1), pinv(P2), 2);
    % Compensate the spectral variability in Ym
    Ym = Ym - Psi_est;
end

end








