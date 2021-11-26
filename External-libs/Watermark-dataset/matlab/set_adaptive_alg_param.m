%% a few possible advancements/tweaks on the code (if desired)
% this is only a few of the possible tweaks.  The adaptive algorithm can be
% tweaked to converge faster, find the dominant Doppler more quickly, have 
% a wider range of compensable Dopppler,...

crr_grd=3;          %we include a few additional carriers before carrie 0 and after carrier K-1
                    %as a guard in case Doppler shifts the frequency more than
                    %the carrier spacing.  (see "demod_PFFT.m" and "detect_PFFT.m")
                    
%% LMS.G contains all the common parameters (among the methods) used by the adaptive algorithm:
    LMS.G.fd_frwd   = 0;    %include a feed-forward compensasion in the adaptive algorithm
                            %(see "detect_PFFT.m" for details
    LMS.G.norm_grad = 0;    %enable normalizatoin of the gradeint (see "stochastic_gradient.m" for details)
    LMS.G.crr_slide = 0;    %enable carrier sliding
    LMS.G.norm_coef = 1;    %enable normalization of the amplitude of the coefficients after each iteration
                            %such normalization makes the plots showing
                            %time evolution of the coefficients more
                            %apealing.

%% LMS.P contains all the parameters used by P-FFT algorithm:
    LMS.P.mu_decdir = 0.03;  %step size of the adaptive algorithm
    LMS.P.mu_PIL    = 0.08;  %step size over pilot carriers (typically larger than mu)
    LMS.P.Thrs_Ek   = 1;     %the threshold on Ek (symol detection error) beyond which the symbol will not participate in the adaptive algorithm
    LMS.P.Thrs_g    = P;     %the threshold on g (the gradient vector) beyond which the gradient is ignored
    LMS.P.crr_slide = 0;     %crr_slide holds the value of carrier shifts throughout the program
    LMS.P.is_pilot  = 1;     %when ever is_pilot=1, the receiver uses pilots.  is_pilot will be reset after PIL carriers
    LMS.P.a         = ones(P,M);     %vector of combiner weights
    
%% LMS.S contains all the parameters used by S-FFT algorithm:
    LMS.S.mu_decdir = 0.05;  
    LMS.S.mu_PIL    = 0.1;  
    LMS.S.Thrs_Ek   = 1;
    LMS.S.Thrs_g    = S;
    LMS.S.crr_slide = 0;
    LMS.S.is_pilot  = 1;
    LMS.S.a         = ones(S,M);     %vector of combiner weights
    
%% LMS.F contains all the parameters used by F-FFT algorithm:
    LMS.F.mu_decdir = 0.04; 
    LMS.F.mu_PIL    = 0.08; 
    LMS.F.Thrs_Ek   = 1;
    LMS.F.Thrs_g    = 2;
    LMS.F.crr_slide = 0;
    LMS.F.is_pilot  = 1;
    F2=floor((F-1)/2);
    LMS.F.a         = [zeros(F2,M); ones(1,M); zeros(F2,M)]; 
%% LMS.T contains all the parameters used by T-FFT algorithm:
    LMS.T.mu_decdir = 0.005;
    LMS.T.mu_PIL    = 0.02;
    LMS.T.Thrs_Ek   = 1;
    LMS.T.Thrs_g    = 2;
    LMS.T.crr_slide = 0;
    LMS.T.is_pilot  = 1;
    LMS.T.a         = [ones(1,M); zeros(T-1,M)];
%% LMS.E contains all the parameters used for conventional equalization:
    LMS.E.mu_decdir = 0.005;
    LMS.E.mu_PIL    = 0.04; 
    LMS.E.Thrs_Ek   = 1;
    LMS.E.Thrs_g    = 2;
    LMS.E.crr_slide = 0;
    LMS.E.is_pilot  = 1;
    LMS.E.a         = 0;     %We will define the vector "a" inside detect_MFFT for conventional receiver

    

