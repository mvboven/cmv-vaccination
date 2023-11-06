/* Inference for CMV cross-sectional serological data using transmission model, assuming       */
/* endemic equilibrium and short disease approximation. Initial code has been drafted by Chris */
/* van Dorp (CvD)(no sex), and subsequently adapted by Sophia de Jong and Michiel van Boven    */
/* (MvB). Stan/HMC yields identical estimates as earlier Mathematica program based on MH       */
/* When replacing the transmission model with splines for the prevalences, the model yields    */
/* identical results to the earlier JAGS model; cf van Boven et al (2017) PCB (mvb17)          */
/* 12/2018: added code/lik contributions for data on congenital infections, ensuring that      */
/* fraction that is infected congenitally is approximately 0.5% - see PhD thesis M Korndewal   */ 
/* Copyright by CvD and MvB and licensed under BSD (3-clause) (2017-2019)                      */

functions {
  // construct b-spline basis functions; taken from Milad Kharratzadeh's (MK) example 
  // see https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order); 
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piecewise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  int<lower=0> N;                                     // number of subjects
  int<lower=1> DeltaA;                                // length of age intervals
  int<lower=1> A;                                     // number of age intervals
  matrix<lower=0>[A, A] Contact_MM;                   // gender- and age-specific contact matrices
  matrix<lower=0>[A, A] Contact_FM;                   // male to female contact matrix
  matrix<lower=0>[A, A] Contact_MF;                   // female to male contact matirx
  matrix<lower=0>[A, A] Contact_FF;                   // female to female
  int<lower=0, upper=DeltaA*A> Ages[N];               // subject ages (rounded from months (LFTB2 in P2) to years)
  real Titers[N];                                     // antibody titers
  int Censor[N];                                      // 0 = normal, 1 = censored, 2 = spike (see mvb17)
  real RightCensor;                                   // titers above this value are right-censored
  real MuS;                                           // mean of classification mixture (S) (estimated in mvb17)
  real MuL;                                           // mean of classification mixture (L)
  real MuB;                                           // mean of classification mixture (B)
  real<lower=0> SigmaS;                               // standard deviations of the classification mixture; uninfected component
  real<lower=0> SigmaL;                               // standard deviations of the classification mixture; unfected component
  real<lower=0> SigmaB;                               // standard deviations of the classification mixture; infected with raised antibodies
  int<lower=1> numbertestedinfants;                   // number of infants tested for congenital CMV (cCMV) - 31484; Epi Inf (2016) 7, 1520-1527
  int<lower=1> numbercCMVinfants;                     // number of infants positive for cCMV - 154 (approx 0.5%); Epi Inf (2016) 7, 1520-1527
  real<lower=0> Penalty;                              // accuracy for estimation of the fois/S0 using LHS ~ N(RHS,1/Penalty); future: use solver
  int<lower=0, upper=1> Gender[N];                    // 0 = female, 1 = male
  vector[A] BirthContribution;                        // prob distribution of ages of mothers with newborns in 2006 (from http://www.cbs.nl) 
  int num_knots;                                      // number of spline knots
  vector[num_knots] knots;                            // the sequence of knots 
  int spline_degree;                                  // the spline degree (order - 1)	
  real ts[DeltaA*A];                                  // ages at which splines are calculated
  real<lower=0, upper=1> reducinf;                    // infectivity reduction in L compared to B
  int<lower=0, upper=1> mode;                         // 0 = regular sampling, 1 = sampling to compute WBIC
}

transformed data {
  real<lower=0, upper=1> watanabe_beta;               // WBIC indicator
    
  /* auxiliary variables for vectorisation */
  vector[1] Zero;
  vector[A] LongZeros;
  vector[DeltaA*A] LongOnes;
  vector[DeltaA*A+1] LlongOnes;
    
  /* indices for transforming lambda into longLambda and vice versa */
  int<lower=1, upper=A> ExtendIdxs[DeltaA*A];
  int<lower=1, upper=DeltaA*A+1> ReduceIdxs[A];
  int<lower=1, upper=DeltaA*A+1> ReduceIdxsRightShift[A];

  /* splines, after MK */
  int num_basis = num_knots + spline_degree - 1;        // number of B-splines
  matrix[num_basis, DeltaA*A] B;                        // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2 * spline_degree + num_knots] ext_knots;      // extended knots
  
  /* spline construction; after MK */	
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis){
    B[ind,:] = to_row_vector(build_b_spline(ts, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, DeltaA*A] = 1;                                      
  
  /* Sampling temperature; see http://www.jmlr.org/papers/volume14/watanabe13a/watanabe13a.pdf for theory   */
  /* and for application in Stan http://tbz533.blogspot.com/2017/04/how-to-calculate-wbic-with-stan.html    */
  if ( mode == 0 ) { // normal sampling
      watanabe_beta = 1.0;
  }
  else { // mode == 1, WBIC sampling
      watanabe_beta = 1.0/log(N); // needs checking that N is indeed the number of observations
  }
	
  /* preliminaries for vectorisation and to improve readibility */
  Zero = rep_vector(0.0, 1);                          // a vector-version of 0 (for cumulative sums starting at 0)
  LongZeros = rep_vector(0.0, A);
  LongOnes = rep_vector(1.0, DeltaA*A);               // DeltaA*A ones (same length as longLambda and longPi)
  LlongOnes = rep_vector(1.0, DeltaA*A+1);            // DeltaA*A+1 ones (same length as S, L, and B)
  for ( aa in 1:DeltaA*A ) {                          // aa-1 = 0,1,2,3,4,5,6,7,8,9,10,...
      int j;
      j = (aa-1)/DeltaA + 1;                          // j = 1,1,1,1,1, 2,2,2,2,2, 3,...
      ExtendIdxs[aa] = j;
  }
  for ( j in 1:A ) {                                  // j = 1,2,3,...
      ReduceIdxs[j] = 1 + (j-1)*DeltaA;               // x-1 = 0,5,10,...,DeltaA*(A-1)
      ReduceIdxsRightShift[j] = 1 + j*DeltaA;         // x-1 = 5,10,15,...,DeltaA*A
  }
}

parameters {
  real<lower=0> beta1;                                // infectivity after primary infection
  real<lower=0> beta2;                                // infectivity after reactivation (and/or re-infection) in L or L/B (model-dependent)
  real<lower=0, upper=1> z;                           // reduction in susceptibility to reinfection
  real<lower=0, upper=1> probLtoB;                    // probability that reactivation or reinfection in L leads to antibody boosting (transfer to B)
  real<lower=0> S0;                                   // fraction of the population not vertically infected
  real<lower=0> nu;                                   // probability of vertical transmission
  real<lower=0, upper=1> qcCMV;                       // prob of congenital infection during acute infection of mother (possibly times infectious period)
  matrix<lower=0>[2, num_basis] a_raw;                // spline basis functions; 1 = female, 2 = male ; NB check whether <lower=0> is needed - depends on model
  vector<lower=0>[A] lambda_f;                        // forces of infection on the age-intervals in females
  vector<lower=0>[A] lambda_m;                        // forces of infection on the age-intervals in males
}

transformed parameters {
  /* reactivation rates */
  vector<lower=0>[DeltaA*A] rho_f;                    // reactivation rate in females on intervals            
  vector<lower=0>[DeltaA*A] rho_m;                    // reactivation rate in females on intervals
  
  /* prevalences */
  vector<lower=0, upper=1>[DeltaA*A+1] S_f;           // susceptible prevalence at points of intervals (female)
  vector<lower=0, upper=1>[DeltaA*A+1] S_m;           // susceptible prevalence at points of intervals (male)
  vector<lower=0, upper=1>[DeltaA*A+1] L_f;           // latently infected (female)
  vector<lower=0, upper=1>[DeltaA*A+1] L_m;           // latently infected (male)
  vector<lower=0, upper=1>[DeltaA*A+1] B_f;           // infected with boosted titers (female)
  vector<lower=0, upper=1>[DeltaA*A+1] B_m;           // infected with boosted titers (male)
  	
  /* auxiliary vectors to efficiently solve ODEs for S, L, and B */
  vector<lower=0, upper=1>[DeltaA*A+1] X_f;           // latently infected at birth (female)
  vector<lower=0, upper=1>[DeltaA*A+1] X_m;           // latently infected at birth (male)
  vector<lower=0>[DeltaA*A+1] Y_f;                    // =(L_f-X_f)/X_f (ratio of L's infected horizontally and vertically)
  vector<lower=0>[DeltaA*A+1] Y_m;                    // =(L_m-X_m)/X_m
	
  /* lambda hat (i.e. the force of infection) should be very similar to lambda */
  vector<lower=0>[A] lambda_hat_f;
  vector<lower=0>[A] lambda_hat_m;
  
  /*  perinatal and congenital infection */
  vector<lower = 0>[A] S_birth;
  real<lower = 0> frac_S_birth;
  real<lower = 0> pcCMV;
  
  /* long versions of lambda, pi, and others */ 
  vector<lower=0>[DeltaA*A] longLambda_f;
  vector<lower=0>[DeltaA*A] longLambda_m;
  vector<lower=0>[DeltaA*A] longPi_f;
  vector<lower=0>[DeltaA*A] longPi_m;

  /* cumulated prevalences */ 
  vector<lower=0>[A] aggr_S_f;
  vector<lower=0>[A] aggr_L_f;
  vector<lower=0>[A] aggr_L_m;
  vector<lower=0>[A] aggr_B_f;
  vector<lower=0>[A] aggr_B_m;
    
  /* (cubic) B-splines for reactivation rates                                                                        */
  /* P-splines but computationally too costly here */
  rho_f = to_vector(a_raw[1]*B);                      
  rho_m = to_vector(a_raw[2]*B);                      
   
  /* make long versions of lambda and pi; rho is already long (dim: DeltaA*A)    */
  /* for notational convenience, define pi = rho + z * lambda                    */
  longLambda_f = lambda_f[ExtendIdxs];
  longLambda_m = lambda_m[ExtendIdxs];
  longPi_f = rho_f + z * longLambda_f;                
  longPi_m = rho_m + z * longLambda_m;
    
  /* solution of the ODEs S, L, B, and intermediates X and Y in terms of the foi (lambda) in females and males */
  /* X : perinatally infected and still in L                                                                   */
  /* Y : ratio of persons in L that have been infected after birth (L-X) over those infected perinatally (X)   */
  S_f = S0 * exp(-cumulative_sum(append_row(Zero, longLambda_f)));
  S_m = S0 * exp(-cumulative_sum(append_row(Zero, longLambda_m)));

  X_f = (1.0 - S0) * exp(-cumulative_sum(append_row(Zero, probLtoB * longPi_f)));
  X_m = (1.0 - S0) * exp(-cumulative_sum(append_row(Zero, probLtoB * longPi_m)));
  
  Y_f = cumulative_sum(append_row(Zero, longLambda_f .* (S_f[:DeltaA*A] ./ X_f[:DeltaA*A]) .* (LongOnes - exp(-(longLambda_f - probLtoB * longPi_f))) ./ (longLambda_f - probLtoB * longPi_f)));
  Y_m = cumulative_sum(append_row(Zero, longLambda_m .* (S_m[:DeltaA*A] ./ X_m[:DeltaA*A]) .* (LongOnes - exp(-(longLambda_m - probLtoB * longPi_m))) ./ (longLambda_m - probLtoB * longPi_m)));

  L_f = X_f .* (Y_f + LlongOnes);
  L_m = X_m .* (Y_m + LlongOnes);

  B_f = LlongOnes - S_f - L_f;
  B_m = LlongOnes - S_m - L_m;

  /* new model (compared to mvb17) that splits between infectiousness from L vs B and enables estimation             */
  /* of all parameters. It still may not entirely be biologically intuitive. A full model which distinguishes        */
  /* between infectiousness after reactivation and re-infection and allows cycling in L and B would need at least    */
  /* five infected classes.                                                                                          */
  for (a in 1 : A) {
    aggr_S_f[a] = sum(longLambda_f[1+DeltaA*(a-1):DeltaA*a] .* S_f[1+DeltaA*(a-1):DeltaA*a]); // needed for cCMV
    aggr_L_f[a] = sum((rho_f[1+DeltaA*(a-1):DeltaA*a] + z * longLambda_f[1+DeltaA*(a-1):DeltaA*a]) .* L_f[1+DeltaA*(a-1):DeltaA*a]); 
    aggr_L_m[a] = sum((rho_m[1+DeltaA*(a-1):DeltaA*a] + z * longLambda_m[1+DeltaA*(a-1):DeltaA*a]) .* L_m[1+DeltaA*(a-1):DeltaA*a]);
    aggr_B_f[a] = sum((rho_f[1+DeltaA*(a-1):DeltaA*a] + z * longLambda_f[1+DeltaA*(a-1):DeltaA*a]) .* B_f[1+DeltaA*(a-1):DeltaA*a]); 
    aggr_B_m[a] = sum((rho_m[1+DeltaA*(a-1):DeltaA*a] + z * longLambda_m[1+DeltaA*(a-1):DeltaA*a]) .* B_m[1+DeltaA*(a-1):DeltaA*a]);
  }
  lambda_hat_f = Contact_FF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + reducinf * beta2 * aggr_L_f + beta2 * aggr_B_f) 
               + Contact_FM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + reducinf * beta2 * aggr_L_m + beta2 * aggr_B_m);
  lambda_hat_m = Contact_MM * (beta1 * (S_m[ReduceIdxs] - S_m[ReduceIdxsRightShift]) + reducinf * beta2 * aggr_L_m + beta2 * aggr_B_m) 
               + Contact_MF * (beta1 * (S_f[ReduceIdxs] - S_f[ReduceIdxsRightShift]) + reducinf * beta2 * aggr_L_f + beta2 * aggr_B_f);
  
  /* vertical transmission - i.e. serologically positive at 6 months of age */
  for (a in 1 : A) { 
	S_birth[a] = BirthContribution[a] * (0.5*(S_f[(a-1)*DeltaA+1]+S_f[a*DeltaA+1]) + (1-nu)*(1-0.5*(S_f[(a-1)*DeltaA+1]+S_f[a*DeltaA+1])));
  }
  frac_S_birth = sum(S_birth);
  
  /* congenital infection - assumed to born by acutely infected mothers */
  /* acute infection can be caused by primary or secondary infection,   */
  /* or by reactivation                                                 */
  pcCMV = qcCMV * sum(BirthContribution .* ((aggr_S_f + aggr_L_f + aggr_B_f) / DeltaA)); 
}

model {    
  /* spline weights for the reactivation rates */
  /* As discussed in Rozhnova et al (2019), priors for the weight do have a strong impact on the parameter estimates */
  /* in the new, flexible model that allows for multiple re-infection/reactivation events over a person's life.      */
  /* For instance, with the alternative weights for the priors, the estimated reactivation rates are much higher     */
  /* Hoewever, there is a concomittant reduction in the infectivity of re-infection/reactivation, such that overall  */
  /* estimated age-specific prevalence in females and males remain virtually identical. As shown in Rozhnova et al   */
  /* (2019) the impact of interventions (e.g., vaccination) is also braodly unaffected                               */
  for (s in 1:num_basis) {
	 a_raw[,s] ~  gamma(2, 50);                       // default scenario based on mvb17 estimates and the premise that reactivation is rare
	 //a_raw[,s] ~  gamma(10, 20);                    // alternative scenario, assuming that reactivation may occur frequently
  }
  
  /* penalise the difference between lambda and lambda_hat and S0 and S0_hat to solve the nonlinear equations           */
  /* TODO: lambda_f - lambda_hat_f ~ normal(0,1/Penalty) does not work: do something like target += log(det(Jacobian))  */
  /* future: use solver. I have tried this but now still seems prohibitively slow                                       */
  /* the formulation below can be viewed as a model in itself (a priori we would like to select parameter               */
  /* values that are compatible with a transmission model)                                                              */
  lambda_f ~ normal(lambda_hat_f, 1/Penalty);                
  lambda_m ~ normal(lambda_hat_m, 1/Penalty);                
  S0 ~ normal(frac_S_birth, 1/Penalty);                      
  
  /* data set 1: serological data from PIENTER2 study */
  for ( i in 1 : N ) {// loop over subjects
     int aa;
     real pS; real pL; real pB;
        
     // improve readability
	 aa = Ages[i] + 1; // the index for S, L and B
        
     // compute the compartment-probabilities given the subjects' age
     if ( Gender[i] == 0 ) { // 0 = female
        pS = S_f[aa];
        pL = L_f[aa];
        pB = B_f[aa];
     }
     else { // 1 = male
        pS = S_m[aa];
        pL = L_m[aa];
        pB = B_m[aa];
     }
        
     // likelihood contributions
     if ( Censor[i] == 0 ) { // normal data
       target += watanabe_beta * log( pS * exp(normal_lpdf(Titers[i] | MuS, SigmaS)) +
                                      pL * exp(normal_lpdf(Titers[i] | MuL, SigmaL)) +
                                      pB * exp(normal_lpdf(Titers[i] | MuB, SigmaB)) ); //or use log_exp_sum 
     }
     else if ( Censor[i] == 1 ) { // right censored
            target += watanabe_beta * log( pS * exp(normal_lccdf(RightCensor | MuS, SigmaS)) +
                                           pL * exp(normal_lccdf(RightCensor | MuL, SigmaL)) +
                                           pB * exp(normal_lccdf(RightCensor | MuB, SigmaB)) );
          }
          else if ( Censor[i] == 2 ) { // spike
                 target += watanabe_beta * log(pS);
               }
  }
  
  /* data set 2: cCMV cases from the CROCUS study (PhD thesis M Korndewal) */
  target += watanabe_beta * binomial_lpmf(numbercCMVinfants | numbertestedinfants, pcCMV); 
}

generated quantities {
    vector[N] log_lik;
    real log_like;

    /* for WAIC/WBIC calculation */
	for ( i in 1:N ) {
		int aa;
        real pS; real pL; real pB;

        // make the code a little bit more readable
		aa = Ages[i] + 1; // the index for S, L and B
        
        // compute the compartment-probabilities given the subjects' age
        if ( Gender[i] == 0 ) { // 0 = female
            pS = S_f[aa];
            pL = L_f[aa];
            pB = B_f[aa];
        }
        else{ // 1 = male
            pS = S_m[aa];
            pL = L_m[aa];
            pB = B_m[aa];
        }
        
        // define the likelihood
        if ( Censor[i] == 0 ) { // normal data
            log_lik[i] = log( pS * exp(normal_lpdf(Titers[i] | MuS, SigmaS)) +
                              pL * exp(normal_lpdf(Titers[i] | MuL, SigmaL)) +
                              pB * exp(normal_lpdf(Titers[i] | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 1 ) { // right censored
            log_lik[i] = log( pS * exp(normal_lccdf(RightCensor | MuS, SigmaS)) +
                              pL * exp(normal_lccdf(RightCensor | MuL, SigmaL)) +
                              pB * exp(normal_lccdf(RightCensor | MuB, SigmaB)) );
        }
        else if ( Censor[i] == 2 ) { // spiked
            log_lik[i] = log(pS);
        }
    }
    log_like = sum(log_lik);
}
