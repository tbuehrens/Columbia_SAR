data{
  int data_yrs; //years of PIT tag data
  int forecast_yrs; //years to forecast SAR
  int pops; //number of different populations we have PIT tag data for
  int obs; //number of PIT tag SAR observations 
  int obs_yr[obs]; //year index for PIT tag SAR observations
  int obs_age[obs]; //age index for PIT tag SAR observations
  int obs_pop[obs]; //population index for PIT tag SAR observations
  int det[obs]; //number of smolt PIT tags detected passing BON going downstream
  int rel[obs]; //number of adult PIT tags detected passing BON going upstream
  int ry_adults_obs[data_yrs-1,pops]; //Bonneville Dam adult count
  vector[data_yrs + forecast_yrs] NPGO; //NPGO for outmigration years and forecast years
  vector[data_yrs + forecast_yrs] X2; //second covariate
}
transformed data{
  int yrs = data_yrs + forecast_yrs;
}
parameters{
  vector[pops] b1;
  vector[pops] b2;

  vector<lower=0,upper=1>[pops] mu_SAR;
  vector<lower=-1,upper=1>[pops] phi_SAR;
  vector<lower=0>[pops] sigma_SAR;
  matrix[yrs,pops] eps_SAR;

  vector<lower=0,upper=1>[pops] mu_p_OA1;
  vector<lower=-1,upper=1>[pops] phi_p_OA1;
  vector<lower=0>[pops] sigma_p_OA1;
  matrix[yrs,pops] eps_p_OA1;
  
  vector<lower=0>[pops] smolts_0;
  matrix[yrs-1,pops] eps_smolts;
  
  vector<lower=0>[pops]sigma_tot;
  vector<lower=0,upper=1>[pops] p_smolts;

}
transformed parameters{
  matrix<lower=0,upper=1>[yrs,pops] SAR;
  matrix[yrs,pops] SAR_resid;
  matrix<lower=0,upper=1>[yrs,pops] p_OA1;
  matrix[yrs,pops] p_OA1_resid;
  matrix<lower=0>[yrs,pops] smolts;
  matrix<lower=0>[yrs-1,pops] ry_adults;
  matrix<lower=0,upper=1>[yrs-1,pops] ry_p_OA1;
  vector<lower=0>[pops] sigma_disp;
  vector<lower=0>[pops] sigma_smolts;
  
  for(p in 1:pops){
    sigma_smolts[p] = sqrt(p_smolts[p] * square(sigma_tot[p]));
    sigma_disp[p] = sqrt((1-p_smolts[p]) * square(sigma_tot[p]));
  }
  
  SAR_resid[1,] =  eps_SAR[1,] .* to_row_vector(sigma_SAR[1:pops]) ./ sqrt(1-to_row_vector(square(phi_SAR)));
  p_OA1_resid[1,] =  eps_p_OA1[1,] .* to_row_vector(sigma_p_OA1[1:pops]) ./ sqrt(1-to_row_vector(square(phi_p_OA1)));
  smolts[1,] = to_row_vector(smolts_0 .* rep_vector(1e7,pops));
  
  for(y in 2:yrs){
    SAR_resid[y,1:pops] = to_row_vector(phi_SAR[1:pops]) .* SAR_resid[1,1:pops] + eps_SAR[y,1:pops] .* to_row_vector(sigma_SAR[1:pops]);
    p_OA1_resid[y,1:pops] = to_row_vector(phi_p_OA1[1:pops]) .* p_OA1_resid[1,1:pops] + eps_p_OA1[y,1:pops] .* to_row_vector(sigma_p_OA1[1:pops]);
    smolts[y,] = exp(log(smolts[y-1,]) + eps_smolts[y-1,] .* to_row_vector(sigma_smolts));
  }
  
  for(p in 1:pops){
    SAR_resid[1:yrs,p] = SAR_resid[1:yrs,p] - mean(SAR_resid[1:yrs,p]);//constrain SAR resids to have zero mean
    p_OA1_resid[1:yrs,p] = p_OA1_resid[1:yrs,p] - mean(p_OA1_resid[1:yrs,p]);//constrain p_OA1 resids to have zero mean

    SAR[1:yrs,p] = inv_logit(logit(mu_SAR[p]) + b1[p] * NPGO[1:yrs] + b2[p] * X2[1:yrs] + SAR_resid[1:yrs,p]);
    p_OA1[1:yrs,p] = inv_logit(logit(mu_p_OA1[p]) + p_OA1_resid[1:yrs,p]);
  }
  
  for(y in 1:(yrs-1)){
    for(p in 1:pops){
     ry_adults[y,p] = SAR[y+1,p] * p_OA1[y+1,p] * smolts[y+1,p] + SAR[y,p] * (1-p_OA1[y,p]) * smolts[y,p];
     ry_p_OA1[y,p] = SAR[y+1,p] * p_OA1[y+1,p] * smolts[y+1,p] / ry_adults[y,p];
    }
  }
}
model{
  //priors
  b1 ~ normal(0,10);
  b2 ~ normal(0,10);


  mu_SAR ~ beta(1,1);
  phi_SAR ~ std_normal();
  sigma_SAR ~ std_normal();
  to_vector(eps_SAR) ~ std_normal();

  mu_p_OA1 ~ beta(1,1);
  phi_p_OA1 ~ std_normal();
  sigma_p_OA1 ~ std_normal();
  to_vector(eps_p_OA1) ~ std_normal();
  
  smolts_0 ~ lognormal(0,0.5);
  to_vector(eps_smolts) ~ std_normal();
  
  sigma_tot ~ std_normal();
  p_smolts ~ beta(1,1);
  //likelihoods
  for(i in 1:obs){
    if(obs_age[i]==1){
      det[i] ~ poisson(rel[i] * SAR[obs_yr[i],obs_pop[i]] * p_OA1[obs_yr[i],obs_pop[i]]);
    }
    if(obs_age[i]==2){
      det[i] ~ poisson(rel[i] * SAR[obs_yr[i],obs_pop[i]] * (1-p_OA1[obs_yr[i],obs_pop[i]]));
    }
    if(obs_age[i]==99){
      det[i] ~ poisson(rel[i] * (1-SAR[obs_yr[i],obs_pop[i]]));
    }
  }
  for(y in 1:(data_yrs-1)){
    for(p in 1:pops){
      //ry_adults_obs[y,p] ~ neg_binomial(ry_adults[y,p],1/square(sigma_disp));
      ry_adults_obs[y,p] ~ lognormal(log(ry_adults[y,p]),sigma_disp[p]);
      //ry_adults_obs[y,p] ~ poisson(ry_adults[y,p]);
    }
  }
}
generated quantities{
  //int ry_adults_pred[yrs-1,pops];
  matrix[yrs-1,pops] ry_adults_pred;
  for(y in 1:(yrs-1)){
    for(p in 1:pops){
      //ry_adults_pred[y,p] = neg_binomial_2_rng(ry_adults[y,p],1/square(sigma_disp));
      ry_adults_pred[y,p] = lognormal_rng(log(ry_adults[y,p]),sigma_disp[p]);
      //ry_adults_pred[y,p] = poisson_rng(ry_adults[y,p]);
    }
  }
}
