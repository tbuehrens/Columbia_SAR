data{
  int data_yrs;
  int forecast_yrs;
  int pops; 
  int obs;
  int obs_yr[obs];
  int obs_age[obs];
  int obs_pop[obs];
  int det[obs];
  int rel[obs];
}
transformed data{
    int yrs = data_yrs + forecast_yrs;
}
parameters{
  vector<lower=0,upper=1>[pops] SAR_0;
  real<lower=0> sigma_SAR;
  matrix[yrs-1,pops] eps_SAR;
  vector<lower=0,upper=1>[pops] p_OA1_0;
  real<lower=0> sigma_p_OA1;
  matrix[yrs-1,pops] eps_p_OA1;
}
transformed parameters{
  matrix<lower=0,upper=1>[yrs,pops] SAR;
  matrix<lower=0,upper=1>[yrs,pops] p_OA1;
  SAR[1,] = to_row_vector(SAR_0);
  p_OA1[1,] = to_row_vector(p_OA1_0);
  
  for(y in 2:yrs){
    SAR[y,] = inv_logit(logit(SAR[y-1,]) + eps_SAR[y-1,] * sigma_SAR);
    p_OA1[y,] = inv_logit(logit(p_OA1[y-1,]) + eps_p_OA1[y-1,] * sigma_p_OA1);
  }
}
model{
  //priors
  SAR_0 ~ beta(1,1);
  sigma_SAR ~ std_normal();
  to_vector(eps_SAR) ~ std_normal();
  p_OA1_0 ~ beta(1,1);
  sigma_p_OA1 ~ std_normal();
  to_vector(eps_p_OA1) ~ std_normal();
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
}