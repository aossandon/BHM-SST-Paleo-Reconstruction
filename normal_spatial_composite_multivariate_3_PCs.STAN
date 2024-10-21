functions {
  
  matrix exp_cov_matrix (matrix x, real sill, real range, real nugget){

    int S;
    matrix[dims(x)[1],dims(x)[1]] cov;
    real neg_inv_range;

    S = dims(x)[1];
    neg_inv_range = -1.0/range;

    // off-diagonal elements
    for (i in 1:(S-1)) {
      cov[i,i] = sill + nugget;
      for (j in (i+1):S) {
      
        cov[i,j] = sill * exp(neg_inv_range * x[i,j]);
        cov[j,i] = cov[i,j];
      }
    }
    cov[S,S] = sill + nugget;

    return cov;

  }
 
   matrix exp_cov_matrix_2 (matrix x, vector scale, vector nugget, real range){

    int S;
    real neg_inv_range;
    matrix[dims(x)[1],dims(x)[1]] cov;

    S = dims(x)[1];

    neg_inv_range = -1.0/range;

    // off-diagonal elements
    for (i in 1:(S-1)) {
      cov[i,i]=square(scale[i])+square(nugget[i]);
      for (j in (i+1):S) {

        cov[i,j] = scale[i]*scale[j]*exp(neg_inv_range * x[i,j]);
        cov[j,i] = cov[i,j];
      }
    }
    cov[S,S] = square(scale[S])+square(nugget[S]);

    return cov;

  } 
  
 

  
  // real eliptical_copula_log_v1 (matrix covar, matrix aux_mu, matrix aux_data){
  //   vector[dims(aux_data)[1]] lp;
  // 
  //   int T; 
  //   T = dims(aux_data)[1]; 
  // 
  //   for(t in 1:T) lp[t] = multi_normal_cholesky_lpdf(aux_data[t,] | aux_mu[t,], covar);
  //   return sum(lp);  
  // }

}

data {
  int<lower=1> S;   // number of locations (grid points) for fitting the model
  int<lower=1> S_t;   // number of total locations (grid points)
  int<lower=1> T;   // number of data points (modern years) for the calibration
  int<lower=1> T_t;   // number of data points (modern years) for the validation
  int<lower=1> T_p;   // number of data points (paleo years) for the paleo reconstruction
  int<lower=1> G;   // number of total composite groups
  int<lower=1> R;   // number of elements in each group
  int<lower=1> P;  // number of space predictors (lat and long)
  int<lower=1> K;   // number of regression coefficient knots
  int<lower=1> PCS ; //number of PCs


  matrix[T,S] y;    // observations for the calibration
  matrix[S,S] dist; //distance matrix for grid points for fitting the model
  matrix[S,K] dist_knots_to_points; 
  matrix[S_t,S_t] dist_t; //distance matrix for grid points (all)
  matrix[S_t,K] dist_knots_to_points_t;// (all)

  matrix[S,P] covars; // space predictors
  matrix[S_t,P] covars_t; // space predictors
  matrix[T,PCS] covarsNS; // PCs for the calibration
  matrix[T_t,PCS] covarsNS_t; // PCs for the validation
  matrix[T_p,PCS] covarsNSp; // PCs for the paleo reconstruction
  
}

transformed data {

  matrix[P,1] ones_vector_p;

  int start;
  int end;
  int a;
  int b;

  array[G] matrix[R,R] dist_mats;
  matrix[S,K] dist_knots_to_points_sq;
  matrix[S_t,K] dist_knots_to_points_sq_t;
  
  ones_vector_p[,1]=rep_vector(1.0,P);

  // compute smaller distance matricies based on groups
  for(g in 1:G){
    start = (g-1) * R + 1;
    end = g * R;
    a = 0;
    for(i in start:end){
      a = a + 1;
      b = 0;
      for(j in start:end){
        b = b + 1;
        dist_mats[g,a,b] = dist[i,j];
      }
    }
  }
  
  
  for(k in 1:K){
    dist_knots_to_points_sq[,k] = col(dist_knots_to_points,k) .* col(dist_knots_to_points,k);
    dist_knots_to_points_sq_t[,k] = col(dist_knots_to_points_t,k) .* col(dist_knots_to_points_t,k);
  }

}

parameters {

  vector<lower=-50,upper=50>[S] loc;
  vector<lower=-50,upper=50>[S] loc1;
  vector<lower=-50,upper=50>[S] loc2;
  vector<lower=-50,upper=50>[S] loc3;
  vector<lower=-10,upper=10>[S] log_scale;
  vector<lower=-10,upper=10>[S] log_nugget;
  

  matrix<lower=-50,upper=50>[P,K] basis_loc;
  matrix<lower=-50,upper=50>[P,K] basis_loc1;
  matrix<lower=-50,upper=50>[P,K] basis_loc2;
  matrix<lower=-50,upper=50>[P,K] basis_loc3;
  matrix<lower=-1,upper=1>[P,K] basis_log_scale;
  matrix<lower=-1,upper=1>[P,K] basis_log_nugget;

  real<lower=-100,upper=100> intercept_loc;
  real<lower=-100,upper=100> intercept_loc1;
  real<lower=-100,upper=100> intercept_loc2;
  real<lower=-100,upper=100> intercept_loc3;
  real<lower=-1,upper=1> intercept_log_scale;
  real<lower=-1,upper=1> intercept_log_nugget;

  real<lower=0,upper=1> sill_loc;
  real<lower=0,upper=1> sill_loc1;
  real<lower=0,upper=1> sill_loc2;
  real<lower=0,upper=1> sill_loc3;
  real<lower=0,upper=1> sill_log_scale;
  real<lower=0,upper=1> sill_log_nugget;

  real<lower=0,upper=1> nugget_loc;
  real<lower=0,upper=1> nugget_loc1;
  real<lower=0,upper=1> nugget_loc2;
  real<lower=0,upper=1> nugget_loc3;
  real<lower=0,upper=1> nugget_log_scale;
  real<lower=0,upper=1> nugget_log_nugget;

  real<lower=0> range;
  real<lower=0> range_loc;
  real<lower=0> range_loc1;
  real<lower=0> range_loc2;
  real<lower=0> range_loc3;
  real<lower=0> range_log_scale;
  real<lower=0> range_log_nugget;

  matrix<lower=0>[P,K] range_basis_loc;
  matrix<lower=0>[P,K] range_basis_loc1;
  matrix<lower=0>[P,K] range_basis_loc2;
  matrix<lower=0>[P,K] range_basis_loc3;
  matrix<lower=0>[P,K] range_basis_log_scale;
  matrix<lower=0>[P,K] range_basis_log_nugget;

} 


model {    


  real LogLik;
  matrix[R,R] covariance;
  matrix[T,S] mu;
  matrix[T,R] aux_mu;
  matrix[T,R] aux_data;

  //vector of LogLIk
  vector[T] lp;
  vector[G] log_lik;

  array[R] int ind;


  vector[R] loc_mean_group;
  vector[R] loc1_mean_group;
  vector[R] loc2_mean_group;
  vector[R] loc3_mean_group;
  vector[R] log_scale_mean_group;
  vector[R] log_nugget_mean_group;

  matrix[S,1] loc_mean;
  matrix[S,1] loc1_mean;
  matrix[S,1] loc2_mean;
  matrix[S,1] loc3_mean;
  matrix[S,1] log_scale_mean;
  matrix[S,1] log_nugget_mean;

  matrix[S,P] beta_loc;
  matrix[S,P] beta_loc1;
  matrix[S,P] beta_loc2;
  matrix[S,P] beta_loc3;
  matrix[S,P] beta_log_scale;
  matrix[S,P] beta_log_nugget;

  matrix[R,R] cov_loc;
  matrix[R,R] cov_loc1;
  matrix[R,R] cov_loc2;
  matrix[R,R] cov_loc3;
  matrix[R,R] cov_log_scale;
  matrix[R,R] cov_log_nugget;

  int start2;
  int end2;
  int m;
  matrix[P,K] neg_inv_range_loc;
  matrix[P,K] neg_inv_range_loc1;
  matrix[P,K] neg_inv_range_loc2;
  matrix[P,K] neg_inv_range_loc3;
  matrix[P,K] neg_inv_range_log_scale;
  matrix[P,K] neg_inv_range_log_nugget;

  // Priors, fairly diffuse
  intercept_loc ~   normal(0, 10);
  intercept_loc1 ~   normal(0, 10);
  intercept_loc2 ~   normal(0, 10);
  intercept_loc3 ~   normal(0, 10);
  intercept_log_scale ~ normal(0, 1);
  intercept_log_nugget ~ normal(0, 1);

  range ~ uniform(0, 1000);
  range_loc ~    uniform(0, 1000);
  range_loc1 ~    uniform(0, 1000);
  range_loc2 ~    uniform(0, 1000);
  range_loc3 ~    uniform(0, 1000);
  range_log_scale ~  uniform(0, 1000);
  range_log_nugget ~  uniform(0, 1000);

  sill_loc ~   normal(0, 1);
  sill_loc1 ~   normal(0, 1);
  sill_loc2 ~   normal(0, 1);
  sill_loc3 ~   normal(0, 1);
  sill_log_scale ~ normal(0, 1);
  sill_log_nugget ~ normal(0, 1);

  nugget_loc ~   normal(0, 1);
  nugget_loc1 ~   normal(0, 1);
  nugget_loc2 ~   normal(0, 1);
  nugget_loc3 ~   normal(0, 1);
  nugget_log_scale ~ normal(0, 1);
  nugget_log_nugget ~ normal(0, 1);
  
  for (k in 1:K){
    for(p in 1:P){
      basis_loc[p,k] ~   normal(0, 10);
      basis_loc1[p,k] ~   normal(0, 10);
      basis_loc2[p,k] ~   normal(0, 10);
      basis_loc3[p,k] ~   normal(0, 10);
      basis_log_scale[p,k] ~ normal(0, 1);
      basis_log_nugget[p,k] ~ normal(0, 1);
  
      range_basis_loc[p,k] ~   uniform(0, 1000);
      range_basis_loc1[p,k] ~   uniform(0, 1000);
      range_basis_loc2[p,k] ~   uniform(0, 1000);
      range_basis_loc3[p,k] ~   uniform(0, 1000);
      range_basis_log_scale[p,k] ~ uniform(0, 1000);
      range_basis_log_nugget[p,k] ~ uniform(0, 1000);
    }
  }

  
  for(k in 1:K){
      neg_inv_range_loc[,k] = -1.0./(col(range_basis_loc,k) .* col(range_basis_loc,k));
      neg_inv_range_loc1[,k] = -1.0./(col(range_basis_loc1,k) .* col(range_basis_loc1,k));
      neg_inv_range_loc2[,k] = -1.0./(col(range_basis_loc2,k) .* col(range_basis_loc2,k));
      neg_inv_range_loc3[,k] = -1.0./(col(range_basis_loc3,k) .* col(range_basis_loc3,k));
      neg_inv_range_log_scale[,k] = -1.0./(col(range_basis_log_scale,k) .* col(range_basis_log_scale,k));
      neg_inv_range_log_nugget[,k] = -1.0./(col(range_basis_log_nugget,k) .* col(range_basis_log_nugget,k));
  }

  // compute regression coefficients from basis functions
  for(s in 1:S){
    for(p in 1:P){
      // zero out betas
      beta_loc[s,p] = 0.0;
      beta_loc1[s,p] = 0.0;
      beta_loc2[s,p] = 0.0;
      beta_loc3[s,p] = 0.0;
      beta_log_scale[s,p] = 0.0;
      beta_log_nugget[s,p] = 0.0;

      // add up all the basis funcitons (gaussian kernals)
      for(k in 1:K){
        beta_loc[s,p] = beta_loc[s,p] + basis_loc[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_loc[p,k]);
        beta_loc1[s,p] = beta_loc1[s,p] + basis_loc1[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_loc1[p,k]);
        beta_loc2[s,p] = beta_loc2[s,p] + basis_loc2[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_loc2[p,k]);
        beta_loc3[s,p] = beta_loc3[s,p] + basis_loc3[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_loc3[p,k]);
        beta_log_scale[s,p] = beta_log_scale[s,p] + basis_log_scale[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_log_scale[p,k]);
        beta_log_nugget[s,p] = beta_log_nugget[s,p] + basis_log_nugget[p,k] * exp(dist_knots_to_points_sq[s,k] * neg_inv_range_log_nugget[p,k]);
      }
    }
  }

  // right multiplying by a vector of ones will sum the rows  
  loc_mean = intercept_loc + (covars .* beta_loc) * ones_vector_p;
  loc1_mean = intercept_loc1 + (covars .* beta_loc1) * ones_vector_p;
  loc2_mean = intercept_loc2 + (covars .* beta_loc2) * ones_vector_p;
  loc3_mean = intercept_loc3 + (covars .* beta_loc3) * ones_vector_p;
  log_scale_mean = intercept_log_scale + (covars .* beta_log_scale) * ones_vector_p;
  log_nugget_mean = intercept_log_nugget + (covars .* beta_log_nugget) * ones_vector_p;
  // loop over groups
  for(g in 1:G){

    start2 = (g-1) * R + 1;
    end2 = g * R;

    m = 0;
    for(i in start2:end2){
      m = m + 1;
      ind[m] = i;

      loc_mean_group[m] = loc_mean[i,1];
      loc1_mean_group[m] = loc1_mean[i,1];
      loc2_mean_group[m] = loc2_mean[i,1];
      loc3_mean_group[m] = loc3_mean[i,1];
      log_scale_mean_group[m] = log_scale_mean[i,1];
      log_nugget_mean_group[m] = log_nugget_mean[i,1];
     
      // loc_group[m] = loc[i];
      // loc1_group[m] = loc1[i];
      // loc2_group[m] = loc2[i];
      // log_scale_group[m] = log_scale[i];
    }

    cov_loc = exp_cov_matrix(dist_mats[g], sill_loc, range_loc, nugget_loc);
    cov_loc1 = exp_cov_matrix(dist_mats[g], sill_loc1, range_loc1, nugget_loc1);
    cov_loc2 = exp_cov_matrix(dist_mats[g], sill_loc2, range_loc2, nugget_loc2);
    cov_loc3 = exp_cov_matrix(dist_mats[g], sill_loc3, range_loc3, nugget_loc3);
    cov_log_scale = exp_cov_matrix(dist_mats[g], sill_log_scale, range_log_scale, nugget_log_scale);
    cov_log_nugget = exp_cov_matrix(dist_mats[g], sill_log_nugget, range_log_nugget, nugget_log_nugget);

    loc[start2:end2] ~ multi_normal_cholesky(loc_mean_group, cholesky_decompose(cov_loc));
    loc1[start2:end2] ~ multi_normal_cholesky(loc1_mean_group, cholesky_decompose(cov_loc1));
    loc2[start2:end2] ~ multi_normal_cholesky(loc2_mean_group, cholesky_decompose(cov_loc2));
    loc3[start2:end2] ~ multi_normal_cholesky(loc3_mean_group, cholesky_decompose(cov_loc3));
    log_scale[start2:end2] ~ multi_normal_cholesky(log_scale_mean_group, cholesky_decompose(cov_log_scale));
    log_nugget[start2:end2] ~ multi_normal_cholesky(log_nugget_mean_group, cholesky_decompose(cov_log_nugget));
    

    // loop over stations in the group
    for(r in 1:R){

      mu[,ind[r]] = loc[ind[r]] + loc1[ind[r]]*col(covarsNS,1) + loc2[ind[r]]*col(covarsNS,2) + loc3[ind[r]]*col(covarsNS,3); 

      aux_mu[,r] = mu[,ind[r]];
      aux_data[,r] = y[,ind[r]];
    }

    
    covariance = exp_cov_matrix_2(dist_mats[g],exp(log_scale[start2:end2]),exp(log_nugget[start2:end2]), range);
    for(t in 1:T){
      lp[t] = multi_normal_cholesky_lpdf(aux_data[t,] | aux_mu[t,], cholesky_decompose(covariance));
    } 
    log_lik[g]=sum(lp);
  }
  
  target += sum(log_lik);
  
}


generated quantities{

  matrix [T_t,S_t] y_rep;
  matrix [T_p,S_t] y_repP;
  
  
  

  
  vector[T] log_lik;

  { matrix[T_t,S_t] mu_t;
    matrix[S_t,S_t] covariance_t;
    matrix[T_p,S_t] mu_p;
    matrix[S,S] covariance;
    matrix[T,S] mu;
    
    
    vector[S_t] loc_t;
    vector[S_t] loc1_t;
    vector[S_t] loc2_t;
    vector[S_t] loc3_t;
    vector[S_t] log_scale_t;
    vector[S_t] log_nugget_t;

    matrix[S_t,1] loc_mean_t;
    matrix[S_t,1] loc1_mean_t;
    matrix[S_t,1] loc2_mean_t;
    matrix[S_t,1] loc3_mean_t;
    matrix[S_t,1] log_scale_mean_t;
    matrix[S_t,1] log_nugget_mean_t;

    matrix[S_t,P] beta_loc_t;
    matrix[S_t,P] beta_loc1_t;
    matrix[S_t,P] beta_loc2_t;
    matrix[S_t,P] beta_loc3_t;
    matrix[S_t,P] beta_log_scale_t;
    matrix[S_t,P] beta_log_nugget_t;

    matrix[S_t,S_t] cov_loc_t;
    matrix[S_t,S_t] cov_loc1_t;
    matrix[S_t,S_t] cov_loc2_t;
    matrix[S_t,S_t] cov_loc3_t;
    matrix[S_t,S_t] cov_log_scale_t;
    matrix[S_t,S_t] cov_log_nugget_t;

    matrix[P,K] neg_inv_range_loc;
    matrix[P,K] neg_inv_range_loc1;
    matrix[P,K] neg_inv_range_loc2;
    matrix[P,K] neg_inv_range_loc3;
    matrix[P,K] neg_inv_range_log_scale;
    matrix[P,K] neg_inv_range_log_nugget;

    for(k in 1:K){
      neg_inv_range_loc[,k] = -1.0./(col(range_basis_loc,k) .* col(range_basis_loc,k));
      neg_inv_range_loc1[,k] = -1.0./(col(range_basis_loc1,k) .* col(range_basis_loc1,k));
      neg_inv_range_loc2[,k] = -1.0./(col(range_basis_loc2,k) .* col(range_basis_loc2,k));
      neg_inv_range_loc3[,k] = -1.0./(col(range_basis_loc3,k) .* col(range_basis_loc3,k));
      neg_inv_range_log_scale[,k] = -1.0./(col(range_basis_log_scale,k) .* col(range_basis_log_scale,k));
      neg_inv_range_log_nugget[,k] = -1.0./(col(range_basis_log_nugget,k) .* col(range_basis_log_nugget,k));
    }

    // compute regression coefficients from basis functions
    for(s in 1:S_t){
      for(p in 1:P){
        // zero out betas
        beta_loc_t[s,p] = 0.0;
        beta_loc1_t[s,p] = 0.0;
        beta_loc2_t[s,p] = 0.0;
        beta_loc3_t[s,p] = 0.0;
        beta_log_scale_t[s,p] = 0.0;
        beta_log_nugget_t[s,p] = 0.0;

        // add up all the basis funcitons (gaussian kernals)
        for(k in 1:K){
          beta_loc_t[s,p] = beta_loc_t[s,p] + basis_loc[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_loc[p,k]);
          beta_loc1_t[s,p] = beta_loc1_t[s,p] + basis_loc1[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_loc1[p,k]);
          beta_loc2_t[s,p] = beta_loc2_t[s,p] + basis_loc2[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_loc2[p,k]);
          beta_loc3_t[s,p] = beta_loc3_t[s,p] + basis_loc3[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_loc3[p,k]);
          beta_log_scale_t[s,p] = beta_log_scale_t[s,p] + basis_log_scale[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_log_scale[p,k]);
          beta_log_nugget_t[s,p] = beta_log_nugget_t[s,p] + basis_log_nugget[p,k] * exp(dist_knots_to_points_sq_t[s,k] * neg_inv_range_log_nugget[p,k]);
        }
      }
    }

    // right multiplying by a vector of ones will sum the rows
    loc_mean_t = intercept_loc + (covars_t .* beta_loc_t) * ones_vector_p;
    loc1_mean_t = intercept_loc1 + (covars_t .* beta_loc1_t) * ones_vector_p;
    loc2_mean_t = intercept_loc2 + (covars_t .* beta_loc2_t) * ones_vector_p;
    loc3_mean_t = intercept_loc3 + (covars_t .* beta_loc3_t) * ones_vector_p;
    log_scale_mean_t = intercept_log_scale + (covars_t .* beta_log_scale_t) * ones_vector_p;
    log_nugget_mean_t = intercept_log_nugget + (covars_t .* beta_log_nugget_t) * ones_vector_p;

    cov_loc_t = exp_cov_matrix(dist_t, sill_loc, range_loc, nugget_loc);
    cov_loc1_t = exp_cov_matrix(dist_t, sill_loc1, range_loc1, nugget_loc1);
    cov_loc2_t = exp_cov_matrix(dist_t, sill_loc2, range_loc2, nugget_loc2);
    cov_loc3_t = exp_cov_matrix(dist_t, sill_loc3, range_loc3, nugget_loc3);
    cov_log_scale_t = exp_cov_matrix(dist_t, sill_log_scale, range_log_scale, nugget_log_scale);
    cov_log_nugget_t = exp_cov_matrix(dist_t, sill_log_nugget, range_log_nugget, nugget_log_nugget);

    loc_t = multi_normal_cholesky_rng(loc_mean_t[,1], cholesky_decompose(cov_loc_t));
    loc1_t = multi_normal_cholesky_rng(loc1_mean_t[,1], cholesky_decompose(cov_loc1_t));
    loc2_t = multi_normal_cholesky_rng(loc2_mean_t[,1], cholesky_decompose(cov_loc2_t));
    loc3_t = multi_normal_cholesky_rng(loc3_mean_t[,1], cholesky_decompose(cov_loc3_t));
    log_scale_t = multi_normal_cholesky_rng(log_scale_mean_t[,1], cholesky_decompose(cov_log_scale_t));
    log_nugget_t = multi_normal_cholesky_rng(log_nugget_mean_t[,1], cholesky_decompose(cov_log_nugget_t));
    
    for(r in 1:S_t){
      mu_t[,r] = loc_t[r] + loc1_t[r]*col(covarsNS_t,1) + loc2_t[r]*col(covarsNS_t,2) + loc3_t[r]*col(covarsNS_t,3);
      mu_p[,r] = loc_t[r] + loc1_t[r]*col(covarsNSp,1) + loc2_t[r]*col(covarsNSp,2) + loc3_t[r]*col(covarsNSp,3);
      }
    covariance_t = exp_cov_matrix_2(dist_t,exp(log_scale_t),exp(log_nugget_t), range);
    
    for(r in 1:S){
      mu[,r] = loc[r] + loc1[r]*col(covarsNS,1) + loc2[r]*col(covarsNS,2) + loc3[r]*col(covarsNS,3);
      }
    covariance = exp_cov_matrix_2(dist,exp(log_scale),exp(log_nugget), range);
    for(t in 1:T){
      log_lik[t] = multi_normal_cholesky_lpdf(y[t,] | mu[t,], cholesky_decompose(covariance));
    }
    for(t in 1:T_t){
      y_rep[t,] = multi_normal_cholesky_rng(mu_t[t,], cholesky_decompose(covariance_t))';
    }
    for(t in 1:T_p){
      y_repP[t,] = multi_normal_cholesky_rng(mu_p[t,], cholesky_decompose(covariance_t))';
    }
  }// loop over stations in the group


}

