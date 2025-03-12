functions{
  matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) { 
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-(square(x[i,j])*(sq_rho/3)));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
}
data{
  int N;
  int N2;
  int N_pred;
  int NFam;
  int NParam;
  
  int PI_notmis_len;
  int GE_notmis_len;
  int AP_notmis_len;
  int DI_notmis_len;
  int TT_notmis_len;
  int CS_notmis_len;

  int PI_mis_len;
  int GE_mis_len;
  int AP_mis_len;
  int DI_mis_len;
  int TT_mis_len;
  int CS_mis_len;

  array[N] int F;
  matrix[N2,N2] M;

  //means
  array[5] real means;
  array[5] real sds;
  
  //Counts
  array[N] int GR;
  array[N] int PL;
  array[N] int PI;
  array[N] int CL;
  array[N] int GE;
  array[N] int AP;
  array[N] int QU;
  array[N] int AN;
  array[N] int TV;
  array[N] int DI;
  array[N] int TT;
  array[N] int IS;
  array[N] int FS;

  array[PI_notmis_len] int PI_notmis;
  array[GE_notmis_len] int GE_notmis;
  array[AP_notmis_len] int AP_notmis;
  array[DI_notmis_len] int DI_notmis;
  array[TT_notmis_len] int TT_notmis;
  array[CS_notmis_len] int CS_notmis;

  array[PI_mis_len] int PI_mis;
  array[GE_mis_len] int GE_mis;
  array[AP_mis_len] int AP_mis;
  array[DI_mis_len] int DI_mis;
  array[TT_mis_len] int TT_mis;
  array[CS_mis_len] int CS_mis;
  
  //Binomial
  array[N] int CS;

  array[N] int ObsID;
  array[N_pred] int PredID;
}
parameters{
  vector<lower=0>[NParam] alphasq;
  vector<lower=0>[NParam] rhosq;
  matrix[NFam, NParam] a_z;
  vector[NParam] a_mu;
  vector<lower=0>[NParam] a_sigma;
  matrix[N2, NParam] eta;
  cholesky_factor_corr[NParam] SIGMA_MV;
}
transformed parameters{
  matrix[NFam, NParam] a;
  for(i in 1:NParam){
    a[,i] = a_mu[i] + a_z[,i] .* a_sigma[i];
  }
  
  matrix[N2, NParam] f;
  matrix[N, NParam] lambda;
  matrix[N_pred, NParam] lambda_pred;
{
  array[NParam] matrix[N2,N2] SIGMA;
  array[NParam] matrix[N2, N2] L_K;
  for(i in 1:NParam){
    SIGMA[i] = cov_GPL2(M, alphasq[i], rhosq[i], 0.001);
    L_K[i] = cholesky_decompose(SIGMA[i]);
    f[,i] = L_K[i] * eta[,i];
  }
}
  for(i in 1:13){  
    lambda[,i] = exp(a[F,i] + f[ObsID,i]);
    lambda_pred[,i] = exp(a_mu[i] + f[PredID,i]);
  }
  lambda[,14] = a[F,14] + f[ObsID,14];
  lambda_pred[,14] = a_mu[14] + f[PredID,14];
}
model{
  rhosq ~ exponential(2);
  alphasq ~ exponential(0.2);
  for(i in 1:N2){
    eta[i,] ~ multi_normal_cholesky(rep_vector(0,NParam), diag_pre_multiply(rep_vector(1,NParam), SIGMA_MV));
  }

  SIGMA_MV ~ lkj_corr_cholesky(2);


  to_vector(a_z) ~ normal(0,1);
  a_mu ~ normal(0,1);
  a_sigma ~ exponential(1);


  GR ~ poisson(lambda[1:N,1]);
  PL ~ poisson(lambda[1:N,2]);
  PI[PI_notmis] ~ poisson(lambda[PI_notmis,3]);
  CL ~ poisson(lambda[1:N,4]);
  GE[GE_notmis] ~ poisson(lambda[GE_notmis,5]);
  AP[AP_notmis] ~ poisson(lambda[AP_notmis,6]);
  QU ~ poisson(lambda[1:N,7]);
  AN ~ poisson(lambda[1:N,8]); 
  TV ~ poisson(lambda[1:N,9]);
  DI[DI_notmis] ~ poisson(lambda[DI_notmis,10]);
  TT[TT_notmis] ~ poisson(lambda[TT_notmis,11]);
  IS ~ poisson(lambda[1:N,12]);
  FS ~ poisson(lambda[1:N,13]);
  CS[CS_notmis] ~ bernoulli_logit(lambda[CS_notmis,14]);

  lambda[PI_mis,3] ~ normal(means[1], sds[1]);
  lambda[GE_mis,5] ~ normal(means[2], sds[2]);
  lambda[AP_mis,6] ~ normal(means[3], sds[3]);
  lambda[DI_mis,10] ~ normal(means[4], sds[4]);
  lambda[TT_mis,11] ~ normal(means[5], sds[5]);
  lambda[CS_mis,14] ~ normal(0, 3);

      
}
