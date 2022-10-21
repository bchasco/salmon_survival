
// Anisotropic version of "spde.cpp".
#include <TMB.hpp>
#include <fenv.h> // Extra line needed

template <class Type>
Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); // Extra line needed
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // Indices
  DATA_INTEGER(j_flag);
  DATA_INTEGER(l_flag);
  DATA_INTEGER(t_flag);
  DATA_INTEGER(z_jl_flag);
  DATA_INTEGER(z_jlt_flag);
  DATA_INTEGER(sim); //The data likelihood
  DATA_INTEGER(cond_sim);//0 means include the random effects
  DATA_INTEGER( n_i );         // Total number of observations
  DATA_INTEGER( n_s_jl );        //Number of vertices
  // Data
  DATA_IVECTOR( x_s_jl );	      // Association of each station with a given vertex in SPDE mesh
  DATA_IVECTOR( s_i_jl );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IARRAY( j_shift );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IARRAY( l_shift );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IARRAY( management );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IVECTOR( t_i );	      // year
  DATA_IVECTOR( l_i );	      // length
  DATA_IVECTOR( j_i );	      // day
  DATA_IVECTOR( a_i );	      // bypass
  DATA_IVECTOR( j_im7 );	      // day
  DATA_IVECTOR( j_ip7 );	      // day
  DATA_IVECTOR( l_im10 );	      // length
  DATA_IVECTOR( l_ip10 );	      // length
  
  
  DATA_MATRIX(env);
  DATA_VECTOR(total); //total number of fish of observed
  DATA_VECTOR(surv); //data number of fish that survived
  DATA_STRUCT(spde_jl,spde_aniso_t);
  DATA_VECTOR(sim_lre);
  DATA_VECTOR(sim_lre_sd);
  DATA_VECTOR(sim_jre);
  DATA_VECTOR(sim_jre_sd);
  
  
  PARAMETER(mu);
  PARAMETER(bypass);
  PARAMETER(log_tau_jl);
  PARAMETER(log_tau_jl2);
  PARAMETER(log_kappa_jl);
  PARAMETER_VECTOR(ln_H_input_jl);
  PARAMETER_VECTOR(z_jl);
  PARAMETER_ARRAY(z_jlt);
  PARAMETER_ARRAY(y_re);
  PARAMETER_ARRAY(l_re);
  PARAMETER_ARRAY(j_re);
  PARAMETER(ln_sigl);
  PARAMETER(ln_sigy);
  PARAMETER(ln_sigj);
  PARAMETER(f_phiy);
  PARAMETER(f_phil);
  PARAMETER(f_phij);
  PARAMETER(f_psiy);
  PARAMETER(f_psil);
  PARAMETER(f_psij);
  PARAMETER(f_psijl);
  
  
  Type tau_jl = exp(log_tau_jl);
  Type tau_jl2 = exp(log_tau_jl2);
  Type kappa_jl = exp(log_kappa_jl);
  
  vector<Type> nll(10);
  nll.setZero();
  // Need to parameterize H matrix such that det(H)=1 (preserving volume) 
  // Note that H appears in (20) in Lindgren et al 2011
  matrix<Type> H_jl(2,2);

  if(z_jlt_flag==1 || z_jl_flag==1){
    H_jl(0,0) = exp(ln_H_input_jl(0));
    H_jl(1,0) = ln_H_input_jl(1);
    H_jl(0,1) = ln_H_input_jl(1);
    H_jl(1,1) = (1+ln_H_input_jl(1)*ln_H_input_jl(1)) / exp(ln_H_input_jl(0));
  }
  REPORT(H_jl);
  
  //Transform the AR1 parameters
  Type sig_l = exp(ln_sigl);
  Type sig_j = exp(ln_sigj);
  Type sig_y = exp(ln_sigy);
  Type phi_l = 1/(1+exp(f_phil));
  Type phi_j = 1/(1+exp(f_phij));
  Type phi_y = 1/(1+exp(f_phiy));
  Type psi_l = 1/(1+exp(f_psil));
  Type psi_j = 1/(1+exp(f_psij));
  Type psi_y = 1/(1+exp(f_psiy));
  Type psi_jl = 1/(1+exp(f_psijl));
  
  
  //2-D smoother
  SparseMatrix<Type> Q_jl = Q_spde(spde_jl,kappa_jl,H_jl);

  if(z_jl_flag==1){
    nll(0) += GMRF(Q_jl)(z_jl);
    if(cond_sim==0){
      GMRF(Q_jl).simulate(z_jl);
    }
  }
  vector<Type> ln_zexp_sp_jl( n_s_jl);
  ln_zexp_sp_jl.setZero();
  ln_zexp_sp_jl = z_jl/tau_jl;
  int n_t = 22;  

  //Scale the 2D variance  
  matrix<Type> ln_zexp_sp_jlt( n_s_jl, n_t);
  matrix<Type> jl_cov(2,2);
  jl_cov(0,0) = 1.;
  jl_cov(1,1) = 1.;
  jl_cov(0,1) = psi_jl;
  jl_cov(1,0) = psi_jl;
  MVNORM_t<Type> jl_dnorm(jl_cov); 
  if(z_jlt_flag==1){
    for(int t=0; t<n_t; t++){
      // nll(2) += GMRF(Q_jl)(z_jlt.col(t));
      nll(2) += SCALE(SEPARABLE(jl_dnorm, GMRF(Q_jl)), 1/tau_jl2)( z_jlt.col(t));
      // ln_zexp_sp_jlt.col(t) = z_jlt.col(t)/tau_jl2;
    }
    // jnll_pointer += SCALE(SEPARABLE(MVNORM(Cov_cc), gmrf_Q), exp(-logtau))( diff_gmrf_sc );
    if(cond_sim==0){
      // SIMULATE{
      //   vector<Type> omega_t(n_s_jl);
      //   for(int t = 0; t < n_t; t++){
      //     GMRF(Q_jl).simulate(omega_t);
      //     z_jlt.col(t) = omega_t;
      //     ln_zexp_sp_jlt.col(t) = z_jlt.col(t)/tau_jl2;
      //   }
      // }
    }
  }
  
  
  matrix<Type> y_cov(2,2);
  y_cov(0,0) = sig_y * sig_y;
  y_cov(1,1) = sig_y * sig_y;
  y_cov(0,1) = psi_y * sig_y * sig_y;
  y_cov(1,0) = psi_y * sig_y * sig_y;
  MVNORM_t<Type> y_dnorm(y_cov); 
  if(t_flag==1){
    nll(3) += AR1(phi_y,y_dnorm)(y_re.transpose());
  }   
  
  matrix<Type> l_cov(2,2);
  l_cov(0,0) = sig_l * sig_l;
  l_cov(1,1) = sig_l * sig_l;
  l_cov(0,1) = psi_l * sig_l * sig_l;
  l_cov(1,0) = psi_l * sig_l * sig_l;
  MVNORM_t<Type> l_dnorm(l_cov); 
  if(l_flag==1){
    nll(4) += AR1(phi_l,l_dnorm)(l_re.transpose());
  }
  
  
  matrix<Type> j_cov(2,2);
  j_cov(0,0) = sig_j * sig_j;
  j_cov(1,1) = sig_j * sig_j;
  j_cov(0,1) = psi_j * sig_j * sig_j;
  j_cov(1,0) = psi_j * sig_j * sig_j;
  MVNORM_t<Type> j_dnorm(j_cov); 
  if(j_flag==1){
    nll(5) += AR1(phi_j,j_dnorm)(j_re.transpose());
  }
  
  vector<Type> nu_i(n_i);
  vector<Type> eta_i(n_i);
  // matrix<Type> tmp(n_i,9);
  // tmp.setZero();
  // matrix<Type> tmp2(n_i,9);
  // tmp.setZero();
  // matrix<Type> proj_surv(y_re.size(),9);
  // matrix<Type> proj_surv_nox(y_re.size(),9);
  // vector<Type> proj_surv_ag(9);
  // proj_surv.setZero();
  // proj_surv_nox.setZero();
  // proj_surv_ag.setZero();
  // vector<Type> t_total(y_re.size());
  // t_total.setZero();
  // vector<Type> j_total(j_re.size());
  // j_total.setZero();
  // vector<Type> l_total(l_re.size());
  // l_total.setZero();
  // vector<Type> j_est(j_re.size());
  // j_est.setZero();
  // vector<Type> l_est(l_re.size());
  // l_est.setZero();
  // vector<Type> y_est(y_re.size());
  // y_est.setZero();
  // array<Type> proj_ef(n_i,4,9);
  
  for(int i=0; i<n_i; i++){
    eta_i(i) =  mu +
      bypass * a_i(i) +
      y_re(t_i(i),a_i(i)) +
      l_re(l_i(i),a_i(i)) +
      j_re(j_i(i),a_i(i)) +
      ln_zexp_sp_jl(x_s_jl(s_i_jl(i))) +  //day X length
      // ln_zexp_sp_jlt(x_s_jl(s_i_jl(i)),t_i(i)); //day X length X year
      z_jlt(x_s_jl(s_i_jl(i)), a_i(i), t_i(i)); //day X length X year
      
    // for(int jj = 0; jj<9;jj++){
      // tmp(i,jj) = mu +
      //   bypass * a_i(i) +
      //   y_re(t_i(i)) + //Marginal y
      //   j_re(j_shift(i,jj)) + //marginal day
      //   l_re(l_shift(i,jj),a_i(i)) + //marginal length
      //   ln_zexp_sp_jl(x_s_jl(management(i,jj))) +  //day X length
      //   ln_zexp_sp_jlt(x_s_jl(management(i,jj)),t_i(i)); //day X length X year
      // 
      // proj_ef(i,0,jj) = y_re(t_i(i));
      // proj_ef(i,1,jj) = j_re(j_shift(i,jj));
      // proj_ef(i,2,jj) = l_re(l_shift(i,jj),a_i(i));
      // proj_ef(i,3,jj) = ln_zexp_sp_jlt(x_s_jl(management(i,jj)),t_i(i));
      // 
      // tmp2(i,jj) = mu +
      //   y_re(t_i(i)) + //Marginal y
      //   j_re(j_shift(i,jj)) + //marginal day
      //   l_re(l_shift(i,jj),a_i(i)); //marginal length
      // 
      // proj_surv(t_i(i),jj) += invlogit(tmp(i,jj)) * total(i);
      // proj_surv_nox(t_i(i),jj) += invlogit(tmp2(i,jj)) * total(i);
      // proj_surv_ag(jj) += invlogit(tmp(i,jj)) * total(i);
    // }

    // t_total(t_i(i)) += total(i);
    // j_total(j_i(i)) += total(i);
    // l_total(l_i(i)) += total(i);
    nu_i(i) = exp(eta_i(i))/(1+exp(eta_i(i)));
    // j_est(j_i(i)) += invlogit(eta_i(i))*total(i);
    // l_est(l_i(i)) += invlogit(eta_i(i))*total(i);
    // y_est(t_i(i)) += invlogit(eta_i(i))*total(i);
    
    nll(6) -= dbinom(surv(i),total(i),nu_i(i), true );
  }

  // if(sim){
  //   SIMULATE{
  //     for(int i=0; i<n_i; i++){
  //       surv(i) = rbinom(total(i),nu_i(i));
  //     }
  //     REPORT(surv);
  //   }
  // }
  

  // for(int i=0;i<j_re.size();i++){
  //   j_est(i) /= j_total(i);
  // }
  // for(int i=0;i<l_re.size();i++){
  //   l_est(i) /= l_total(i);
  // }
  // for(int i=0;i<y_re.size();i++){
  //   y_est(i) /= t_total(i);
  // }
  
  
  // for(int jj = 0; jj<9;jj++){
  //   proj_surv_ag(jj) = log(proj_surv_ag(jj)/total.sum());
  //   for(int i=0; i<y_re.size(); i++){
  //     proj_surv(i,jj) =  log(proj_surv(i,jj)/t_total(i));
  //     proj_surv_nox(i,jj) =  log(proj_surv_nox(i,jj)/t_total(i));
  //     
  //   }
  // }

  REPORT(mu);
  REPORT(n_t);
  REPORT(z_jl);
  REPORT(nu_i);
  REPORT(phi_y);
  REPORT(phi_l);
  REPORT(phi_j);
  // REPORT(ln_j);
  // REPORT(ln_y);
  // REPORT(ln_l);
  // REPORT(l_seq);
  // REPORT(j_seq);
  // REPORT(beta);
  // REPORT(s);
  REPORT(ln_zexp_sp_jl);
  REPORT(ln_zexp_sp_jlt);
  REPORT(l_re);
  REPORT(y_re);
  REPORT(j_re);
  REPORT(z_jlt);
  REPORT(eta_i);
  // REPORT(tmp);
  // REPORT(tmp2);
  // REPORT(proj_ef);
  // REPORT(proj_surv_ag);
  // REPORT(l_est);
  // REPORT(j_est);
  REPORT(sig_l);
  REPORT(sig_j);
  REPORT(sig_y);

  
  Type Range_raw_jl = sqrt(8.0) / exp( log_kappa_jl );

  REPORT(Range_raw_jl);
  
  // vector<Type> mar_l(l_re.size());
  // for(int i=0; i<l_re.size();i++){
  //   mar_l(i) = (mu + l_re(i));
  // }
  // vector<Type> mar_j(j_re.size());
  // for(int i=0; i<j_re.size();i++){
  //   mar_j(i) = (mu + j_re(i));
  // }
  // 
  // vector<Type> mar_y(y_re.size());
  // for(int i=0; i<y_re.size();i++){
  //   mar_y(i) = (mu + y_re(i));
  // }

  // REPORT(mar_l);
  // REPORT(mar_j);
  // REPORT(mar_y);
  // REPORT(proj_surv);
  // REPORT(proj_surv_nox);
  // 
  // ADREPORT(l_est);
  // ADREPORT(j_est);
  // ADREPORT(y_est);
  
  // REPORT(j_total);
  // REPORT(l_total);
  // REPORT(t_total);
  
  // if(l_flag){
    // ADREPORT(mar_l);
    ADREPORT(l_re);
  // }
  // if(j_flag){
    // ADREPORT(mar_j);
    ADREPORT(j_re);
  // }
  // if(t_flag){
    // ADREPORT(mar_y);
    ADREPORT(y_re);
  // }
  // if(ad_proj){
    // ADREPORT(proj_surv);
    // ADREPORT(proj_surv_ag);
  // }

  REPORT(bypass);
  REPORT(nll);
  ADREPORT(mu);
  ADREPORT(bypass);
  
  return nll.sum();
}
