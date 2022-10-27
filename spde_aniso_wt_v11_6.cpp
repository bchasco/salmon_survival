
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
  
  DATA_INTEGER(t_bypass);
  DATA_INTEGER(j_bypass);
  DATA_INTEGER(l_bypass);
  DATA_INTEGER(jlt_bypass);

  DATA_INTEGER(t_AR);
  DATA_INTEGER(j_AR);
  DATA_INTEGER(l_AR);
  
  DATA_INTEGER(jlt_flag);
  DATA_INTEGER( n_i );         // Total number of observations
  // DATA_INTEGER( n_s_jl );        //Number of vertices
  DATA_IVECTOR( x_s_jl );	      // Association of each station with a given vertex in SPDE mesh
  DATA_IVECTOR( s_i_jl );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IVECTOR( t_i );	      // year
  DATA_IVECTOR( l_i );	      // length
  DATA_IVECTOR( j_i );	      // day
  DATA_IVECTOR( a_i );	      // bypass

  DATA_VECTOR(total); //total number of fish of observed
  DATA_VECTOR(surv); //data number of fish that survived
  DATA_STRUCT(spde_jl,spde_aniso_t);

  
  PARAMETER(mu);
  PARAMETER(bypass);
  PARAMETER(log_kappa_jl);
  PARAMETER(log_tau_jl2);
  PARAMETER_VECTOR(ln_H_input_jl);
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
  
  
  Type tau_jl2 = exp(log_tau_jl2);
  Type kappa_jl = exp(log_kappa_jl);
  
  vector<Type> nll(10);
  nll.setZero();
  // Need to parameterize H matrix such that det(H)=1 (preserving volume) 
  // Note that H appears in (20) in Lindgren et al 2011
  matrix<Type> H_jl(2,2);
  if(jlt_flag==1){
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

  int n_t = 22;  
  //Scale the 2D variance  
  // matrix<Type> ln_zexp_sp_jlt( n_s_jl, n_t);
  matrix<Type> jl_cov(2,2);
  jl_cov(0,0) = 1.;
  jl_cov(1,1) = 1.;
  jl_cov(0,1) = psi_jl;
  jl_cov(1,0) = psi_jl;
  MVNORM_t<Type> jl_dnorm(jl_cov); 
  if(jlt_flag==1){
    for(int t=0; t<n_t; t++){
      if(jlt_bypass == 0){
        nll(2) += SCALE(GMRF(Q_jl), 1/tau_jl2)( z_jlt.col(t));
      }
      if(jlt_bypass == 1){
        nll(2) += SCALE(SEPARABLE(jl_dnorm, GMRF(Q_jl)), 1/tau_jl2)( z_jlt.col(t));
      }
    }
  }
  
  
  matrix<Type> y_cov(2,2);
  y_cov(0,0) = sig_y * sig_y;
  y_cov(1,1) = sig_y * sig_y;
  y_cov(0,1) = psi_y * sig_y * sig_y;
  y_cov(1,0) = psi_y * sig_y * sig_y;
  if(t_flag==1){
    if(t_AR == 1){
      if(t_bypass == 1){
        nll(3) += AR1(phi_y,MVNORM(y_cov))(y_re.transpose());
      }
      if(t_bypass == 0){
        nll(3) += SCALE(AR1(phi_y),sig_y)(y_re.col(0));
      }
    }
    if(t_AR == 0){
      if(t_bypass == 1){
        for(int i = 0; i < n_t; i++){
          nll(3) += MVNORM(y_cov)(y_re.transpose().col(i));
        }
      }
      if(t_bypass == 0){
        for(int i = 0; i < n_t; i++){
          nll(3) -= dnorm(y_re(i,0),Type(0.),sig_y,true);
        }
      }
    }
    
  }
  
  matrix<Type> l_cov(2,2);
  l_cov(0,0) = sig_l * sig_l;
  l_cov(1,1) = sig_l * sig_l;
  l_cov(0,1) = psi_l * sig_l * sig_l;
  l_cov(1,0) = psi_l * sig_l * sig_l;
  if(l_flag==1){
    if(l_AR == 1){
      if(l_bypass == 1){
        nll(4) += AR1(phi_l,MVNORM(l_cov))(l_re.transpose());
      }
      if(l_bypass == 0){
        nll(4) += SCALE(AR1(phi_l),sig_l)(l_re.col(0));
      }
    }
    if(l_AR == 0){
      if(l_bypass == 1){
        for(int i = 0; i < l_re.dim(0); i++){
          nll(4) += MVNORM(l_cov)(l_re.transpose().col(i));
        }
      }
      if(l_bypass == 0){
        for(int i = 0; i < l_re.dim(0); i++){
          nll(4) -= dnorm(l_re(i,0),Type(0.),sig_y,true);
        }
      }
    }
  }
  
  
  matrix<Type> j_cov(2,2);
  j_cov(0,0) = sig_j * sig_j;
  j_cov(1,1) = sig_j * sig_j;
  j_cov(0,1) = psi_j * sig_j * sig_j;
  j_cov(1,0) = psi_j * sig_j * sig_j;
  MVNORM_t<Type> j_dnorm(j_cov); 
  if(j_flag==1){
    if(j_AR == 1){
      if(j_bypass == 1){
        nll(5) += AR1(phi_j,j_dnorm)(j_re.transpose());
      }
      if(j_bypass == 0){
        nll(5) += SCALE(AR1(phi_j),sig_j)(j_re.col(0));
      }
    }
    if(j_AR == 0){
      if(j_bypass == 1){
        for(int i = 0; i < j_re.dim(0); i++){
          nll(5) += MVNORM(j_cov)(j_re.transpose().col(i));
        }
      }
      if(j_bypass == 0){
        for(int i = 0; i < j_re.dim(0); i++){
          nll(5) -= dnorm(j_re(i,0),Type(0.),sig_j,true);
        }
      }
    }
  }
  
  vector<Type> nu_i(n_i);
  vector<Type> eta_i(n_i);
  for(int i=0; i<n_i; i++){

    eta_i(i) =  mu +
      a_i(i) * bypass +
      y_re(t_i(i) , a_i(i) * t_bypass) +
      l_re(l_i(i) , a_i(i) * l_bypass) +
      j_re(j_i(i) , a_i(i) * j_bypass) +
      // z_jlt(x_s_jl(s_i_jl(i)), a_i(i) * jlt_bypass, t_i(i)); //day X length X year
    z_jlt(s_i_jl(i), a_i(i) * jlt_bypass, t_i(i)); //day X length X year
    
    nu_i(i) = exp(eta_i(i))/(1+exp(eta_i(i)));

    nll(6) -= dbinom(surv(i),total(i),nu_i(i), true );
  }


  REPORT(mu);
  REPORT(bypass);
  REPORT(nll);
  REPORT(n_t);
  REPORT(nu_i);
  REPORT(phi_y);
  REPORT(phi_l);
  REPORT(phi_j);
  REPORT(l_re);
  REPORT(y_re);
  REPORT(j_re);
  REPORT(z_jlt);
  REPORT(eta_i);
  REPORT(sig_l);
  REPORT(sig_j);
  REPORT(sig_y);

  
  Type Range_raw_jl = sqrt(8.0) / exp( log_kappa_jl );

  REPORT(Range_raw_jl);
  array<Type> mar_l(l_re.dim(0),l_re.dim(1)); mar_l = l_re + j_re(40) + mu;
  array<Type> mar_j(j_re.dim(0),j_re.dim(1)); mar_j = j_re + l_re(22) + mu;
  array<Type> mar_y(y_re.dim(0),y_re.dim(1)); mar_y = y_re + l_re(22) + j_re(40) + mu;
  ADREPORT(mar_l);
  ADREPORT(mar_j);
  ADREPORT(mar_y);
  ADREPORT(mu);
  ADREPORT(bypass);
  
  return nll.sum();
}
