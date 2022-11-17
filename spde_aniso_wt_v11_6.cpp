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
  DATA_INTEGER(jl_bypass);
  DATA_INTEGER(jlt_bypass);

  DATA_INTEGER(t_AR);
  DATA_INTEGER(j_AR);
  DATA_INTEGER(l_AR);
  
  DATA_INTEGER(jlt_flag);
  DATA_INTEGER(jl_flag);
  
  DATA_INTEGER( n_i );         // Total number of observations
  // DATA_INTEGER( n_s_jl );        //Number of vertices
  DATA_IVECTOR( x_s_jl );	      // Association of each station with a given vertex in SPDE mesh
  DATA_IVECTOR( s_i_jl );	      // Association of each station with a given vertex in SPDE projection mesh
  DATA_IVECTOR( t_i );	      // year
  DATA_IVECTOR( l_i );	      // length
  DATA_IVECTOR( j_i );	      // day
  DATA_IVECTOR( a_i );	      // bypass
  DATA_IARRAY(m);   //management projections
  DATA_IVECTOR(m_imv)
  DATA_VECTOR(total); //total number of fish of observed
  DATA_VECTOR(surv); //data number of fish that survived
  DATA_STRUCT(spde_jl,spde_aniso_t);
  DATA_SCALAR(sim_size);
  DATA_INTEGER(proj_sim);
  DATA_VECTOR(proj_H);
  
  
  PARAMETER(mu);
  PARAMETER(bypass);
  PARAMETER(log_kappa_jl);
  PARAMETER(log_tau_jl);
  PARAMETER(log_tau_jl2);
  PARAMETER_VECTOR(ln_H_input_jl);
  PARAMETER_ARRAY(z_jlt);
  PARAMETER_ARRAY(z_jl);
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
  PARAMETER(f_psijlt);
  
  
  Type tau_jl2 = exp(log_tau_jl2);
  Type tau_jl = exp(log_tau_jl);
  Type kappa_jl = exp(log_kappa_jl);
  
  vector<Type> nll(10);
  nll.setZero();
  // Need to parameterize H matrix such that det(H)=1 (preserving volume) 
  // Note that H appears in (20) in Lindgren et al 2011
  matrix<Type> H_jl(2,2);
  if((jlt_flag == 1) | (jl_flag == 1)){
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
  Type psi_jlt = 1/(1+exp(f_psijlt));
  
  
  //2-D smoother
  SparseMatrix<Type> Q_jl = Q_spde(spde_jl,kappa_jl,H_jl);

  int n_t = 22;  
  //Scale the 2D variance  
  matrix<Type> jl_cov(2,2);
  jl_cov(0,0) = 1.;
  jl_cov(1,1) = 1.;
  jl_cov(0,1) = psi_jl;
  jl_cov(1,0) = psi_jl;
  MVNORM_t<Type> jl_dnorm(jl_cov); 
  if(jl_flag==1){
    if(jl_bypass == 0){
      nll(1) += SCALE(GMRF(Q_jl), 1/tau_jl)(z_jl.col(0));
    }
    if(jl_bypass == 1){
      nll(2) += SCALE(SEPARABLE(jl_dnorm, GMRF(Q_jl)), 1/tau_jl)( z_jl );
    }
  }
  
  matrix<Type> jlt_cov(2,2);
  jlt_cov(0,0) = 1.;
  jlt_cov(1,1) = 1.;
  jlt_cov(0,1) = psi_jlt;
  jlt_cov(1,0) = psi_jlt;
  matrix<Type> z_sim(z_jlt.dim[0], z_jlt.dim[1] );
  MVNORM_t<Type> jlt_dnorm(jlt_cov); 
  if(jlt_flag==1){
    for(int t=0; t<n_t; t++){
      if(jlt_bypass == 0){
        nll(2) += SCALE(GMRF(Q_jl), 1/tau_jl2)( z_jlt.col(t));
        if(proj_sim){
          SIMULATE{
            z_jlt.col(t) = GMRF(Q_jl).simulate();
            z_jlt.col(t) /= tau_jl2;
          } 
        }
      }
      if(jlt_bypass == 1){
        nll(2) += SCALE(SEPARABLE(jlt_dnorm, GMRF(Q_jl)), 1/tau_jl2)( z_jlt.col(t));
        if(proj_sim){
          SIMULATE{
            // z_jlt.col(t) =SEPARABLE(jl_dnorm, GMRF(Q_jl)).simulate();
            // z_jlt.col(t) /= tau_jl2;
          } 
        }
      }
    }
  }


  
  matrix<Type> y_cov(2,2);
  y_cov(0,0) = sig_y * sig_y;
  y_cov(1,1) = sig_y * sig_y;
  y_cov(0,1) = psi_y * sig_y * sig_y;
  y_cov(1,0) = psi_y * sig_y * sig_y;
  MVNORM_t<Type> y_dnorm(y_cov); 
  if(t_flag==1){
    if(t_AR == 2){
      if(t_bypass == 1){
        nll(5) += y_dnorm(y_re.transpose().col(0));
        for(int i = 1; i < y_re.dim(0); i++){
          nll(5) += y_dnorm(y_re.transpose().col(i) - y_re.transpose().col(i-1));
        }
      }
      if(t_bypass == 0){
        nll(5) -= dnorm(y_re(0,0),Type(0.),Type(10.),true);
        for(int i = 1; i < y_re.dim(0); i++){
          nll(5) -= dnorm(y_re(i,0),y_re(i-1,0),sig_y,true);
        }
      }
    }
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
  MVNORM_t<Type> l_dnorm(l_cov); 
  if(l_flag==1){
    if(l_AR == 2){
      if(l_bypass == 1){
        nll(4) += l_dnorm(l_re.transpose().col(0));
        for(int i = 1; i < l_re.dim(0); i++){
          nll(4) += l_dnorm(l_re.transpose().col(i) - l_re.transpose().col(i-1));
        }
      }
      if(l_bypass == 0){
        nll(4) -= dnorm(l_re(0,0),Type(0.),Type(10.),true);
        for(int i = 1; i < l_re.dim(0); i++){
          nll(4) -= dnorm(l_re(i,0),l_re(i-1,0),sig_l,true);
        }
      }
    }
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
          nll(4) -= dnorm(l_re(i,0),Type(0.),sig_l,true);
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
    if(j_AR == 2){
      if(j_bypass == 1){
        nll(5) += j_dnorm(j_re.transpose().col(0));
        for(int i = 1; i < j_re.dim(0); i++){
          nll(5) += j_dnorm(j_re.transpose().col(i) - j_re.transpose().col(i-1));
        }
      }
      if(j_bypass == 0){
        nll(5) -= dnorm(j_re(0,0),Type(0.),Type(10.),true);
        for(int i = 1; i < j_re.dim(0); i++){
          nll(5) -= dnorm(j_re(i,0),j_re(i-1,0),sig_j,true);
        }
      }
    }
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
  matrix<Type> proj(3,m_imv(1));
  proj.setZero();
  array<Type> proj_y(3,m_imv(1),n_t);
  proj_y.setZero();
  
  vector<Type> a_total(2);
  a_total = 0;
  matrix<Type> ay_total(2,n_t);
  ay_total.setZero();
  
  for(int i = 0; i < n_i; i++){

    eta_i(i) =  mu +
      a_i(i) * bypass +
      y_re(t_i(i) , a_i(i) * t_bypass) +
      l_re(l_i(i) , a_i(i) * l_bypass) +
      j_re(j_i(i) , a_i(i) * j_bypass) +
      z_jl(s_i_jl(i), a_i(i) * jl_bypass) + //day X length X year
      z_jlt(s_i_jl(i), a_i(i) * jlt_bypass, t_i(i)); //day X length X year
    
    nu_i(i) = exp(eta_i(i))/(1+exp(eta_i(i)));

    //Do the projections
    for(int mm = 0; mm < m_imv(1); mm++){
      
      int ss = m(i,mm,2);
      int jj = m(i,mm,1);
      int ll = m(i,mm,0);
      
      Type wt = total(i) * invlogit( mu +
        a_i(i) * bypass +
        y_re(t_i(i) , a_i(i) * t_bypass) +
        l_re(ll , a_i(i) * l_bypass) +
        j_re(jj , a_i(i) * j_bypass) +
        z_jl(ss, a_i(i) * jl_bypass) + //day X length X year
        z_jlt(ss, a_i(i) * jlt_bypass, t_i(i))); //day X length X year
      
      proj(a_i(i),mm) +=  wt;     // 
      proj_y(a_i(i),mm,t_i(i)) +=  wt;     // 
        // 
      proj(2,mm) += (eta_i(i)) * total(i);
    }
    a_total(a_i(i)) += total(i);
    ay_total(a_i(i),t_i(i)) += total(i);
    
    nll(6) -= dbinom(surv(i),total(i),nu_i(i), true );
  }
  
  SIMULATE{
    for(int i = 0; i < n_i; i++){
      Type tmp_total = total(i) * sim_size;
      surv(i) = rbinom( total(i), nu_i(i) );
    }
    REPORT(surv);
  }
  
  
  for(int mm = 0; mm < m_imv(1); mm++){
    for(int aa = 0; aa < 2; aa++){
      proj(aa,mm) /= a_total(aa);
      for(int y = 0; y < n_t; y++){
        proj_y(aa,mm,y) /= ay_total(aa,y);
      }
    }
    proj(2,mm) /= total.sum();
  }
  
  for(int mm = 0; mm < m_imv(1); mm++){
    for(int aa = 0; aa < 2; aa++){
      for(int y = 0; y < n_t; y++){
        if(mm != 4){
          proj_y(aa,mm,y) = (proj_y(aa,mm,y) - proj_y(aa,4,y))/proj_y(aa,4,y) * 100;
        }
      }
    }
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
  REPORT(proj);
  REPORT(a_total);
  REPORT(ay_total);
  REPORT(proj_y);

  
  Type Range_raw_jl = sqrt(8.0) / exp( log_kappa_jl );

  REPORT(Range_raw_jl);
  array<Type> mar_j(j_re.dim(0),3); 
  vector<int> j_int(3);
  j_int << 14,40,67;
  for(int i = 0; i < 3; i++ ){
    mar_j.col(i) = j_re.col(0) + j_re(j_int(i)) + mu;
  }

  array<Type> mar_l(l_re.dim(0),3); 
  vector<int> l_int(3);
  l_int << 38,25,5;
  for(int i = 0; i < 3; i++ ){
    mar_l.col(i) = l_re.col(0) + l_re(l_int(i)) + mu;
  }

  array<Type> mar_y(y_re.dim(0),3); 
  for(int i = 0; i < 3; i++ ){
    mar_y.col(i) = y_re.col(0) + 
    l_re(l_int(i)) +
    j_re(j_int(i)) +
    mu;
  }
  REPORT(mar_l);
  REPORT(mar_j);
  REPORT(mar_y);
  ADREPORT(mar_l);
  ADREPORT(mar_j);
  ADREPORT(mar_y);
  ADREPORT(mu);
  ADREPORT(mu+bypass);
  ADREPORT(proj);
  ADREPORT(proj_y);
  
  return nll.sum();
}
