// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include "rcpp_pgdraw.h"
using namespace Rcpp;

// generate a number from the inverse gaussian distribution
double rinvgauss(double mu, double lambda) {
  if(arma::is_finite(mu)) {
    double z,y,x,u;
    z=R::rnorm(0,1);
    y=z*z;
      x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4.0*mu*lambda*y+mu*mu*y*y);
    u=R::runif(0,1);
    if(u <= mu/(mu+x)){
      return x;
    }else{
      return mu*mu/x;
    }
  } else{
      return 1.0/pow(R::rnorm(0,1),2.0)*lambda;
}
}

// MCMC algorihtm
// y: response, vector of length n
// X: predictors, n by p matrix
// A_tilde: the existing knowledge graph, p by p matrix
// A_star: the graph estimated by graphical lasso, p by p matrix
// mu, nu : parameters of the prior distribution of theta, double precision
// a_sigma, b_sigma: parameters of the prior distribution of sigma^2, double precision
// a_xi, b_xi: parameters of the prior distribution of X_i, double precision
// a_q, b_q: parameters of the prior distribution of the true graph, double precision
// u00_mean_prior, sigma_u00: parameters of the prior distribution of $\tilde s_0$, double precision
// u10_mean_prior, sigma_u10: parameters of the prior distribution of $\tilde s_1$, double precision
// u00_star_mean_prior, sigma_u00_star: parameters of the prior distribution of $s*_0$, double precision
// u10_star_mean_prior, sigma_u10_star: parameters of the prior distribution of $s*_1$, double precision
// r: length of latent scales $\tilde u_{0j}$ and $\tilde u_{1j}$, integer
// r_star: length of latent scales $u*_{0j}$ and $u*_{1j}$, integer
// uj1: the first element of every latent scale, double precision
// u0_mean_prior: mean of the prior distribution of (\tilde u_{0j}[1:r-1],u*_{0j}[1:r_star-1]), vector of length r+r_star-2
// u1_mean_prior: mean of the prior distribution of (\tilde u_{1j}[1:r-1],u*_{1j}[1:r_star-1]), vector of length r+r_star-2
// Sigma0: variance of the prior distribution of (\tilde u_{0j}[1:r-1],u*_{0j}[1:r_star-1]), r+r_star-2 by r+r_star-2 matrix
// Sigma1: variance of the prior distribution of (\tilde u_{1j}[1:r-1],u*_{1j}[1:r_star-1]), r+r_star-2 by r+r_star-2 matrix
// n_iter: number of iterations, integer
// k: interval length for computing the acceptance rate of Metropolis/MH algorithm, integer. (set k < n_iter if we want to update proposal variance)

// [[Rcpp::export]]
List MCMC(arma::vec y, arma::mat X, arma::mat A_tilde,arma::mat A_star, int n_iter,double mu,double nu, double a_sigma,double b_sigma,double a_xi,double b_xi,double a_q, double b_q, 
                double u00_mean_prior,double sigma_u00, double u10_mean_prior, double sigma_u10,double u00_star_mean_prior,double sigma_u00_star, double u10_star_mean_prior, double sigma_u10_star, 
                int r, int r_star, double uj1,arma::vec u0_mean_prior,arma::vec u1_mean_prior, arma::mat Sigma0,arma::mat Sigma1, int k) {
  int p = X.n_cols, n = X.n_rows, r_total= r+r_star;
  arma::vec beta_sd_prop(p,arma::fill::ones), theta_sd_prop(p,arma::fill::ones);
  double beta_mean_prior = 0.0, theta_mean_prior= mu;
  
 arma::mat BETA( p,n_iter, arma::fill::zeros), THETA( p, n_iter,arma::fill::none);
 arma::vec SIGMA2(n_iter, arma::fill::ones);

  arma::vec beta = BETA.col(0), theta(p,arma::fill::none);
  theta.fill(mu);
  arma::mat a =  A_tilde;
    double sigma2 = 1.0, q = a_q/(a_q+b_q);
  
  arma::mat acs_theta(p,n_iter,arma::fill::zeros), acs_beta(p,n_iter,arma::fill::zeros);
  
  arma::mat r_theta(p,floor(n_iter/k),arma::fill::zeros), r_beta(p,floor(n_iter/k),arma::fill::zeros);
  double r_opt = 0.44, c1 = 0.8;
  
  arma::vec adaptFactor_vec(p,arma::fill::none);
  double lhr, sigma2_prop, gamma1, gamma2, adaptFactor;
  int s;
  
  arma::vec theta_tmp;
  
  arma::mat xx=trans(X)*X, yx=y.t()*X, mean_tmp(n,1,arma::fill::none), theta_diff(p,p,arma::fill::none),
    prob_ratio(p,p,arma::fill::none), prob_A(p,p,arma::fill::none);
  double betaj, thetaj, beta_diff, beta_prop, theta_prop, A_tilde_sum=arma::accu(A_tilde)/2.0,a_sum,
    a_q_pos,b_q_pos, A_tilde_sum0, a_alpha0_pos, b_alpha0_pos, A_tilde_sum1, a_alpha1_pos,
    b_alpha1_pos, prob_tmp;
  arma::uvec upper_index = trimatu_ind(size(a),1);
  
  double gamma_a_xi = std::tgamma(a_xi);
  double value_tmp,theta_j2,A_tilde_tmp;
    
    arma::uvec indices_r(r-1), indices_r_star(r_star-1);
    
    for(int l=0;l<(r-1);++l){
      indices_r(l) = l+1;
    }
    for(int l=0;l<(r_star-1);++l){
      indices_r_star(l) = l+1;
    }
 
   double u00 = u00_mean_prior, u10 = u10_mean_prior,u00_star = u00_star_mean_prior, u10_star = u10_star_mean_prior;
   arma::mat u0(p,r,arma::fill::randn), u1(p,r,arma::fill::randn), u0_star(p,r_star,arma::fill::randn), u1_star(p,r_star,arma::fill::randn);
   arma::vec u0_mean_prior_tmp(r),u1_mean_prior_tmp(r),u0_star_mean_prior_tmp(r_star),u1_star_mean_prior_tmp(r_star);
   u0_mean_prior_tmp.fill(uj1),u1_mean_prior_tmp.fill(uj1),u0_star_mean_prior_tmp.fill(uj1),u1_star_mean_prior_tmp.fill(uj1);
   u0_mean_prior_tmp(indices_r) = u0_mean_prior.subvec(0,r-2);
   u1_mean_prior_tmp(indices_r) = u1_mean_prior.subvec(0,r-2);
   u0_star_mean_prior_tmp(indices_r_star) = u0_mean_prior.subvec(r-1,r_total-3);
   u1_star_mean_prior_tmp(indices_r_star) = u1_mean_prior.subvec(r-1,r_total-3);
   
   u0.each_row()= arma::conv_to< arma::rowvec >::from(u0_mean_prior_tmp);
   u1.each_row()= arma::conv_to< arma::rowvec >::from(u1_mean_prior_tmp);
   u0_star.each_row()= arma::conv_to< arma::rowvec >::from(u0_star_mean_prior_tmp);
   u1_star.each_row()= arma::conv_to< arma::rowvec >::from(u1_star_mean_prior_tmp);
    
    
  arma::mat uu1, uu0,uu1_star, uu0_star,adaptFactor_mat;
  arma::uvec a_zero, a_nonzero, index_tmp;
  arma::vec  uu0_part, uu0_part_prop, A_tilde_part, uu1_part, uu1_part_prop, uu0_star_part, uu0_star_part_prop, A_star_part, uu1_star_part, uu1_star_part_prop;
  arma::rowvec a_tildej,u0j, a_tildej_long, u0j_prop, u0_part, u0_part_prop,u1j, u1j_prop, u1_part, u1_part_prop, u0j_star,u0j_star_prop, u0_star_part, u0_star_part_prop,u1j_star, u1j_star_prop, u1_star_part, u1_star_part_prop;
  
  double lhr0,lhr1,lhr_org,lhr_prop,lhr_org1,lhr_prop1, lhr_star,lhr0_star,lhr1_star,lhr_star_org,lhr_star_prop,lhr_star_org1,lhr_star_prop1;
  arma::mat u0_rows,u1_rows,u0_star_rows,u1_star_rows;
    
    double FNrate=0.0,FPrate=0.0,TP=0.0,Accu=0.0, A0_sum1,
      Alpha0Est=0.0,Alpha1Est=0.0;
  arma::vec beta_postmed(p,arma::fill::none);
  //
  arma::vec tau2_inv(p,arma::fill::none), xy=trans(X)*y;
  tau2_inv = exp(2.0*theta)/2.0;
  arma::mat B_inv(p,p,arma::fill::none),D_tau2_inv(p,p,arma::fill::zeros), x_t=trans(X),W0_tilde(p,p,arma::fill::ones),W1_tilde(p,p,arma::fill::ones),W0_star(p,p,arma::fill::ones),W1_star(p,p,arma::fill::ones),
    Kappa_tilde = A_tilde-0.5, Kappa_star=A_star-0.5, Sigma0_inv=arma::inv_sympd(Sigma0), Sigma1_inv=arma::inv_sympd(Sigma1);
  D_tau2_inv = arma::diagmat(tau2_inv);
    double tau2_inv_u00_tilde, tau2_inv_u10_tilde,tau2_inv_u00_star, tau2_inv_u10_star, u00_tilde_mean_post_top, u10_tilde_mean_post_top, u00_star_mean_post_top, u10_star_mean_post_top, u00_tilde_mean_post_bottom, u10_tilde_mean_post_bottom, u00_star_mean_post_bottom, u10_star_mean_post_bottom;
  tau2_inv_u00_tilde = 1.0/2.0/pow(sigma_u00,2.0), tau2_inv_u10_tilde = 1.0/2.0/pow(sigma_u10,2.0), tau2_inv_u00_star = 1.0/2.0/pow(sigma_u00_star,2.0), tau2_inv_u10_star = 1.0/2.0/pow(sigma_u10_star,2.0);
     
    arma::mat U0j, U0j_tilde,U0j_star, U1j, U1j_tilde,U1j_star, W0j, W0j_tilde,W0j_star, W1j, W1j_tilde,W1j_star,V1j(r_total,r_total),V1j_tilde(r-1,r-1),V1j_star(r_star-1,r_star-1), V0j(r_total,r_total),V0j_tilde(r-1,r-1),V0j_star(r_star-1,r_star-1),V1j_inv(r_total,r_total), V1j_tilde_inv(r-1,r-1), V1j_star_inv(r_star-1,r_star-1), V0j_inv(r_total,r_total),V0j_tilde_inv(r-1,r-1),V0j_star_inv(r_star-1,r_star-1),u0j_total(r_total,1), u1j_total(r_total,1);
    arma::vec Z0j, Z1j,w0j,w1j,kappa0j,kappa1j,Z0j_tilde, Z1j_tilde,w0j_tilde,w1j_tilde,kappa0j_tilde,kappa1j_tilde,Z0j_star, Z1j_star,w0j_star,w1j_star,kappa0j_star,kappa1j_star;
    int m;
    arma::mat Sigma0_11=Sigma0.submat(0,0,r-2,r-2), Sigma0_12=Sigma0.submat(0,r-1,r-2,r_total-3), Sigma0_21=Sigma0.submat(r-1,0,r_total-3,r-2),Sigma0_22=Sigma0.submat(r-1,r-1,r_total-3,r_total-3),Sigma1_11=Sigma1.submat(0,0,r-2,r-2), Sigma1_12=Sigma1.submat(0,r-1,r-2,r_total-3), Sigma1_21=Sigma1.submat(r-1,0,r_total-3,r-2),Sigma1_22=Sigma1.submat(r-1,r-1,r_total-3,r_total-3);
    arma::vec u0_mean_prior_1= u0_mean_prior.subvec(0,r-2), u0_mean_prior_2= u0_mean_prior.subvec(r-1,r_total-3),u1_mean_prior_1= u1_mean_prior.subvec(0,r-2), u1_mean_prior_2= u1_mean_prior.subvec(r-1,r_total-3);
    arma::mat Sigma0_part1 = Sigma0_12*arma::inv_sympd(Sigma0_22), Sigma0_part2 = Sigma0_21*arma::inv_sympd(Sigma0_11),Sigma1_part1 = Sigma1_12*arma::inv_sympd(Sigma1_22), Sigma1_part2 = Sigma1_21*arma::inv_sympd(Sigma1_11), Sigma0_tilde = Sigma0_11-Sigma0_part1*Sigma0_21, Sigma0_tilde_inv=arma::inv_sympd(Sigma0_tilde), Sigma1_tilde=Sigma1_11-Sigma1_part1*Sigma1_21,Sigma1_tilde_inv=arma::inv_sympd(Sigma1_tilde), Sigma0_star=Sigma0_22-Sigma0_part2*Sigma0_12, Sigma0_star_inv=arma::inv_sympd(Sigma0_star), Sigma1_star=Sigma1_22-Sigma1_part2*Sigma1_12, Sigma1_star_inv=arma::inv_sympd(Sigma1_star);
    double A_total=0.0, A_impimp=0.0, A_impnoimp=0.0,A_noimpnoimp=0.0,a_impimp, a_impnoimp,a_noimpnoimp;
    
  for(int i=0;i<n_iter; ++i) {

     // update beta
      B_inv =arma::inv(xx+D_tau2_inv);
      beta = B_inv*(xy+sqrt(sigma2)*(x_t*arma::randn<arma::vec>(n)+sqrt(tau2_inv)%arma::randn<arma::vec>(p)));
      BETA.col(i) = beta;

      // update tau
     for (int j=0; j<p; ++j) {
         tau2_inv[j] = rinvgauss(sqrt(sigma2)*exp(theta[j])/std::abs(beta[j]),exp(2.0*theta[j]) );
     }
     D_tau2_inv = arma::diagmat(tau2_inv);

     // update sigma2
     mean_tmp = y-X*beta;
     sigma2 = 1.0/R::rgamma( a_sigma+n/2.0+p/2.0, 1.0/(b_sigma+as_scalar(mean_tmp.t()*mean_tmp)/2.0+arma::sum(pow(beta,2)%tau2_inv)/2.0) );
     SIGMA2[i] = sigma2;

     // update theta
     for (int j=0; j<p; ++j) {
       thetaj = theta[j];
       theta_prop = R::rnorm(thetaj,theta_sd_prop[j]);
       theta_tmp = theta.elem(find(a.col(j)==1) );
       lhr = 2.0*theta_prop - 2.0*thetaj +
         (exp(2.0*thetaj)-exp(2.0*theta_prop))/2.0/tau2_inv[j]+
         - pow(theta_prop-theta_mean_prior,2.0)/2.0/nu + pow(thetaj-theta_mean_prior,2.0)/2.0/nu +
         -a_xi*arma::sum(log(b_xi+pow(theta_prop-theta_tmp,2.0)/2.0/nu ))-
         -a_xi* arma::sum(log(b_xi+pow(thetaj-theta_tmp,2.0)/2.0/nu ));
         if(lhr > log(R::runif(0,1)) ) { theta[j] = theta_prop; acs_theta(j,i)++;}
     }
    THETA.col(i) = theta;

    // update A, u00, u10, u00_star, u10_star, W
    uu1 = u10+u1*u1.t();
    uu0 = u00+u0*u0.t();
    uu1_star = u10_star+u1_star*u1_star.t();
    uu0_star = u00_star+u0_star*u0_star.t();
    a.fill(0.0);a_sum=0.0;A0_sum1=0.0; a_impimp=0.0; a_impnoimp=0.0; a_noimpnoimp=0.0;
    u00_tilde_mean_post_top = u00_mean_prior*tau2_inv_u00_tilde;
    u10_tilde_mean_post_top = u10_mean_prior*tau2_inv_u10_tilde;
    u00_star_mean_post_top = u00_star_mean_prior*tau2_inv_u00_star;
    u10_star_mean_post_top = u10_star_mean_prior*tau2_inv_u10_star;
    u00_tilde_mean_post_bottom = tau2_inv_u00_tilde;
    u10_tilde_mean_post_bottom = tau2_inv_u10_tilde;
    u00_star_mean_post_bottom = tau2_inv_u00_star;
    u10_star_mean_post_bottom = tau2_inv_u10_star;


    for(int j2=1; j2<p; ++j2){
      theta_j2 = theta(j2);
      for(int j1=0;j1<j2;++j1) {

          prob_tmp = q/(1.0-q)*gamma_a_xi/pow(b_xi+(pow(theta(j1)-theta_j2,2.0))/nu/2.0,a_xi)*exp(Kappa_tilde(j1,j2)*(uu1(j1,j2)-uu0(j1,j2))-W1_tilde(j1,j2)*pow(uu1(j1,j2),2.0)/2.0+ W0_tilde(j1,j2)*pow(uu0(j1,j2),2.0)/2.0+ Kappa_star(j1,j2)*(uu1_star(j1,j2)-uu0_star(j1,j2))-W1_star(j1,j2)*pow(uu1_star(j1,j2),2.0)/2.0+ W0_star(j1,j2)*pow(uu0_star(j1,j2),2.0)/2.0);
        prob_tmp /= (1.0+prob_tmp);
        value_tmp = R::rbinom(1,prob_tmp);
        if(value_tmp==1.0) {
          a(j1,j2) = value_tmp; a(j2,j1) = value_tmp;  a_sum++;
           W0_tilde(j1,j2) = samplepg(0.0); W0_tilde(j2,j1) = W0_tilde(j1,j2);
           W0_star(j1,j2) = samplepg(0.0); W0_star(j2,j1) = W0_star(j1,j2);
            W1_tilde(j1,j2) =samplepg(uu1(j1,j2)), W1_star(j1,j2)=samplepg(uu1_star(j1,j2)),W1_star(j2,j1) = W1_star(j1,j2),W1_star(j2,j1) = W1_star(j1,j2);
          u10_tilde_mean_post_top+= Kappa_tilde(j1,j2)-W1_tilde(j1,j2)*(uu1(j1,j2)-u10);
          u10_star_mean_post_top+= Kappa_star(j1,j2)-W1_star(j1,j2)*(uu1_star(j1,j2)-u10_star);
          u10_tilde_mean_post_bottom +=  W1_tilde(j1,j2);
          u10_star_mean_post_bottom +=  W1_star(j1,j2);
          if(A_tilde(j1,j2)==1.0) A0_sum1++;
        } else {
          W0_tilde(j1,j2) = samplepg(uu0(j1,j2)); W0_tilde(j2,j1) = W0_tilde(j1,j2);
          W0_star(j1,j2) = samplepg(uu0_star(j1,j2));W0_star(j2,j1) = W0_star(j1,j2);
           W1_tilde(j1,j2) =samplepg(0.0), W1_star(j1,j2)=samplepg(0.0),W1_star(j2,j1) = W1_star(j1,j2),W1_star(j2,j1) = W1_star(j1,j2);
          u00_tilde_mean_post_top+= Kappa_tilde(j1,j2)-W0_tilde(j1,j2)*(uu0(j1,j2)-u00);
          u00_star_mean_post_top+= Kappa_star(j1,j2)-W0_star(j1,j2)*(uu0_star(j1,j2)-u00_star);
          u00_tilde_mean_post_bottom+= W0_tilde(j1,j2);
          u00_star_mean_post_bottom+= W0_star(j1,j2);

        }
      }
    }
  

      u10 = R::rnorm(u10_tilde_mean_post_top/u10_tilde_mean_post_bottom,1.0/sqrt(u10_tilde_mean_post_bottom));
      u10_star = R::rnorm(u10_star_mean_post_top/u10_star_mean_post_bottom,1.0/sqrt(u10_star_mean_post_bottom));
      u00 = R::rnorm(u00_tilde_mean_post_top/u00_tilde_mean_post_bottom,1.0/sqrt(u00_tilde_mean_post_bottom));
      u00_star = R::rnorm(u00_star_mean_post_top/u00_star_mean_post_bottom,1.0/sqrt(u00_star_mean_post_bottom));



    // update q
    a_q_pos = a_q+ a_sum;
    b_q_pos = b_q+ p*(p-1.0)/2.0- a_sum;
    q = R::rbeta(a_q_pos,b_q_pos);

   // update tau2_inv_u00
    tau2_inv_u00_tilde = rinvgauss(1.0/sigma_u00/std::abs(u00-u00_mean_prior),1.0/pow(sigma_u00,2.0));
    tau2_inv_u10_tilde = rinvgauss(1.0/sigma_u10/std::abs(u10-u10_mean_prior),1.0/pow(sigma_u10,2.0));
    tau2_inv_u00_star = rinvgauss(1.0/sigma_u00_star/std::abs(u00_star-u00_star_mean_prior),1.0/pow(sigma_u00_star,2.0));
    tau2_inv_u10_star = rinvgauss(1.0/sigma_u10_star/std::abs(u10_star-u10_star_mean_prior),1.0/pow(sigma_u10_star,2.0));

    // update u0
      a.diag().ones();
      for (int j=0; j<p; ++j) {
          index_tmp = find(a.col(j)==0.0);
          if(index_tmp.n_elem > 0 ) {
              U0j_tilde = u0.submat(index_tmp,indices_r);
              U0j_tilde = U0j_tilde.tail_cols(r-1);
              w0j_tilde = W0_tilde.unsafe_col(j).elem(index_tmp);
              kappa0j_tilde = Kappa_tilde.unsafe_col(j).elem(index_tmp);
              W0j_tilde = arma::diagmat(w0j_tilde);
              Z0j_tilde = kappa0j_tilde/w0j_tilde+u00;
              V0j_tilde = U0j_tilde.t()*W0j_tilde*U0j_tilde+Sigma0_tilde_inv;
              V0j_tilde_inv=arma::inv_sympd(V0j_tilde);
              U0j_tilde = arma::mvnrnd(V0j_tilde_inv*(U0j_tilde.t()*W0j_tilde*Z0j_tilde+Sigma0_tilde_inv*(u0_mean_prior_1+Sigma0_part1*(arma::conv_to< arma::vec >::from(u0_star.submat(j,1,j,r_star-1)) -u0_mean_prior_2) ) ), V0j_tilde_inv);
              u0.submat(j,1,j,r-1) = arma::conv_to< arma::rowvec >::from(U0j_tilde);
              } else {
             u0.submat(j,1,j,r-1) = arma::conv_to< arma::rowvec >::from(arma::mvnrnd(u0_mean_prior_1+Sigma0_part1*(arma::conv_to< arma::vec >::from(u0_star.submat(j,1,j,r_star-1)) -u0_mean_prior_2) ,Sigma0_tilde));
          }
      }
    
    // update u0_star
      for (int j=0; j<p; ++j) {
          index_tmp = find(a.col(j)==0.0);
          if(index_tmp.n_elem > 0 ) {
              U0j_star = u0_star(index_tmp,indices_r_star);
              w0j_star = W0_star.unsafe_col(j).elem(index_tmp);
              kappa0j_star = Kappa_star.unsafe_col(j).elem(index_tmp);
              W0j_star = arma::diagmat(w0j_star);
              Z0j_star = kappa0j_star/w0j_star+u00_star;
              V0j_star = U0j_star.t()*W0j_star*U0j_star+Sigma0_star_inv;
              V0j_star_inv=arma::inv_sympd(V0j_star);
              U0j_star = arma::mvnrnd(V0j_star_inv*(U0j_star.t()*W0j_star*Z0j_star+Sigma0_star_inv*(u0_mean_prior_2+Sigma0_part2*(arma::conv_to< arma::vec >::from(u0.submat(j,1,j,r-1)) -u0_mean_prior_1) ) ), V0j_star_inv);
              u0_star.submat(j,1,j,r_star-1) = arma::conv_to< arma::rowvec >::from(U0j_star);
              } else {
             u0_star.submat(j,1,j,r_star-1) = arma::conv_to< arma::rowvec >::from(arma::mvnrnd(u0_mean_prior_2+Sigma0_part2*(arma::conv_to< arma::vec >::from(u0.submat(j,1,j,r-1)) -u0_mean_prior_1) ,Sigma0_star));
          }
      }
      a.diag().zeros();
      
       // update u1
      for (int j=0; j<p; ++j) {
        index_tmp = find(a.col(j)==1.0);
        if(index_tmp.n_elem > 0 ) {
          U1j_tilde = u1(index_tmp,indices_r);
          w1j_tilde = W1_tilde.unsafe_col(j).elem(index_tmp);
          kappa1j_tilde = Kappa_tilde.unsafe_col(j).elem(index_tmp);
          W1j_tilde = arma::diagmat(w1j_tilde);
          Z1j_tilde = kappa1j_tilde/w1j_tilde+u10;
          V1j_tilde = U1j_tilde.t()*W1j_tilde*U1j_tilde+Sigma1_tilde_inv;
          V1j_tilde_inv=arma::inv_sympd(V1j_tilde);
          U1j_tilde = arma::mvnrnd(V1j_tilde_inv*(U1j_tilde.t()*W1j_tilde*Z1j_tilde+Sigma1_tilde_inv*(u1_mean_prior_1+Sigma1_part1*(arma::conv_to< arma::vec >::from(u1_star.submat(j,1,j,r_star-1)) -u1_mean_prior_2) ) ), V1j_tilde_inv);
          u1.submat(j,1,j,r-1) = arma::conv_to< arma::rowvec >::from(U1j_tilde);
        } else {
          u1.submat(j,1,j,r-1) = arma::conv_to< arma::rowvec >::from(arma::mvnrnd(u1_mean_prior_1+Sigma1_part1*(arma::conv_to< arma::vec >::from(u1_star.submat(j,1,j,r_star-1)) -u1_mean_prior_2) ,Sigma1_tilde));
        }
      }
      
       // update u1_star
      for (int j=0; j<p; ++j) {
        index_tmp = find(a.col(j)==1.0);
        if(index_tmp.n_elem > 0 ) {
          U1j_star = u1_star(index_tmp,indices_r_star);
          w1j_star = W1_star.unsafe_col(j).elem(index_tmp);
          kappa1j_star = Kappa_star.unsafe_col(j).elem(index_tmp);
          W1j_star = arma::diagmat(w1j_star);
          Z1j_star = kappa1j_star/w1j_star+u10_star;
          V1j_star = U1j_star.t()*W1j_star*U1j_star+Sigma1_star_inv;
          V1j_star_inv=arma::inv_sympd(V1j_star);
          U1j_star = arma::mvnrnd(V1j_star_inv*(U1j_star.t()*W1j_star*Z1j_star+Sigma1_star_inv*(u1_mean_prior_2+Sigma1_part2*(arma::conv_to< arma::vec >::from(u1.submat(j,1,j,r-1)) -u1_mean_prior_1) ) ), V1j_star_inv);
          u1_star.submat(j,1,j,r_star-1) = arma::conv_to< arma::rowvec >::from(U1j_star);
        } else {
          u1_star.submat(j,1,j,r_star-1) = arma::conv_to< arma::rowvec >::from(arma::mvnrnd(u1_mean_prior_2+Sigma1_part2*(arma::conv_to< arma::vec >::from(u1.submat(j,1,j,r-1)) -u1_mean_prior_1) ,Sigma1_star));
        }
      }
      
    // update proposal variance
    if((i+1)%k==0){
      s = (i+1)/k;
      gamma1 = 1.0/pow(s+3,c1);
      gamma2 = 10.0 * gamma1;

      r_theta.col(s-1) = arma::sum(acs_theta.cols(s*k-k,s*k-1),1)/k;
      adaptFactor_vec = exp(gamma2 * (r_theta.col(s-1) - r_opt));
      theta_sd_prop = theta_sd_prop%adaptFactor_vec;
    }

  }

  return List::create(Named("BETA")=BETA.t(),_["SIGMA2"]=SIGMA2,_["THETA"]=THETA.t());
}
