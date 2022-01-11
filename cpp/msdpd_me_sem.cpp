#define ARMA_64BIT_WORD
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

//fill in
sp_mat btri_mat(const sp_mat& all_matrices){
  const int p = all_matrices.n_cols;
  const int tp = all_matrices.n_rows;
  const int t = tp/p;
  //make container
  sp_mat all_col = sp_mat(tp, tp);
  for (int i = 0; i < t; ++i){
    all_col.submat(i*p, i*p, tp - 1, (i + 1)*p - 1) = all_matrices.head_rows((t-i)*p);
  }
  return all_col;
}

//fill in
field<sp_mat> make_A_df(const int& p, const int& t, const sp_mat& dthetas, const sp_mat& lws, const sp_mat& iirws){
  //make containers
  field<sp_mat> output = field<sp_mat>(2);
  sp_mat all_matrices = sp_mat((t + 1)*p, p);
  //build the long block
  const sp_mat diag_p = speye(p, p);
  const sp_mat c_mat = iirws*(dthetas + lws);
  sp_mat tmp_mat = diag_p - c_mat;
  //fill the containers
  all_matrices.rows(0, p - 1) = diag_p;
  all_matrices.rows(p, 2*p - 1) = c_mat - 2*diag_p;
  tmp_mat = tmp_mat * tmp_mat;
  all_matrices.rows(2*p, 3*p - 1) = tmp_mat;
  for (int i = 3; i < t + 1; ++i){
    tmp_mat = c_mat * tmp_mat;
    all_matrices.rows(i*p, (i+1)*p - 1) = tmp_mat;
  }
  const sp_mat mat_proto = btri_mat(all_matrices);
  output(0) = mat_proto.submat(p, 0, (t + 1)*p - 1, t*p - 1);
  output(1) = mat_proto.submat(0, 0, t*p - 1, t*p - 1);
  return output;
}

//transform wide_mat to field
field<sp_mat> wmat2field_sp(const sp_mat& input){
  const int n = input.n_rows;
  const int q = input.n_cols/n;
  field<sp_mat> output(q);
  for (int i=0; i< q; ++i){
    output(i) = input.cols(i*n, (i+1)*n - 1);
  }
  return output;
}

//Use field
sp_mat mat_exp(const field<sp_mat>& mat_arr, const double& par_ast, const int& p, const double& zero_th){
  const int q = mat_arr.n_elem;
  sp_mat output = speye(p,p);
  for(int i=1; i <= q; ++i){
    output += mat_arr(i-1)*pow(par_ast, i)/tgamma(i+1);
  }
  return output.clean(zero_th);
}

// [[Rcpp::export]]
Rcpp::RObject msdpd_me_sem_aqs(const arma::vec& para, 
                        const arma::mat& x_, 
                        const arma::vec& y,
                        const arma::vec& y1, 
                        const arma::sp_mat& w_er,
                        const arma::sp_mat& me_alp,
                        const arma::sp_mat& inv_c, 
                        const bool& correction,
                        const double& zero_th,
                        const bool& sq = true
                        ){
  const int tp = y.size();
  const int p = w_er.n_cols;
  const int t = tp/p;
  const double alp_ast = para(0);
  const double theta = para(1);
  const sp_mat diag_p = speye(p, p);
  const sp_mat diag_t = speye(t, t);
  const sp_mat bdinv_c = kron(inv_c, diag_p);
  const field<sp_mat> me_alp_fld = wmat2field_sp(me_alp);
  const sp_mat iaws = mat_exp(me_alp_fld, alp_ast, p, zero_th);
  const sp_mat iaw = kron(diag_t, iaws);
  const sp_mat diaws = iaws.t()*iaws;
  const sp_mat iaaw = kron(inv_c, diaws);
  const vec beta = (x_.t()*iaaw*x_).i()*x_.t()*iaaw*(y - theta*y1);
  const vec k_ast = y - x_*beta - theta*y1;
  const double sigs = as_scalar(k_ast.t()*iaaw*k_ast)/tp;
  const sp_mat bdw_er = kron(diag_t, w_er);
  const rowvec tmp_mat_1 = k_ast.t()*iaaw;
  vec eq = vec(2);
  if (correction) {
    field<sp_mat> A_mats = make_A_df(p, t, theta*speye(p, p), sp_mat(p, p), speye(p, p));
    eq(0) = 1/sigs*as_scalar(tmp_mat_1*y1) + trace(bdinv_c*A_mats(1));
    eq(1) = 0.5/sigs*as_scalar(k_ast.t()*kron(inv_c, diaws*(w_er+w_er.t()))*k_ast);
  } else {
    eq(0) = 1/sigs*as_scalar(tmp_mat_1*y1);
    eq(1) = 0.5/sigs*as_scalar(k_ast.t()*kron(inv_c, iaws*(w_er+w_er.t()))*k_ast);
  }
  if (sq){
    return Rcpp::wrap(sum(square(eq)));
  } else {
    return Rcpp::wrap(eq);
  }
}
