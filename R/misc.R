input_order = function(mat){
  mat[,1] = factor(mat[,1])
  mat = mat[order(mat[,1]),]
  if (ncol(mat) > 2){
    mat = mat[order(mat[,2]),]
  }
  return(mat)
}


btri_mat_sp = function(mat_list){
  t = length(mat_list)
  tmp_long_mat = do.call(rbind, mat_list)
  n = dim(tmp_long_mat)[2]
  output = tmp_long_mat
  for (i in 2:t){
    output = cbind(output, rbind(Matrix(0, nrow = (i - 1)*n, ncol = n), tmp_long_mat[1:(n*(t - (i - 1))),]))
  }
  return(output)
}

c_mat_pows_sp = function(p, t, iirws, dthetas, lws){
  c_mat = iirws %*% (dthetas + lws)
  output = list()
  temp_mat = Diagonal(p)
  for (i in 0:t){
    output[[as.character(i)]] = temp_mat
    temp_mat = temp_mat %*% c_mat
  }
  return(output)
}

psi_plus_all_sp = function(p, t, psi){
  output = Matrix(0, nrow = p*t, ncol = p)
  for (i in 1:t) {
    output = output + psi[, ((i - 1)*p + 1):(i*p)]
  }
  return(output)
}

make_A_df_sp = function(p, t, dthetas, lws, iirws){
  c_mat = iirws %*% (dthetas + lws)
  m_list = list()
  m_list[[1]] = Diagonal(p)
  for (i in 2:(t+1)){
    if (i == 2){
      m_list[[i]] = c_mat - 2*Diagonal(p)
    }else if(i == 3){
      tmp_mat = Diagonal(p) - c_mat
      tmp_mat = tmp_mat %*% tmp_mat
      m_list[[i]] = tmp_mat
    } else{
      tmp_mat = c_mat %*% tmp_mat
      m_list[[i]] = tmp_mat
    }
  }
  mat_proto = btri_mat_sp(m_list)
  mat_A = mat_proto[-c(1:p), -c((t*p+1):((t+1)*p))]
  mat_A1 = mat_proto[-c((t*p+1):((t+1)*p)), -c((t*p+1):((t+1)*p))]
  return(list("A" = mat_A, 
              "A_1" = mat_A1
  )
  )
}

make_A_deriv_df_sp = function(p, t, w, w_lam, dthetas, lws, iirws, mode){
  c_mat = iirws %*% (dthetas + lws)
  switch(mode,
         "rho" = {
           deri_c_mat = -iirws %*% w %*% (dthetas + lws)
         },
         "theta" = {
           deri_c_mat = iirws
         },
         "lambda" = {
           deri_c_mat = iirws %*% w_lam
         },
         stop("please choose from rho, theta and lambda"))
  sq_mat = (Diagonal(p) - c_mat)%*%(Diagonal(p) - c_mat)
  deri_sq_mat = -((Diagonal(p) - c_mat)%*%deri_c_mat + deri_c_mat %*%(Diagonal(p) - c_mat))
  m_list = list()
  m_list[[1]] = Matrix(0, p, p)
  for (i in 2:(t+1)){
    if (i == 2){
      m_list[[i]] = deri_c_mat
    }else if(i == 3){
      m_list[[i]] = deri_sq_mat
    } else if(i == 4){
      tmp_deri_mat = deri_c_mat
      tmp_pow_mat = c_mat
      m_list[[i]] = tmp_pow_mat%*%deri_sq_mat + tmp_deri_mat%*%sq_mat
    } else{
      tmp_deri_mat = c_mat%*%tmp_deri_mat + deri_c_mat%*%tmp_pow_mat
      tmp_pow_mat = tmp_pow_mat%*%c_mat
      m_list[[i]] = tmp_pow_mat%*%deri_sq_mat + tmp_deri_mat%*%sq_mat
    }
  }
  mat_deriv_proto = btri_mat_sp(m_list)
  mat_A_deriv = mat_deriv_proto[-c(1:p), -c((t*p+1):((t+1)*p))]
  mat_A1_deriv = mat_deriv_proto[-c((t*p+1):((t+1)*p)), -c((t*p+1):((t+1)*p))]
  return(list("A_deriv" = mat_A_deriv,
              "A_deriv_1" = mat_A1_deriv
  )
  )
}

make_R_mats_sp = function(list_pow_c_mat){
  mat_R = bdiag(list_pow_c_mat[-1])
  mat_R_1 = bdiag(list_pow_c_mat[-length(list_pow_c_mat)])
  return(list("R" = mat_R,
              "R_1" = mat_R_1))
}

make_G_mats_sp = function(list_pow_c_mat){
  t = length(list_pow_c_mat) - 1
  p = dim(list_pow_c_mat[["0"]])[1]
  list_pow_c_mat[[t + 1]] = NULL
  list_pow_c_mat = c(list(Matrix(0, p, p)), list_pow_c_mat)
  G_proto = btri_mat_sp(list_pow_c_mat)
  mat_G = G_proto[-c(1:p), -c((t*p+1):((t+1)*p))]
  mat_G_1 = G_proto[-c((t*p+1):((t+1)*p)), -c((t*p+1):((t+1)*p))]
  return(list("G" = mat_G,
              "G_1" = mat_G_1))
}

#slice horizotally to get g1i
make_g1i_all_sp = function(p, t, pi, v_vec){
  output = Matrix(0, nrow = p, ncol = dim(pi)[2])
  temp = pi*v_vec
  for (j in 1:t) {
    output = output + temp[((j-1)*p+1):(j*p),]
  }
  return(output)
}

#contain p elements, i = 1, 2,..., p
make_g2i_all_sp = function(p, t, phi, v_vec, sigs){
  slu_phi = phi + t(phi)
  #first term
  temp1 = numeric(p*t)
  for (i in 1:t){
    temp1_mat_1 = Matrix(0, nrow = p, ncol = t*p)
    for (j in 1:t){
      temp1_mat_2 = slu_phi[((i-1)*p + 1):(i*p),((j-1)*p + 1):(j*p)]
      l_temp1_mat_2 = tril(temp1_mat_2, k = -1)
      temp1_mat_1[, ((j-1)*p + 1):(j*p)] = l_temp1_mat_2
    }
    temp1[((i-1)*p + 1):(i*p)] = temp1_mat_1%*%v_vec
  }
  #second term
  temp2 = numeric(p*t)
  for (i in 1:t){
    temp2_mat_1 = Matrix(0, nrow = p, ncol = t*p)
    for (j in 1:t){
      temp2_mat_2 = phi[((i-1)*p + 1):(i*p),((j-1)*p + 1):(j*p)]
      temp2_mat_1[, ((j-1)*p + 1):(j*p)] = Diagonal(p, diag(temp2_mat_2))
    }
    temp2[((i-1)*p + 1):(i*p)] = temp2_mat_1%*%v_vec
  }
  #third term
  temp3 = -sigs*diag(kronecker(bandSparse(t, k = c(-1,0,1), diagonals = list(rep(-1, t-1), rep(2, t), rep(-1, t-1))), Diagonal(p))%*%phi)
  return(rowSums(Matrix(v_vec*(temp1 + temp2) + temp3, ncol = t)))
}

#contain n elements, i = 1, 2,..., n
make_g3i_all_sp = function(p, t, y0, y0_ast, psi, v_vec, iirws, iiaws, sigs){
  all_psi = psi_plus_all_sp(p, t, psi)
  list_all_psi = list()
  for (i in 1:t){
    list_all_psi[[i]] = all_psi[((i-1)*p + 1):(i*p),]
  }
  Theta = list_all_psi[[1]]%*%iirws%*%iiaws
  u_Theta = triu(Theta,k=1)
  l_Theta = tril(Theta,k=-1)
  y0_psi = as.vector(bdiag(list_all_psi)%*%y0)
  temp1 = v_vec[1:p]*((l_Theta + t(u_Theta))%*%y0_ast) + diag(Theta)*(v_vec[1:p]*y0_ast + sigs)
  temp2 = rowSums(Matrix((v_vec*y0_psi)[-c(1:p)], nrow = p))
  return(temp1 + temp2)
}


df_data = function(dat){
  index = dat[,c(1,2)]
  n = length(unique(index[,1]))
  t = length(unique(index[,2]))
  output = dat[(n+1):(n*t), -c(1,2), drop = F] - dat[1:(n*(t-1)), -c(1,2), drop = F]
  return(cbind(index[-c(1:n),], output))
}

mat_exp_arr = function(mat, q, sp_mode = T){
  mat_dim = dim(mat)
  if (mat_dim[1]==mat_dim[2]){
    output = list()
  }else{
    stop("non-square matrix")
  }
  if (sp_mode){
    temp_mat = Diagonal(mat_dim[1])
  } else {
    temp_mat = diag(mat_dim[1])
  }
  for (i in 1:q){
    temp_mat = temp_mat%*%mat
    output[[i]] = temp_mat
  }
  return(output)
}


mat_exp = function(mat_arr, rho_ast, sp_mode = T, zero_th = 1e-12){
  n = dim(mat_arr[[1]])[1]
  q = length(mat_arr)
  if (sp_mode){
    output = Diagonal(n)
  } else {
    output = diag(n)
  }
  for (i in 1:q){
    output = output + mat_arr[[i]]*rho_ast^i/factorial(i)
  }
  output[abs(output)<zero_th] = 0
  return(output)
}

c_inv = function(n){
  output = matrix(rep(seq(1,n), n), nrow = n, ncol = n)
  output = output%*%diag(seq(n,1))
  output = as.matrix(forceSymmetric(output))
  output = output/(n+1)
  return(output)
}

full_aqs_me = function(para, x_, y, y1, w, w_er, w_lam, me_rho, me_alp, inv_c, correction, mode = "normal", hessian_er, sp_mode, zero_th, sq = T){
  x = x_
  k = dim(x)[2]
  tp = length(y)
  p = dim(w)[1]
  t = tp/p
  rho_ast = para[1]
  alp_ast = para[2]
  theta = para[3]
  lam = para[4]
  eq = numeric(4)
  xt = t(x)
  iaws = mat_exp(me_alp, alp_ast, sp_mode, zero_th)
  iiaws = mat_exp(me_alp, -alp_ast, sp_mode, zero_th)
  irws = mat_exp(me_rho, rho_ast, sp_mode, zero_th)
  iirws = mat_exp(me_rho, -rho_ast, sp_mode, zero_th)
  lws = lam * w_lam
  if (sp_mode){
    bdinv_c = kronecker(inv_c, Diagonal(p))
    irw = kronecker(Diagonal(t), irws)
    iaw = kronecker(Diagonal(t), iaws)
    lw = kronecker(Diagonal(t), lws)
    iirw = kronecker(Diagonal(t), iirws)
    iiaw = kronecker(Diagonal(t), iiaws)
    bdw_lam = kronecker(Diagonal(t), w_lam)
    bdw = kronecker(Diagonal(t), w)
    bddw = irw%*%bdw
  } else {
    bdinv_c = kronecker(inv_c, diag(p))
    irw = kronecker(diag(t), irws)
    iaw = kronecker(diag(t), iaws)
    lw = kronecker(diag(t), lws)
    iirw = kronecker(diag(t), iirws)
    iiaw = kronecker(diag(t), iiaws)
    bdw_lam = kronecker(diag(t), w_lam)
    bdw = kronecker(diag(t), w)
    bddw = irw%*%bdw
  }
  #
  dw_er = t(w_er) + w_er
  diaws = t(iaws)%*%iaws
  iaaw = kronecker(inv_c, diaws)
  beta = solve(xt%*%iaaw%*%x)%*%xt%*%iaaw%*%(irw%*%y - theta*y1 - lw%*%y1)
  k_ast = irw %*% y - x%*%beta - theta*y1 - lw%*%y1
  sigs = as.numeric(t(k_ast) %*% iaaw %*% k_ast/tp)
  switch(mode,
         "normal" = {
           tmp_mat_1 = t(k_ast) %*% iaaw
           if (correction){
             if (sp_mode){
               A_mats = make_A_df_sp(p, t, theta*Diagonal(p), lws, iirws)
             } else {
               A_mats = make_A_df_sp(p, t, diag(theta, p, p), lws, iirws)
             }
             #aqs
             vec_bias_theta = diag(bdinv_c %*% A_mats$A_1 %*% iirw)
             vec_bias_rho = diag(bdinv_c %*% A_mats$A %*% iirw %*% bddw)
             vec_bias_lam = diag(bdinv_c %*% A_mats$A_1 %*% iirw %*% bdw_lam)
             eq[1] = 1/sigs * tmp_mat_1 %*% y1 + sum(vec_bias_theta)
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y + sum(vec_bias_rho)
             eq[3] = 1/sigs * tmp_mat_1 %*% bdw_lam %*% y1 + sum(vec_bias_lam)
             eq[4] = 0.5/sigs * t(k_ast) %*% kronecker(inv_c, diaws%*%dw_er) %*% k_ast
           } else {
             #qs
             eq[1] = 1/sigs * tmp_mat_1 %*% y1
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y
             eq[3] = 1/sigs * tmp_mat_1 %*% bdw_lam %*% y1
             eq[4] = 0.5/sigs * t(k_ast) %*% kronecker(inv_c, diaws%*%dw_er) %*% k_ast
           }
           if (sq) {
             return(sum(eq^2))
           } else {
             return(eq)
           }
         },
         "beta_sigs" = {
           return(list(beta = beta,
                       sigma2 = sigs))
         },
         "opmd" = {
           #calculate second derivative
           sec_deri = Matrix(0, nrow = k + 5, ncol = k + 5)
           #calculate derivative of A_mat
           dthetas = diag(theta, p, p)
           A_mats = make_A_df_sp(p, t, dthetas, lws, iirws)
           A_mats_rho = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "rho")
           A_mats_theta = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "theta")
           A_mats_lam = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "lambda")
           #lines beta
           sec_deri[1:k, 1:k] = -1/sigs*t(x)%*%iaaw%*%x
           sec_deri[1:k, k+1] = -1/sigs^2*t(x)%*%iaaw%*%k_ast
           sec_deri[1:k, k+2] = -1/sigs*t(x)%*%iaaw%*%y1
           sec_deri[1:k, k+3] = 1/sigs*t(x)%*%iaaw%*%bddw%*%y
           sec_deri[1:k, k+4] = -1/sigs*t(x)%*%iaaw%*%bdw_lam%*%y1
           sec_deri[1:k, k+5] = 1/sigs*t(x)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line sigma2
           sec_deri[k+1, 1:k] = t(sec_deri[1:k, k+1])
           sec_deri[k+1, k+1] = -1/sigs^3*t(k_ast)%*%iaaw%*%k_ast + tp/(2*sigs^2)
           sec_deri[k+1, k+2] = -1/sigs^2*t(y1)%*%iaaw%*%k_ast
           sec_deri[k+1, k+3] = 1/sigs^2*t(y)%*%t(bddw)%*%iaaw%*%k_ast
           sec_deri[k+1, k+4] = -1/sigs^2*t(y1)%*%t(bdw_lam)%*%iaaw%*%k_ast
           sec_deri[k+1, k+5] = 1/(2*sigs^2)*t(k_ast)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line theta
           sec_deri[k+2, 1:k] = t(sec_deri[1:k, k+2])
           sec_deri[k+2, k+1] = sec_deri[k+1, k+2]
           sec_deri[k+2, k+2] = -1/sigs*t(y1)%*%iaaw%*%y1 + sum(diag(bdinv_c%*%A_mats_theta$A_deriv_1%*%iirw))
           sec_deri[k+2, k+3] = 1/sigs*t(y1)%*%iaaw%*%bddw%*%y + sum(diag(bdinv_c%*%(A_mats_rho$A_deriv_1%*%iirw - A_mats$A_1%*%iirw%*%bdw)))
           sec_deri[k+2, k+4] = -1/sigs*t(y1)%*%iaaw%*%bdw_lam%*%y1 + sum(diag(bdinv_c%*%A_mats_lam$A_deriv_1%*%iirw))
           sec_deri[k+2, k+5] = 1/sigs*t(y1)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line rho
           sec_deri[k+3, 1:k] = t(sec_deri[1:k, k+3])
           sec_deri[k+3, k+1] = sec_deri[k+1, k+3]
           sec_deri[k+3, k+2] = sec_deri[k+2, k+3]
           sec_deri[k+3, k+3] = -1/sigs*(t(k_ast)%*%iaaw%*%bddw%*%bdw%*%y + t(bddw%*%y)%*%iaaw%*%bddw%*%y) 
                                - sum(diag(bdinv_c%*%((A_mats_rho$A_deriv%*%iirw - A_mats$A%*%iirw%*%bdw)%*%bddw + A_mats$A%*%iirw%*%bddw%*%bdw)))
           sec_deri[k+3, k+4] = 1/sigs*t(y)%*%t(bddw)%*%iaaw%*%bdw_lam%*%y1 - sum(diag(bdinv_c%*%A_mats_lam$A_deriv%*%iirw%*%bddw))
           sec_deri[k+3, k+5] = -1/sigs*t(y)%*%t(bddw)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line lambda
           sec_deri[k+4, 1:k] = t(sec_deri[1:k, k+4])
           sec_deri[k+4, k+1] = sec_deri[k+1, k+4]
           sec_deri[k+4, k+2] = sec_deri[k+2, k+4] 
           sec_deri[k+4, k+3] = sec_deri[k+3, k+4] 
           sec_deri[k+4, k+4] = -1/sigs*t(y1)%*%t(bdw_lam)%*%iaaw%*%bdw_lam%*%y1 + sum(diag(bdinv_c%*%A_mats_lam$A_deriv_1%*%iirw%*%bdw_lam))
           sec_deri[k+4, k+5] = 1/sigs*t(y1)%*%t(bdw_lam)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line alpha
           sec_deri[k+5, 1:k] = t(sec_deri[1:k, k+5])
           sec_deri[k+5, k+1] = sec_deri[k+1, k+5]
           sec_deri[k+5, k+2] = sec_deri[k+2, k+5]
           sec_deri[k+5, k+3] = sec_deri[k+3, k+5]
           sec_deri[k+5, k+4] = sec_deri[k+4, k+5]
           sec_deri[k+5, k+5] = -0.5/sigs*t(k_ast)%*%(kronecker(inv_c, diaws%*%dw_er%*%dw_er))%*%k_ast
           #make list of powers of c_mat
           list_pow_c_mat = c_mat_pows_sp(p,t, iirws, dthetas, lws)
           #make R and G mats
           R_mats = make_R_mats_sp(list_pow_c_mat)
           G_mats = make_G_mats_sp(list_pow_c_mat)
           #make eta, S, pi, phi, psi and v_vec
           ihaws = iaws
           iihaws = iiaws
           ihaw = iaw
           iihaw = iiaw
           iciaw = kronecker(inv_c, ihaws)
           eta = G_mats$G %*% iirw %*% x %*% beta
           eta_1 = G_mats$G_1 %*% iirw %*% x %*% beta
           S = G_mats$G %*% iirw %*% iihaw
           S_1 = G_mats$G_1 %*% iirw %*% iihaw
           pi_1 = 1/sigs*t(iciaw)%*%x
           pi_2 = 1/sigs*t(iciaw)%*%eta_1
           pi_3 = 1/sigs*t(iciaw)%*%bddw%*%eta
           pi_4 = 1/sigs*t(iciaw)%*%bdw_lam%*%eta_1
           phi_1 = 0.5/(sigs^2)*kronecker(inv_c, Diagonal(p))
           phi_2 = 1/sigs*iciaw%*%S_1
           phi_3 = 1/sigs*iciaw%*%bddw%*%S
           phi_4 = 1/sigs*iciaw%*%bdw_lam%*%S_1
           phi_5 = 0.5/(sigs)*kronecker(inv_c, dw_er)
           psi_1 = 1/sigs*iciaw%*%R_mats$R_1
           psi_2 = 1/sigs*iciaw%*%bddw%*%R_mats$R
           psi_3 = 1/sigs*iciaw%*%bdw_lam%*%R_mats$R_1
           v_vec = as.vector(ihaw%*%k_ast)
           #make y0, y0_ast
           y0 = rep(y1[1:p], t)
           y0_ast = ihaws %*% irws %*% y0[1:p]
           #make all gs
           g11 = make_g1i_all_sp(p, t, pi_1, v_vec)
           g12 = make_g1i_all_sp(p, t, pi_2, v_vec)
           g13 = make_g1i_all_sp(p, t, pi_3, v_vec)
           g14 = make_g1i_all_sp(p, t, pi_4, v_vec)
           g21 = make_g2i_all_sp(p, t, phi_1, v_vec, sigs)
           g22 = make_g2i_all_sp(p, t, phi_2, v_vec, sigs)
           g23 = make_g2i_all_sp(p, t, phi_3, v_vec, sigs)
           g24 = make_g2i_all_sp(p, t, phi_4, v_vec, sigs)
           g25 = make_g2i_all_sp(p, t, phi_5, v_vec, sigs)
           g31 = make_g3i_all_sp(p, t, y0, y0_ast, psi_1, v_vec, iirws, iihaws, sigs)
           g32 = make_g3i_all_sp(p, t, y0, y0_ast, psi_2, v_vec, iirws, iihaws, sigs)
           g33 = make_g3i_all_sp(p, t, y0, y0_ast, psi_3, v_vec, iirws, iihaws, sigs)
           #make all g mat
           all_g = Matrix(0, nrow = k + 5, ncol = p)
           all_g[1:k, ] = t(g11)
           all_g[k + 1, ] = g21
           all_g[k + 2, ] = g31 + g12 + g22
           all_g[k + 3, ] = -(g32 + g13 + g23)
           all_g[k + 4, ] = g33 + g14 + g24
           all_g[k + 5, ] = -g25
           inv_sig_mat = solve(-sec_deri)
           gamma_mat = all_g%*%t(all_g)
           if (hessian_er){
             return_list = list(vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat,
                                hes_mat = sec_deri)
             return(return_list)
           } else {
             vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat
             return(list(vc_mat = vc_mat))
           }
         },
         stop("undefined mode"))
}

slm_aqs_me = function(para, x_, y, y1, w, me_rho, inv_c, correction, mode = "normal", hessian_er, sp_mode, zero_th){
  x = x_
  k = dim(x)[2]
  tp = length(y)
  p = dim(w)[1]
  t = tp/p
  rho_ast = para[1]
  theta = para[2]
  eq = numeric(2)
  xt = t(x)
  irws = mat_exp(me_rho, rho_ast, sp_mode, zero_th)
  iirws = mat_exp(me_rho, -rho_ast, sp_mode, zero_th)
  if (sp_mode){
    bdinv_c = kronecker(inv_c, Diagonal(p))
    irw = kronecker(Diagonal(t), irws)
    iirw = kronecker(Diagonal(t), iirws)
    bdw = kronecker(Diagonal(t), w)
    bddw = irw%*%bdw
  } else {
    bdinv_c = kronecker(inv_c, diag(p))
    irw = kronecker(diag(t), irws)
    iirw = kronecker(diag(t), iirws)
    bdw = kronecker(diag(t), w)
    bddw = irw%*%bdw
  }
  iaaw = bdinv_c
  beta = solve(xt%*%iaaw%*%x)%*%xt%*%iaaw%*%(irw%*%y - theta*y1)
  k_ast = irw %*% y - x%*%beta - theta*y1
  sigs = as.numeric(t(k_ast) %*% iaaw %*% k_ast/tp)
  switch(mode,
         "normal" = {
           tmp_mat_1 = t(k_ast) %*% iaaw
           if (correction){
             if (sp_mode){
               A_mats = make_A_df_sp(p, t, theta*Diagonal(p), 0, iirws)
             } else {
               A_mats = make_A_df_sp(p, t, diag(theta, p, p), 0, iirws)
             }
             #aqs
             vec_bias_theta = diag(bdinv_c %*% A_mats$A_1 %*% iirw)
             vec_bias_rho = diag(bdinv_c %*% A_mats$A %*% iirw %*% bddw)
             eq[1] = 1/sigs * tmp_mat_1 %*% y1 + sum(vec_bias_theta)
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y + sum(vec_bias_rho)
           } else {
             #qs
             eq[1] = 1/sigs * tmp_mat_1 %*% y1
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y
           }
           return(sum(eq^2))
         },
         "beta_sigs" = {
           return(list(beta = beta,
                       sigma2 = sigs))
         },
         "opmd" = {
           #calculate second derivative
           sec_deri = Matrix(0, nrow = k + 3, ncol = k + 3)
           #calculate derivative of A_mat
           dthetas = diag(theta, p, p)
           A_mats = make_A_df_sp(p, t, dthetas, 0, iirws)
           A_mats_rho = make_A_deriv_df_sp(p, t, w, 0, dthetas, 0, iirws, mode = "rho")
           A_mats_theta = make_A_deriv_df_sp(p, t, w, 0, dthetas, 0, iirws, mode = "theta")
           #lines beta
           sec_deri[1:k, 1:k] = -1/sigs*t(x)%*%iaaw%*%x
           sec_deri[1:k, k+1] = -1/sigs^2*t(x)%*%iaaw%*%k_ast
           sec_deri[1:k, k+2] = -1/sigs*t(x)%*%iaaw%*%y1
           sec_deri[1:k, k+3] = 1/sigs*t(x)%*%iaaw%*%bddw%*%y
           #line sigma2
           sec_deri[k+1, 1:k] = t(sec_deri[1:k, k+1])
           sec_deri[k+1, k+1] = -1/sigs^3*t(k_ast)%*%iaaw%*%k_ast + tp/(2*sigs^2)
           sec_deri[k+1, k+2] = -1/sigs^2*t(y1)%*%iaaw%*%k_ast
           sec_deri[k+1, k+3] = 1/sigs^2*t(y)%*%t(bddw)%*%iaaw%*%k_ast
           #line theta
           sec_deri[k+2, 1:k] = t(sec_deri[1:k, k+2])
           sec_deri[k+2, k+1] = sec_deri[k+1, k+2]
           sec_deri[k+2, k+2] = -1/sigs*t(y1)%*%iaaw%*%y1 + sum(diag(bdinv_c%*%A_mats_theta$A_deriv_1%*%iirw))
           sec_deri[k+2, k+3] = 1/sigs*t(y1)%*%iaaw%*%bddw%*%y + sum(diag(bdinv_c%*%(A_mats_rho$A_deriv_1%*%iirw - A_mats$A_1%*%iirw%*%bdw)))
           #line rho
           sec_deri[k+3, 1:k] = t(sec_deri[1:k, k+3])
           sec_deri[k+3, k+1] = sec_deri[k+1, k+3]
           sec_deri[k+3, k+2] = sec_deri[k+2, k+3]
           sec_deri[k+3, k+3] = -1/sigs*(t(k_ast)%*%iaaw%*%bddw%*%bdw%*%y + t(bddw%*%y)%*%iaaw%*%bddw%*%y) 
           - sum(diag(bdinv_c%*%((A_mats_rho$A_deriv%*%iirw - A_mats$A%*%iirw%*%bdw)%*%bddw + A_mats$A%*%iirw%*%bddw%*%bdw)))
           #make list of powers of c_mat
           list_pow_c_mat = c_mat_pows_sp(p,t, iirws, dthetas, Matrix(0,p,p))
           #make R and G mats
           R_mats = make_R_mats_sp(list_pow_c_mat)
           G_mats = make_G_mats_sp(list_pow_c_mat)
           #make eta, S, pi, phi, psi and v_vec
           iciaw = bdinv_c
           eta = G_mats$G %*% iirw %*% x %*% beta
           eta_1 = G_mats$G_1 %*% iirw %*% x %*% beta
           S = G_mats$G %*% iirw
           S_1 = G_mats$G_1 %*% iirw
           pi_1 = 1/sigs*t(iciaw)%*%x
           pi_2 = 1/sigs*t(iciaw)%*%eta_1
           pi_3 = 1/sigs*t(iciaw)%*%bddw%*%eta
           phi_1 = 0.5/(sigs^2)*kronecker(inv_c, Diagonal(p))
           phi_2 = 1/sigs*iciaw%*%S_1
           phi_3 = 1/sigs*iciaw%*%bddw%*%S
           psi_1 = 1/sigs*iciaw%*%R_mats$R_1
           psi_2 = 1/sigs*iciaw%*%bddw%*%R_mats$R
           v_vec = as.vector(k_ast)
           #make y0, y0_ast
           y0 = rep(y1[1:p], t)
           y0_ast = irws %*% y0[1:p]
           #make all gs
           g11 = make_g1i_all_sp(p, t, pi_1, v_vec)
           g12 = make_g1i_all_sp(p, t, pi_2, v_vec)
           g13 = make_g1i_all_sp(p, t, pi_3, v_vec)
           g21 = make_g2i_all_sp(p, t, phi_1, v_vec, sigs)
           g22 = make_g2i_all_sp(p, t, phi_2, v_vec, sigs)
           g23 = make_g2i_all_sp(p, t, phi_3, v_vec, sigs)
           g31 = make_g3i_all_sp(p, t, y0, y0_ast, psi_1, v_vec, iirws, Diagonal(p), sigs)
           g32 = make_g3i_all_sp(p, t, y0, y0_ast, psi_2, v_vec, iirws, Diagonal(p), sigs)
           #make all g mat
           all_g = Matrix(0, nrow = k + 3, ncol = p)
           all_g[1:k, ] = t(g11)
           all_g[k + 1, ] = g21
           all_g[k + 2, ] = g31 + g12 + g22
           all_g[k + 3, ] = -(g32 + g13 + g23)
           inv_sig_mat = solve(-sec_deri)
           gamma_mat = all_g%*%t(all_g)
           if (hessian_er){
             return_list = list(vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat,
                                hes_mat = sec_deri)
             return(return_list)
           } else {
             vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat
             return(list(vc_mat = vc_mat))
           }
         },
         stop("undefined mode"))
}

sem_aqs_me = function(para, x_, y, y1, w_er, me_alp, inv_c, correction, mode = "normal", hessian_er, sp_mode, zero_th){
  x = x_
  k = dim(x)[2]
  tp = length(y)
  p = dim(w_er)[1]
  t = tp/p
  alp_ast = para[1]
  theta = para[2]
  eq = numeric(2)
  xt = t(x)
  iaws = mat_exp(me_alp, alp_ast, sp_mode, zero_th)
  iiaws = mat_exp(me_alp, -alp_ast, sp_mode, zero_th)
  if (sp_mode){
    bdinv_c = kronecker(inv_c, Diagonal(p))
    iaw = kronecker(Diagonal(t), iaws)
    iiaw = kronecker(Diagonal(t), iiaws)
  } else {
    bdinv_c = kronecker(inv_c, diag(p))
    iaw = kronecker(diag(t), iaws)
    iiaw = kronecker(diag(t), iiaws)
  }
  #
  dw_er = t(w_er) + w_er
  diaws = t(iaws)%*%iaws
  iaaw = kronecker(inv_c, diaws)
  beta = solve(xt%*%iaaw%*%x)%*%xt%*%iaaw%*%(y - theta*y1)
  k_ast = y - x%*%beta - theta*y1
  sigs = as.numeric(t(k_ast) %*% iaaw %*% k_ast/tp)
  switch(mode,
         "normal" = {
           tmp_mat_1 = t(k_ast) %*% iaaw
           if (correction){
             if (sp_mode){
               A_mats = make_A_df_sp(p, t, theta*Diagonal(p), 0, Diagonal(p))
             } else {
               A_mats = make_A_df_sp(p, t, diag(theta, p, p), 0, diag(p))
             }
             #aqs
             vec_bias_theta = diag(bdinv_c %*% A_mats$A_1)
             eq[1] = 1/sigs * tmp_mat_1 %*% y1 + sum(vec_bias_theta)
             eq[2] = 0.5/sigs * t(k_ast) %*% kronecker(inv_c, diaws%*%dw_er) %*% k_ast
           } else {
             #qs
             eq[1] = 1/sigs * tmp_mat_1 %*% y1
             eq[2] = 0.5/sigs * t(k_ast) %*% kronecker(inv_c, diaws%*%dw_er) %*% k_ast
           }
           return(sum(eq^2))
         },
         "beta_sigs" = {
           return(list(beta = beta,
                       sigma2 = sigs))
         },
         "opmd" = {
           #calculate second derivative
           sec_deri = Matrix(0, nrow = k + 3, ncol = k + 3)
           #calculate derivative of A_mat
           dthetas = diag(theta, p, p)
           A_mats = make_A_df_sp(p, t, dthetas, 0, Diagonal(p))
           A_mats_theta = make_A_deriv_df_sp(p, t, Matrix(0,p,p), Matrix(0,p,p), dthetas, 0, Diagonal(p), mode = "theta")
           #lines beta
           sec_deri[1:k, 1:k] = -1/sigs*t(x)%*%iaaw%*%x
           sec_deri[1:k, k+1] = -1/sigs^2*t(x)%*%iaaw%*%k_ast
           sec_deri[1:k, k+2] = -1/sigs*t(x)%*%iaaw%*%y1
           sec_deri[1:k, k+3] = 1/sigs*t(x)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           # line sigma2
           sec_deri[k+1, 1:k] = t(sec_deri[1:k, k+1])
           sec_deri[k+1, k+1] = -1/sigs^3*t(k_ast)%*%iaaw%*%k_ast + tp/(2*sigs^2)
           sec_deri[k+1, k+2] = -1/sigs^2*t(y1)%*%iaaw%*%k_ast
           sec_deri[k+1, k+3] = 1/(2*sigs^2)*t(k_ast)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line theta
           sec_deri[k+2, 1:k] = t(sec_deri[1:k, k+2])
           sec_deri[k+2, k+1] = sec_deri[k+1, k+2]
           sec_deri[k+2, k+2] = -1/sigs*t(y1)%*%iaaw%*%y1 + sum(diag(bdinv_c%*%A_mats_theta$A_deriv_1%*%Diagonal(t*p)))
           sec_deri[k+2, k+3] = 1/sigs*t(y1)%*%kronecker(inv_c, diaws%*%dw_er)%*%k_ast
           #line alpha
           sec_deri[k+3, 1:k] = t(sec_deri[1:k, k+3])
           sec_deri[k+3, k+1] = sec_deri[k+1, k+3]
           sec_deri[k+3, k+2] = sec_deri[k+2, k+3]
           sec_deri[k+3, k+3] = -0.5/sigs*t(k_ast)%*%(kronecker(inv_c, diaws%*%dw_er%*%dw_er))%*%k_ast
           #make list of powers of c_mat
           list_pow_c_mat = c_mat_pows_sp(p,t, Diagonal(p), dthetas, 0)
           #make R and G mats
           R_mats = make_R_mats_sp(list_pow_c_mat)
           G_mats = make_G_mats_sp(list_pow_c_mat)
           #make eta, S, pi, phi, psi and v_vec
           ihaws = iaws
           iihaws = iiaws
           ihaw = iaw
           iihaw = iiaw
           iciaw = kronecker(inv_c, ihaws)
           eta = G_mats$G %*% x %*% beta
           eta_1 = G_mats$G_1 %*% x %*% beta
           S = G_mats$G %*% iihaw
           S_1 = G_mats$G_1 %*% iihaw
           pi_1 = 1/sigs*t(iciaw)%*%x
           pi_2 = 1/sigs*t(iciaw)%*%eta_1
           phi_1 = 0.5/(sigs^2)*kronecker(inv_c, Diagonal(p))
           phi_2 = 1/sigs*iciaw%*%S_1
           phi_5 = 0.5/(sigs)*kronecker(inv_c, dw_er)
           psi_1 = 1/sigs*iciaw%*%R_mats$R_1
           v_vec = as.vector(ihaw%*%k_ast)
           #make y0, y0_ast
           y0 = rep(y1[1:p], t)
           y0_ast = ihaws %*% y0[1:p]
           #make all gs
           g11 = make_g1i_all_sp(p, t, pi_1, v_vec)
           g12 = make_g1i_all_sp(p, t, pi_2, v_vec)
           g21 = make_g2i_all_sp(p, t, phi_1, v_vec, sigs)
           g22 = make_g2i_all_sp(p, t, phi_2, v_vec, sigs)
           g25 = make_g2i_all_sp(p, t, phi_5, v_vec, sigs)
           g31 = make_g3i_all_sp(p, t, y0, y0_ast, psi_1, v_vec, Diagonal(p), iihaws, sigs)
           #make all g mat
           all_g = Matrix(0, nrow = k + 3, ncol = p)
           all_g[1:k, ] = t(g11)
           all_g[k + 1, ] = g21
           all_g[k + 2, ] = g31 + g12 + g22
           all_g[k + 3, ] = -g25
           inv_sig_mat = solve(-sec_deri)
           gamma_mat = all_g%*%t(all_g)
           if (hessian_er){
             return_list = list(vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat,
                                hes_mat = sec_deri)
             return(return_list)
           } else {
             vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat
             return(list(vc_mat = vc_mat))
           }
         },
         stop("undefined mode"))
}


sltl_aqs_me = function(para, x_, y, y1, w, w_lam, me_rho, inv_c, correction, mode = "normal", hessian_er, sp_mode, zero_th){
  x = x_
  k = dim(x)[2]
  tp = length(y)
  p = dim(w)[1]
  t = tp/p
  rho_ast = para[1]
  theta = para[2]
  lam = para[3]
  eq = numeric(3)
  xt = t(x)
  irws = mat_exp(me_rho, rho_ast, sp_mode, zero_th)
  iirws = mat_exp(me_rho, -rho_ast, sp_mode, zero_th)
  lws = lam * w_lam
  if (sp_mode){
    bdinv_c = kronecker(inv_c, Diagonal(p))
    irw = kronecker(Diagonal(t), irws)
    lw = kronecker(Diagonal(t), lws)
    iirw = kronecker(Diagonal(t), iirws)
    bdw_lam = kronecker(Diagonal(t), w_lam)
    bdw = kronecker(Diagonal(t), w)
    bddw = irw%*%bdw
  } else {
    bdinv_c = kronecker(inv_c, diag(p))
    irw = kronecker(diag(t), irws)
    lw = kronecker(diag(t), lws)
    iirw = kronecker(diag(t), iirws)
    bdw_lam = kronecker(diag(t), w_lam)
    bdw = kronecker(diag(t), w)
    bddw = irw%*%bdw
  }
  iaaw = bdinv_c
  beta = solve(xt%*%iaaw%*%x)%*%xt%*%iaaw%*%(irw%*%y - theta*y1 - lw%*%y1)
  k_ast = irw %*% y - x%*%beta - theta*y1 - lw%*%y1
  sigs = as.numeric(t(k_ast) %*% iaaw %*% k_ast/tp)
  switch(mode,
         "normal" = {
           tmp_mat_1 = t(k_ast) %*% iaaw
           if (correction){
             if (sp_mode){
               A_mats = make_A_df_sp(p, t, theta*Diagonal(p), lws, iirws)
             } else {
               A_mats = make_A_df_sp(p, t, diag(theta, p, p), lws, iirws)
             }
             #aqs
             vec_bias_theta = diag(bdinv_c %*% A_mats$A_1 %*% iirw)
             vec_bias_rho = diag(bdinv_c %*% A_mats$A %*% iirw %*% bddw)
             vec_bias_lam = diag(bdinv_c %*% A_mats$A_1 %*% iirw %*% bdw_lam)
             eq[1] = 1/sigs * tmp_mat_1 %*% y1 + sum(vec_bias_theta)
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y + sum(vec_bias_rho)
             eq[3] = 1/sigs * tmp_mat_1 %*% bdw_lam %*% y1 + sum(vec_bias_lam)
           } else {
             #qs
             eq[1] = 1/sigs * tmp_mat_1 %*% y1
             eq[2] = 1/sigs * tmp_mat_1 %*% bddw %*% y
             eq[3] = 1/sigs * tmp_mat_1 %*% bdw_lam %*% y1
           }
           return(sum(eq^2))
         },
         "beta_sigs" = {
           return(list(beta = beta,
                       sigma2 = sigs))
         },
         "opmd" = {
           #calculate second derivative
           sec_deri = Matrix(0, nrow = k + 4, ncol = k + 4)
           #calculate derivative of A_mat
           dthetas = diag(theta, p, p)
           A_mats = make_A_df_sp(p, t, dthetas, lws, iirws)
           A_mats_rho = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "rho")
           A_mats_theta = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "theta")
           A_mats_lam = make_A_deriv_df_sp(p, t, w, w_lam, dthetas, lws, iirws, mode = "lambda")
           #lines beta
           sec_deri[1:k, 1:k] = -1/sigs*t(x)%*%iaaw%*%x
           sec_deri[1:k, k+1] = -1/sigs^2*t(x)%*%iaaw%*%k_ast
           sec_deri[1:k, k+2] = -1/sigs*t(x)%*%iaaw%*%y1
           sec_deri[1:k, k+3] = 1/sigs*t(x)%*%iaaw%*%bddw%*%y
           sec_deri[1:k, k+4] = -1/sigs*t(x)%*%iaaw%*%bdw_lam%*%y1
           #line sigma2
           sec_deri[k+1, 1:k] = t(sec_deri[1:k, k+1])
           sec_deri[k+1, k+1] = -1/sigs^3*t(k_ast)%*%iaaw%*%k_ast + tp/(2*sigs^2)
           sec_deri[k+1, k+2] = -1/sigs^2*t(y1)%*%iaaw%*%k_ast
           sec_deri[k+1, k+3] = 1/sigs^2*t(y)%*%t(bddw)%*%iaaw%*%k_ast
           sec_deri[k+1, k+4] = -1/sigs^2*t(y1)%*%t(bdw_lam)%*%iaaw%*%k_ast
           #line theta
           sec_deri[k+2, 1:k] = t(sec_deri[1:k, k+2])
           sec_deri[k+2, k+1] = sec_deri[k+1, k+2]
           sec_deri[k+2, k+2] = -1/sigs*t(y1)%*%iaaw%*%y1 + sum(diag(bdinv_c%*%A_mats_theta$A_deriv_1%*%iirw))
           sec_deri[k+2, k+3] = 1/sigs*t(y1)%*%iaaw%*%bddw%*%y + sum(diag(bdinv_c%*%(A_mats_rho$A_deriv_1%*%iirw - A_mats$A_1%*%iirw%*%bdw)))
           sec_deri[k+2, k+4] = -1/sigs*t(y1)%*%iaaw%*%bdw_lam%*%y1 + sum(diag(bdinv_c%*%A_mats_lam$A_deriv_1%*%iirw))
           #line rho
           sec_deri[k+3, 1:k] = t(sec_deri[1:k, k+3])
           sec_deri[k+3, k+1] = sec_deri[k+1, k+3]
           sec_deri[k+3, k+2] = sec_deri[k+2, k+3]
           sec_deri[k+3, k+3] = -1/sigs*(t(k_ast)%*%iaaw%*%bddw%*%bdw%*%y + t(bddw%*%y)%*%iaaw%*%bddw%*%y) 
           - sum(diag(bdinv_c%*%((A_mats_rho$A_deriv%*%iirw - A_mats$A%*%iirw%*%bdw)%*%bddw + A_mats$A%*%iirw%*%bddw%*%bdw)))
           sec_deri[k+3, k+4] = 1/sigs*t(y)%*%t(bddw)%*%iaaw%*%bdw_lam%*%y1 - sum(diag(bdinv_c%*%A_mats_lam$A_deriv%*%iirw%*%bddw))
           #line lambda
           sec_deri[k+4, 1:k] = t(sec_deri[1:k, k+4])
           sec_deri[k+4, k+1] = sec_deri[k+1, k+4]
           sec_deri[k+4, k+2] = sec_deri[k+2, k+4] 
           sec_deri[k+4, k+3] = sec_deri[k+3, k+4] 
           sec_deri[k+4, k+4] = -1/sigs*t(y1)%*%t(bdw_lam)%*%iaaw%*%bdw_lam%*%y1 + sum(diag(bdinv_c%*%A_mats_lam$A_deriv_1%*%iirw%*%bdw_lam))
           #make list of powers of c_mat
           list_pow_c_mat = c_mat_pows_sp(p,t, iirws, dthetas, lws)
           #make R and G mats
           R_mats = make_R_mats_sp(list_pow_c_mat)
           G_mats = make_G_mats_sp(list_pow_c_mat)
           #make eta, S, pi, phi, psi and v_vec
           iciaw = bdinv_c
           eta = G_mats$G %*% iirw %*% x %*% beta
           eta_1 = G_mats$G_1 %*% iirw %*% x %*% beta
           S = G_mats$G %*% iirw
           S_1 = G_mats$G_1 %*% iirw
           pi_1 = 1/sigs*t(iciaw)%*%x
           pi_2 = 1/sigs*t(iciaw)%*%eta_1
           pi_3 = 1/sigs*t(iciaw)%*%bddw%*%eta
           pi_4 = 1/sigs*t(iciaw)%*%bdw_lam%*%eta_1
           phi_1 = 0.5/(sigs^2)*kronecker(inv_c, Diagonal(p))
           phi_2 = 1/sigs*iciaw%*%S_1
           phi_3 = 1/sigs*iciaw%*%bddw%*%S
           phi_4 = 1/sigs*iciaw%*%bdw_lam%*%S_1
           psi_1 = 1/sigs*iciaw%*%R_mats$R_1
           psi_2 = 1/sigs*iciaw%*%bddw%*%R_mats$R
           psi_3 = 1/sigs*iciaw%*%bdw_lam%*%R_mats$R_1
           v_vec = as.vector(k_ast)
           #make y0, y0_ast
           y0 = rep(y1[1:p], t)
           y0_ast = irws %*% y0[1:p]
           #make all gs
           g11 = make_g1i_all_sp(p, t, pi_1, v_vec)
           g12 = make_g1i_all_sp(p, t, pi_2, v_vec)
           g13 = make_g1i_all_sp(p, t, pi_3, v_vec)
           g14 = make_g1i_all_sp(p, t, pi_4, v_vec)
           g21 = make_g2i_all_sp(p, t, phi_1, v_vec, sigs)
           g22 = make_g2i_all_sp(p, t, phi_2, v_vec, sigs)
           g23 = make_g2i_all_sp(p, t, phi_3, v_vec, sigs)
           g24 = make_g2i_all_sp(p, t, phi_4, v_vec, sigs)
           g31 = make_g3i_all_sp(p, t, y0, y0_ast, psi_1, v_vec, iirws, Diagonal(p), sigs)
           g32 = make_g3i_all_sp(p, t, y0, y0_ast, psi_2, v_vec, iirws, Diagonal(p), sigs)
           g33 = make_g3i_all_sp(p, t, y0, y0_ast, psi_3, v_vec, iirws, Diagonal(p), sigs)
           #make all g mat
           all_g = Matrix(0, nrow = k + 4, ncol = p)
           all_g[1:k, ] = t(g11)
           all_g[k + 1, ] = g21
           all_g[k + 2, ] = g31 + g12 + g22
           all_g[k + 3, ] = -(g32 + g13 + g23)
           all_g[k + 4, ] = g33 + g14 + g24
           inv_sig_mat = solve(-sec_deri)
           gamma_mat = all_g%*%t(all_g)
           if (hessian_er){
             return_list = list(vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat,
                                hes_mat = sec_deri)
             return(return_list)
           } else {
             vc_mat = inv_sig_mat%*%gamma_mat%*%inv_sig_mat
             return(list(vc_mat = vc_mat))
           }
         },
         stop("undefined mode"))
}