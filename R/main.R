msdpd_me = function(y,
                 x, 
                 q,
                 w1, 
                 correction = TRUE,
                 hessian_er = FALSE,
                 max_try = 5, 
                 w2 = w1, 
                 w3 = w1, 
                 no_tf = FALSE,
                 model = "full", 
                 rcpp = FALSE, 
                 sp_mode = T,
                 zero_th = -Inf,
                 cma_pop_multi = 1,
                 solver = "rCMA"){
  p = length(w1[1,])
  y = df_data(input_order(y))
  x = df_data(input_order(x))
  tp = dim(y)[1]
  t = tp / p
  if (t < 3) stop("Time period less than 4") 
  t1 = t-1
  tp1 = t1*p
  y_1 = as.matrix(y[-c((tp1+1):tp),-c(1,2)])
  y = as.matrix(y[-c(1:p),-c(1,2)])
  x = as.matrix(x[,-c(1,2)])
  x = x[-c(1:p),,drop = F]
  k_o = dim(x)[2]
  if (!no_tf){
    #individual (time)
    tif = as.matrix(rep(1, p))
    tif = kronecker(diag(t1), tif)
    x = cbind(x, tif)
  }
  k = dim(x)[2]
  mat_c = bandSparse(t1, k = c(-1,0,1), diagonals = list(rep(-1, t1-1), rep(2, t1), rep(-1, t1-1)))
  if (!sp_mode) mat_c = as.matrix(mat_c)
  inv_c = c_inv(t1)
  if (model %in% c("slm", "sltl", "full")) {
    me_rho = mat_exp_arr(w1, q[1], sp_mode)
    if (rcpp) me_rho_cpp = do.call(cbind, me_rho)
  }  
  if (model %in% c("sem", "full"))  {
    me_alp = mat_exp_arr(w3, q[2], sp_mode)
    if (rcpp) me_alp_cpp = do.call(cbind, me_alp)
  }
  const_rcma = function(para){
    rho_ast = para[1]
    alp_ast = para[2]
    theta = para[3]
    lam = para[4]
    if (alp_ast >= 2*log(2)| abs(1 - exp(rho_ast)) + abs(theta) + abs(lam) >= 1){
      return(F)
    } else{
      return(T)
    }
  }
  switch (model,
          "full" = {
            switch(solver,
                   "rCMA" = {
                     if(rcpp){
                       objfun_rcma = function(par) {msdpd_me_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, w_er = w3, w_lam = w2, me_rho = me_rho_cpp, me_alp = me_alp_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th)}
                     }else{
                       objfun_rcma = function(par) {full_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp,inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th)}
                     }
                     for (i in 1:max_try){
                       cma_obj = cmaNew()
                       cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
                       cmaInit(cma_obj, dimension = 4)
                       optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma, maxDimPrint = 4)
                       if (optim_res$bestFitness < 1e-12) break
                     }
                     output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
                     output$coefficient = list(lambda1_ast = optim_res$bestX[1],
                                               lambda3_ast = optim_res$bestX[2],
                                               rho = optim_res$bestX[3],
                                               lambda2 = optim_res$bestX[4]
                     )
                     optim_par = optim_res$bestX
                   },
                   "pracma" = {
                     if(rcpp){
                       objfun_pracma = function(par) {msdpd_me_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, w_er = w3, w_lam = w2, me_rho = me_rho_cpp, me_alp = me_alp_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th, sq = F)}
                     }else{
                       objfun_pracma = function(par) {full_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp,inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, sq = F)}
                     }
                     for (i in 1:max_try){
                       init_val = runif(4, -0.1, 0.1)
                       optim_res = pracma::lsqnonlin(objfun_pracma,x0 = init_val)
                       bool_solve = optim_res$ssq < 1e-12 & const_rcma(optim_res$x)
                       if (bool_solve) break
                     }
                     output = list(solve = bool_solve, coefficient = list())
                     output$coefficient = list(lambda1_ast = optim_res$x[1],
                                               lambda3_ast = optim_res$x[2],
                                               rho = optim_res$x[3],
                                               lambda2 = optim_res$x[4]
                     )
                     optim_par = optim_res$x
                   },
                   stop("invalid solver"))
            
            output$coefficient = c(output$coefficient, full_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, mode = "beta_sigs", zero_th = zero_th))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (correction){
              if (hessian_er){
                all_se = full_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, mode = "opmd", hessian_er = hessian_er, zero_th = zero_th)
                vc_er = sqrt(diag(all_se$vc_mat))
                hes_er = sqrt(diag(solve(-all_se$hes_mat)))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["lambda3_ast_se"]] = vc_er[k+5]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_se"]] = vc_er[k+4]
                output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
                output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
                output$coefficient[["lambda1_ast_se_hes"]] = hes_er[k+3]
                output$coefficient[["lambda3_ast_se_hes"]] = hes_er[k+5]
                output$coefficient[["rho_se_hes"]] = hes_er[k+2]
                output$coefficient[["lambda2_se_hes"]] = hes_er[k+4]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                  output$hessian = all_se$hes_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                  output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }else{
                all_se = full_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, mode = "opmd", hessian_er = hessian_er, zero_th = zero_th)
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["lambda3_ast_se"]] = vc_er[k+5]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_se"]] = vc_er[k+4]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
            } else {
              all_se = full_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, w_lam = w2, sp_mode = sp_mode, mode = "opmd", hessian_er = T, zero_th = zero_th)
              hes_er = sqrt(diag(solve(-all_se$hes_mat)))
              output$coefficient[["beta_se"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se"]] = hes_er[k+1]
              output$coefficient[["lambda1_ast_se"]] = hes_er[k+3]
              output$coefficient[["lambda3_ast_se"]] = hes_er[k+5]
              output$coefficient[["rho_se"]] = hes_er[k+2]
              output$coefficient[["lambda2_se"]] = hes_er[k+4]
              if (k_o == k){
                output$hessian = all_se$hes_mat
              } else {
                output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
          },
          "slm" = {
            const_rcma_slm = function(x){
              pars = numeric(4)
              pars[1] = x[1]
              pars[3] = x[2]
              return(const_rcma(pars))
            }
            switch(solver,
                   "rCMA" = {
                     if(rcpp){
                       objfun_rcma = function(par) {
                         msdpd_me_slm_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, me_rho = me_rho_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th)
                       }
                     }else{
                       objfun_rcma = function(par) {slm_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th)}
                     }

                     for (i in 1:max_try){
                       cma_obj = cmaNew()
                       cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
                       cmaInit(cma_obj, dimension = 2)
                       optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_slm, maxDimPrint = 2)
                       if (optim_res$bestFitness < 1e-12) break
                     }
                     output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
                     output$coefficient = list(lambda1 = optim_res$bestX[1],
                                               rho = optim_res$bestX[2]
                     )
                     optim_par = optim_res$bestX
                   },
                   "pracma" = {
                     if(rcpp){
                       objfun_pracma = function(par) {msdpd_me_slm_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, me_rho = me_rho_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th, sq = F)}
                     }else{
                       objfun_pracma = function(par) {slm_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th, sq = F)}
                     }
                     for (i in 1:max_try){
                       init_val = runif(2, -0.1, 0.1)
                       optim_res = pracma::lsqnonlin(objfun_pracma,x0 = init_val)
                       bool_solve = optim_res$ssq < 1e-12 & const_rcma_slm(optim_res$x)
                       if (bool_solve) break
                     }
                     output = list(solve = bool_solve, coefficient = list())
                     output$coefficient = list(lambda1_ast = optim_res$x[1],
                                               rho = optim_res$x[2]
                     )
                     optim_par = optim_res$x
                   },
                   stop("invalid solver"))
            output$coefficient = c(output$coefficient, slm_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th, mode = "beta_sigs"))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (correction){
              if (hessian_er){
                all_se = slm_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                hes_er = sqrt(diag(solve(-all_se$hes_mat)))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
                output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
                output$coefficient[["lambda1_ast_se_hes"]] = hes_er[k+3]
                output$coefficient[["rho_se_hes"]] = hes_er[k+2]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                  output$hessian = all_se$hes_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                  output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }else{
                all_se = slm_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
            } else {
              all_se = slm_aqs_me(para = optim_par,  y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = T)
              hes_er = sqrt(diag(solve(-all_se$hes_mat)))
              output$coefficient[["beta_se"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se"]] = hes_er[k+1]
              output$coefficient[["lambda1_ast_se"]] = hes_er[k+3]
              output$coefficient[["rho_se"]] = hes_er[k+2]
              if (k_o == k){
                output$hessian = all_se$hes_mat
              } else {
                output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
          },
          "sem" = {
            const_rcma_sem = function(x){
              pars = numeric(4)
              pars[2] = x[1]
              pars[3] = x[2]
              return(const_rcma(pars))
            }
            switch(solver,
                   "rCMA" = {
                     if(rcpp){
                       objfun_rcma = function(par) {
                         msdpd_me_sem_aqs(para = par, x_ = x, y = y, y1 = y_1,  w_er = w3, me_alp = me_alp_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th)
                       }
                     }else{
                       objfun_rcma = function(par) {sem_aqs_me(para = par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th)}
                     }

                     for (i in 1:max_try){
                       cma_obj = cmaNew()
                       cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
                       cmaInit(cma_obj, dimension = 2)
                       optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_sem, maxDimPrint = 2)
                       if (optim_res$bestFitness < 1e-12) break
                     }
                     output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
                     output$coefficient = list(lambda3 = optim_res$bestX[1],
                                               rho = optim_res$bestX[2]
                     )
                     optim_par = optim_res$bestX
                   },
                   "pracma" = {
                     if(rcpp){
                       objfun_pracma = function(par) {msdpd_me_sem_aqs(para = par, x_ = x, y = y, y1 = y_1,  w_er = w3, me_alp = me_alp_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th, sq = F)}
                     }else{
                       objfun_pracma = function(par) {sem_aqs_me(para = par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th, sq = F)}
                     }
                     for (i in 1:max_try){
                       init_val = runif(2, -0.1, 0.1)
                       optim_res = pracma::lsqnonlin(objfun_pracma,x0 = init_val)
                       bool_solve = optim_res$ssq < 1e-12 & const_rcma_sem(optim_res$x)
                       if (bool_solve) break
                     }
                     output = list(solve = bool_solve, coefficient = list())
                     output$coefficient = list(lambda3_ast = optim_res$x[1],
                                               rho = optim_res$x[2]
                     )
                     optim_par = optim_res$x
                   },
                   stop("invalid solver"))
            output$coefficient = c(output$coefficient, sem_aqs_me(para = optim_par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th, mode = "beta_sigs"))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (correction){
              if (hessian_er){
                all_se = sem_aqs_me(para = optim_par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                hes_er = sqrt(diag(solve(-all_se$hes_mat)))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda3_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
                output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
                output$coefficient[["lambda3_ast_se_hes"]] = hes_er[k+3]
                output$coefficient[["rho_se_hes"]] = hes_er[k+2]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                  output$hessian = all_se$hes_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                  output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }else{
                all_se = sem_aqs_me(para = optim_par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda3_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
            } else {
              all_se = sem_aqs_me(para = optim_par, y = y, x_ = x, y1 = y_1, me_alp = me_alp, inv_c = inv_c, correction = correction, w_er = w3, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = T)
              hes_er = sqrt(diag(solve(-all_se$hes_mat)))
              output$coefficient[["beta_se"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se"]] = hes_er[k+1]
              output$coefficient[["lambda3_ast_se"]] = hes_er[k+3]
              output$coefficient[["rho_se"]] = hes_er[k+2]
              if (k_o == k){
                output$hessian = all_se$hes_mat
              } else {
                output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
          },
          "sltl" = {
            const_rcma_sltl = function(x){
              pars = numeric(4)
              pars[1] = x[1]
              pars[3] = x[2]
              pars[4] = x[3]
              return(const_rcma(pars)) 
            }
            switch(solver,
                   "rCMA" = {
                     if(rcpp){
                       objfun_rcma = function(par) {
                         msdpd_me_sltl_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, w_lam = w2, me_rho = me_rho_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th)
                       }
                     }else{
                       objfun_rcma = function(par) {sltl_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th)}
                     }
                     
                     for (i in 1:max_try){
                       cma_obj = cmaNew()
                       cmaSetPopulationSize(cma_obj, cmaGetPopulationSize(cma_obj)*cma_pop_multi)
                       cmaInit(cma_obj, dimension = 3)
                       optim_res = cmaOptimDP(cma_obj, objfun_rcma, isFeasible = const_rcma_sltl, maxDimPrint = 3)
                       if (optim_res$bestFitness < 1e-12) break
                     }
                     optim_par = optim_res$bestX
                     output = list(solve = optim_res$bestFitness < 1e-12, coefficient = list())
                     output$coefficient = list(lambda1 = optim_res$bestX[1],
                                               rho = optim_res$bestX[2],
                                               lambda2 = optim_res$bestX[3]
                     )
                   },
                   "pracma" = {
                     if(rcpp){
                       objfun_pracma = function(par) {
                         msdpd_me_sltl_aqs(para = par, x_ = x, y = y, y1 = y_1, w = w1, w_lam = w2, me_rho = me_rho_cpp, inv_c = as(inv_c,"dgCMatrix"), correction = correction, zero_th = zero_th, sq = F)
                       }
                     }else{
                       objfun_pracma = function(par) {sltl_aqs_me(para = par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, sq = F)}
                     }
                     for (i in 1:max_try){
                       init_val = runif(3, -0.1, 0.1)
                       optim_res = pracma::lsqnonlin(objfun_pracma,x0 = init_val)
                       bool_solve = optim_res$ssq < 1e-12 & const_rcma_sltl(optim_res$x)
                       if (bool_solve) break
                     }
                     output = list(solve = bool_solve, coefficient = list())
                     output$coefficient = list(lambda1_ast = optim_res$x[1],
                                               rho = optim_res$x[2],
                                               lambda2 = optim_res$x[3]
                     )
                     optim_par = optim_res$x
                   },
                   stop("invalid solver"))
            output$coefficient = c(output$coefficient, sltl_aqs_me(para = optim_par, y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, mode = "beta_sigs"))
            output$coefficient$beta = output$coefficient$beta[1:k_o]
            if (correction){
              if (hessian_er){
                all_se = sltl_aqs_me(para = optim_par,  y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                hes_er = sqrt(diag(solve(-all_se$hes_mat)))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_se"]] = vc_er[k+4]
                output$coefficient[["beta_se_hes"]] = hes_er[1:k_o]
                output$coefficient[["sigma2_se_hes"]] = hes_er[k+1]
                output$coefficient[["lambda1_ast_se_hes"]] = hes_er[k+3]
                output$coefficient[["rho_se_hes"]] = hes_er[k+2]
                output$coefficient[["lambda2_se_hes"]] = hes_er[k+4]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                  output$hessian = all_se$hes_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                  output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }else{
                all_se = sltl_aqs_me(para = optim_par,   y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = hessian_er)
                vc_er = sqrt(diag(all_se$vc_mat))
                output$coefficient[["beta_se"]] = vc_er[1:k_o]
                output$coefficient[["sigma2_se"]] = vc_er[k+1]
                output$coefficient[["lambda1_ast_se"]] = vc_er[k+3]
                output$coefficient[["rho_se"]] = vc_er[k+2]
                output$coefficient[["lambda2_se"]] = vc_er[k+4]
                if (k_o == k){
                  output$vc = all_se$vc_mat
                } else {
                  output$vc = all_se$vc_mat[-c((k_o+1):k), -c((k_o+1):k)]
                }
              }
            } else {
              all_se = sltl_aqs_me(para = optim_par,  y = y, x_ = x, w = w1, y1 = y_1, me_rho = me_rho, inv_c = inv_c, correction = correction, w_lam = w2, sp_mode = sp_mode, zero_th = zero_th, mode = "opmd",  hessian_er = T)
              hes_er = sqrt(diag(solve(-all_se$hes_mat)))
              output$coefficient[["beta_se"]] = hes_er[1:k_o]
              output$coefficient[["sigma2_se"]] = hes_er[k+1]
              output$coefficient[["lambda1_ast_se"]] = hes_er[k+3]
              output$coefficient[["rho_se"]] = hes_er[k+2]
              output$coefficient[["lambda2_se"]] = hes_er[k+4]
              if (k_o == k){
                output$hessian = all_se$hes_mat
              } else {
                output$hessian = all_se$hes_mat[-c((k_o+1):k), -c((k_o+1):k)]
              }
            }
          },
          stop("Undefined model")
  )
  output$model = model
  return(output)
}


