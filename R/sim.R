#sub&super_diagonal_sp_matrix_generation
gen_ssd_wmat = function(n, mode = "s&s"){
  m = matrix(0, ncol = n, nrow = n)
  d_vec = rep(1,n-1)
  switch (mode,
          "s&s" = {
            diag(m[-1,]) = d_vec
            diag(m[,-1]) = d_vec
            m = m/rowSums(m)
          },
          "super" = {
            diag(m[,-1]) = d_vec
          },
          "sub" = {
            diag(m[-1,]) = d_vec
          },
          stop("unsupported mode")
  )
  return(m)
}

#make weight matrix from parameter vector
make_w =  function(n, para){
  if (length(para) != n*(n-1)) stop("para length incompatible")
  m1<- matrix(para,ncol=n)
  m2<- matrix(0,n,n)
  indx<- which(m2==0,arr.ind=TRUE)
  m2[indx[indx[,1]!=indx[,2],]]<- as.vector(m1)
  return(t(m2))
}

#random matrix with fixed number of elements per row
gen_rand_mat_1 = function(n, npr, sym = F){
  vec_elem = numeric(n*(n-1))
  for (i in 1:n) {
    vec_elem[((n-1)*(i-1)+1):((n-1)*i)][sample(1:(n-1), npr)] = 1
  }
  out_mat = make_w(n, vec_elem)
  if (sym) {
    out_mat = forceSymmetric(out_mat, "U") + forceSymmetric(out_mat, "L")
    out_mat[out_mat>0] = 1
  }
  out_mat = out_mat/rowSums(out_mat)
  return(out_mat)
}

rmixnorm = function(n, mus, sigs, mix_ratio){
  if (length(mus) != length(sigs)) stop("mus and sigs should have same length")
  index_mix = sample(seq(1, length(mus)), prob = mix_ratio, size = n,replace = T)
  samples = rnorm(n, mean = mus[index_mix], sd = sigs[index_mix])
  return(samples)
}

#simulation
sim_data2 = function(t_r, w, par_x, par_i, o_coef, beta, w_er = w, w_lam = w, err_par = list(type = "gaussian")){
  t = 100
  beta = as.matrix(beta)
  p = dim(w)[1]
  tp = t*p
  dp = diag(p)
  rho = o_coef[1]
  alp = o_coef[2]
  theta = o_coef[3]
  lambda = o_coef[4]
  sig = o_coef[5]
  irw = dp - rho*w
  iaw = dp - alp*w_er
  lw = lambda*w_lam
  k = length(beta)
  x = as.matrix(rep(1, tp))
  for (i in 1:(k-1)){
    x = cbind(x,rnorm(tp, par_x[1], par_x[2]))
  }
  x_o = x[((t - t_r )*p+ 1):tp, -1]
  tif = as.matrix(rep(1, p))
  tif = kronecker(diag(t), tif)
  x = cbind(tif, x)
  if (par_i[1] == 0){
    t_beta = rep(0,t)
  } else {
    t_beta = rnorm(t, 0, par_i[1])
  }
  beta = rbind(as.matrix(t_beta), beta)
  pif = kronecker(as.matrix(rep(1, t)), diag(p))
  x = cbind(pif, x)
  if (par_i[2] == 0){
    p_beta = rep(0,p)
  } else {
    p_beta = rnorm(p, 0, par_i[2])
  }
  beta = rbind(as.matrix(p_beta), beta)
  switch (err_par$type,
          "gaussian" = {
            err_vec = rnorm(t*p, 0, sig)
          },
          "chi_square" = {
            err_vec = rchisq(t*p, df = err_par$df)
            err_vec = scale(err_vec)
          },
          "gaussian_mixture" = {
            err_vec = rmixnorm(n = t*p, mus = err_par$mus, sigs = err_par$sigs, mix_ratio = err_par$mix_ratio)
            err_vec = scale(err_vec)
          },
          stop("invalid error type")
  )
  #first y
  mu = solve(iaw, err_vec[1:p])
  r = x[1:p,]%*%beta + mu
  y = solve(irw, r)
  #other y 
  for (i in 2:t){
    mu = solve(iaw, err_vec[(p*(i-1) + 1):(i*p)])
    r = x[(p*(i-1) + 1):(i*p), ]%*%beta + mu + theta * y[(1+p*(i-2)):(p*(i-1))] + lw%*%y[(1+p*(i-2)):(p*(i-1))]
    y[(1+p*(i-1)):(p*i)] = solve(irw, r)
  }
  ind = matrix(c(rep(seq(1,p), t_r), rep(seq(1,t_r), each = p)), ncol = 2)
  res = list("y" = cbind(ind, y[((t - t_r )*p+ 1):tp]),
             "beta" = beta,
             "o_coef" = o_coef
  )
  res[["x"]] = cbind(ind, x_o)
  return(res)
}