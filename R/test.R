library(Rcpp)
library(spdep)
library(Matrix)
#depend on the solver used
# library(rCMA)
library(pracma)

w_5x10_rook = nb2mat(cell2nb(5,10))

#full===============================
sourceCpp("./cpp/msdpd_me.cpp")

full_sim1 = sim_data2(t_r = 5,
                      w = w_5x10_rook, 
                      par_x = c(1,1),
                      par_i = c(0,1), 
                      o_coef =  c(0.2, 0.2, 0.3, 0.2, 1),
                      beta = c(0,2))

full_res1 = msdpd_me(y = full_sim1$y, 
                     x = full_sim1$x, 
                     q = c(16,16), 
                     w1 = as(w_5x10_rook, "sparseMatrix"), 
                     correction = T, 
                     no_tf = F, 
                     model = "full", 
                     rcpp = T, 
                     sp_mode = T, 
                     zero_th = 1e-12, 
                     solver = "pracma")

#slm=================================
sourceCpp("./cpp/msdpd_me_slm.cpp")
slm_sim1 = sim_data2(t_r = 5,
                     w = w_5x10_rook, 
                     par_x = c(1,1),
                     par_i = c(0,1), 
                     o_coef =  c(0.2, 0, 0.1, 0, 1),
                     beta = c(0,2))

slm_res1 = msdpd_me(y = slm_sim1$y, 
                    x = slm_sim1$x, 
                    q = c(16,16), 
                    w1 = as(w_5x10_rook, "sparseMatrix"), 
                    correction = T, 
                    no_tf = F, 
                    model = "slm", 
                    rcpp = T, 
                    sp_mode = T, 
                    solver = "pracma")

#sem===================================
sourceCpp("./cpp/msdpd_me_sem.cpp")
sem_sim1 = sim_data2(t_r = 5,
                     w = w_5x10_rook, 
                     par_x = c(1,1),
                     par_i = c(0,1), 
                     o_coef =  c(0, 0.6, 0.1, 0, 1),
                     beta = c(0,2))


sem_res1 = msdpd_me(y = sem_sim1$y, 
                    x = sem_sim1$x, 
                    q = c(16,16), 
                    w1 = as(w_5x10_rook, "sparseMatrix"), 
                    correction = T, 
                    no_tf = F, 
                    model = "sem", 
                    rcpp = T, 
                    sp_mode = T, 
                    solver = "pracma")

#sltl=================================
sourceCpp("./cpp/msdpd_me_sltl.cpp")
sltl_sim1 = sim_data2(t_r = 5,
                      w = w_5x10_rook, 
                      par_x = c(1,1),
                      par_i = c(0,1), 
                      o_coef =  c(0.2, 0, 0.2, 0.2, 1),
                      beta = c(0,2))

sltl_res1 = msdpd_me(y = sltl_sim1$y, 
                     x = sltl_sim1$x, 
                     q = c(16,16), 
                     w1 = as(w_5x10_rook, "sparseMatrix"), 
                     correction = T, 
                     no_tf = F, 
                     model = "sltl",
                     rcpp = T, 
                     sp_mode = T,
                     solver = "pracma")


