library(torch)
library(caret)
library(glmnet)

# ---------- Model builders ----------
delta_model <- nn_module(
  initialize = function(input_dim) {
    self$fc1 <- nn_linear(input_dim, 8)
    self$fc3 <- nn_linear(8, 1)
  },
  forward = function(x) {
    x <- self$fc1(x) %>% torch_relu()
    1 + nnf_softplus(self$fc3(x))
  }
)

phi_model <- nn_module(
  initialize = function(input_dim) {
    self$fc1 <- nn_linear(input_dim, 32)
    self$fc3 <- nn_linear(32, 1)
  },
  forward = function(x) {
    x <- self$fc1(x) %>% torch_relu()
    self$fc3(x)
  }
)

eta_model_fn <- nn_module(
  initialize = function(input_dim) {
    self$fc1 <- nn_linear(input_dim, 16)
    self$fc3 <- nn_linear(16, 1)
  },
  forward = function(x) {
    x <- torch_relu(self$fc1(x))
    self$fc3(x)
  }
)

# ---------- Objectives ----------
objective_delta <- function(R_batch, Rb_batch, fX_batch, fZ_batch,
                            delta_mod, phi_mod, Lambda_G, Lambda_H) {
  delta_pred <- delta_mod(fX_batch)
  phi_pred   <- phi_mod(fZ_batch)
  term1 <- (R_batch - Rb_batch * delta_pred) * phi_pred
  term2 <- - R_batch * (1 + Lambda_G) * phi_pred^2
  term3 <-   Lambda_H * Rb_batch * delta_pred^2
  torch_mean(term1 + term2 + term3)
}

objective_nu <- function(Rb_b, fX_b, fZ_b, kappa_b,
                         nu_mod, delta_adv_mod, Lambda_G, Lambda_H) {
  nu_pred    <- nu_mod(fZ_b)
  delta_pred <- delta_adv_mod(fX_b)
  term <- Rb_b * (kappa_b - nu_pred) * delta_pred
  pen_delta <- - Rb_b * (1 + Lambda_H) * (delta_pred^2)
  pen_nu    <-   Lambda_G * Rb_b * (nu_pred^2)
  torch_mean(term + pen_delta + pen_nu)
}

# ---------- Trainers ----------
train_delta_phi_block <- function(label,
                                  R_train_tensor, Rb_train_tensor,
                                  fX_train_tensor, fX_test_tensor,
                                  fZ_train_tensor, fZ_test_tensor,
                                  pfx, pfz,
                                  epochs, batch_size,
                                  ascent_steps, descent_steps,
                                  lr_ascent, lr_descent,
                                  Lambda_G, Lambda_H) {
  delta_mod <- delta_model(pfx)
  phi_mod   <- phi_model(pfz)
  opt_phi   <- optim_adam(phi_mod$parameters,   lr = lr_ascent)
  opt_delta <- optim_adam(delta_mod$parameters, lr = lr_descent)
  
  for (epoch in 1:epochs) {
    epoch_loss <- 0
    # num_batches <- ceiling(length(Rb_train_tensor) / batch_size)
    num_batches <- max(1, length(Rb_train_tensor) %/% batch_size)
    idx <- sample(length(R_train_tensor))
    
    for (b in 1:num_batches) {
      batch_idx <- idx[((b-1)*batch_size + 1):min(b*batch_size, length(R_train_tensor))]
      R_b  <- R_train_tensor[batch_idx, , drop = FALSE]
      Rbb  <- Rb_train_tensor[batch_idx, , drop = FALSE]
      fX_b <- fX_train_tensor[batch_idx, , drop = FALSE]
      fZ_b <- fZ_train_tensor[batch_idx, , drop = FALSE]
      
      # ascend phi
      for (ii in 1:ascent_steps) {
        opt_phi$zero_grad()
        loss_phi <- -objective_delta(R_b, Rbb, fX_b, fZ_b, delta_mod, phi_mod, Lambda_G, Lambda_H)
        loss_phi$backward()
        opt_phi$step()
      }
      # descend delta
      for (ii in 1:descent_steps) {
        opt_delta$zero_grad()
        loss_delta <-  objective_delta(R_b, Rbb, fX_b, fZ_b, delta_mod, phi_mod, Lambda_G, Lambda_H)
        loss_delta$backward()
        opt_delta$step()
      }
      epoch_loss <- epoch_loss + loss_delta$item()
    }
    if (epoch %% 10 == 0 || epoch == epochs) {
      cat(sprintf("%s, Epoch %3d/%3d, Loss: %.6f\n",
                  label, epoch, epochs, epoch_loss / num_batches))
    }
  }
  
  delta_mod$eval()
  with_no_grad({
    delta_hat_train_tensor <- delta_mod(fX_train_tensor)
    delta_hat_test_tensor  <- delta_mod(fX_test_tensor)
  })
  list(train = as_array(delta_hat_train_tensor),
       test  = as_array(delta_hat_test_tensor))
}

train_weighted_regressor <- function(label, model_fn, p_in,
                                     X_train_tensor, y_train_tensor, w_tensor, X_test_tensor,
                                     epochs, batch_size, lr = 1e-3) {
  model <- model_fn(p_in)
  opt   <- optim_adam(model$parameters, lr = lr)
  
  model$train()
  for (epoch in 1:epochs) {
    epoch_loss <- 0
    nrows = as.integer(X_train_tensor$size()[[1]])
    num_batches <- ceiling(nrows / batch_size)
    idx <- sample(nrows)
    for (b in 1:num_batches) {
      batch_idx <- idx[((b-1)*batch_size+1):min(b*batch_size, nrows)]
      X_b <- X_train_tensor[batch_idx, , drop=FALSE]
      y_b <- y_train_tensor[batch_idx, , drop=FALSE]
      w_b <- w_tensor[batch_idx, , drop=FALSE]
      opt$zero_grad()
      pred_b <- model(X_b)
      loss_b <- torch_mean(w_b * (y_b - pred_b)$pow(2))
      loss_b$backward()
      opt$step()
      epoch_loss <- epoch_loss + loss_b$item()
    }
    if (epoch %% 10 == 0 || epoch == epochs) {
      cat(sprintf("%s, Epoch %3d/%3d, Loss: %.6f\n",
                  label, epoch, epochs, epoch_loss / num_batches))
    }
  }
  
  model$eval()
  with_no_grad({
    train_pred_tensor <- model(X_train_tensor)$detach()
    test_pred_tensor  <- model(X_test_tensor)$detach()
  })
  list(model = model,
       train_tensor = train_pred_tensor,
       test_tensor  = test_pred_tensor,
       train = as_array(train_pred_tensor),
       test  = as_array(test_pred_tensor))
}

train_nu_block <- function(label, pfz,
                           fZ_train_tensor, fZ_test_tensor,
                           kappa_train_tensor,
                           fX_train_tensor, Rb_train_tensor, pfx,
                           epochs, batch_size,
                           ascent_steps, descent_steps,
                           lr_ascent, lr_descent,
                           Lambda_G, Lambda_H) {
  
  nu_mod        <- phi_model(pfz)
  delta_adv_mod <- eta_model_fn(pfx)
  opt_nu   <- optim_adam(nu_mod$parameters,        lr = lr_descent)
  opt_dadv <- optim_adam(delta_adv_mod$parameters, lr = lr_ascent)
  
  for (epoch in 1:epochs) {
    epoch_loss <- 0
    # num_batches <- ceiling(length(Rb_train_tensor) / batch_size)
    num_batches <- max(1, length(Rb_train_tensor) %/% batch_size)
    idx <- sample(length(Rb_train_tensor))
    for (b in 1:num_batches) {
      batch_idx <- idx[((b-1)*batch_size + 1):min(b*batch_size, length(Rb_train_tensor))]
      Rb_b    <- Rb_train_tensor[batch_idx, , drop = FALSE]
      fX_b    <- fX_train_tensor[batch_idx, , drop = FALSE]
      fZ_b    <- fZ_train_tensor[batch_idx, , drop = FALSE]
      kappa_b <- kappa_train_tensor[batch_idx, , drop = FALSE]
      
      # ascend delta_adv -> minimize -objective
      for (ii in 1:ascent_steps) {
        opt_dadv$zero_grad()
        loss_d <- -objective_nu(Rb_b, fX_b, fZ_b, kappa_b,
                                nu_mod, delta_adv_mod, Lambda_G, Lambda_H)
        loss_d$backward()
        opt_dadv$step()
      }
      # descend nu -> minimize objective
      for (ii in 1:descent_steps) {
        opt_nu$zero_grad()
        loss_n <-  objective_nu(Rb_b, fX_b, fZ_b, kappa_b,
                                nu_mod, delta_adv_mod, Lambda_G, Lambda_H)
        loss_n$backward()
        opt_nu$step()
      }
      epoch_loss <- epoch_loss + loss_n$item()
    }
    if (epoch %% 10 == 0 || epoch == epochs) {
      cat(sprintf("%s (nu), Epoch %3d/%3d, Obj: %.6f\n",
                  label, epoch, epochs, epoch_loss / num_batches))
    }
  }
  
  nu_mod$eval()
  with_no_grad({
    nu_hat_test_tensor <- nu_mod(fZ_test_tensor)
  })
  as_array(nu_hat_test_tensor)
}

# ---------- Higher-level pieces ----------
estimate_deltas <- function(R1_t, R2_t, R0_t, Rb_t,
                            fX_tr_t, fX_te_t,
                            fZ10_tr_t, fZ10_te_t,
                            fZ01_tr_t, fZ01_te_t,
                            fZ00_tr_t, fZ00_te_t,
                            pfx, pfz10, pfz01, pfz00,
                            epochs, batch_size,
                            ascent_steps, descent_steps,
                            lr_ascent, lr_descent,
                            Lambda_G, Lambda_H) {
  res10 <- train_delta_phi_block("R=(1,0)", R1_t, Rb_t, fX_tr_t, fX_te_t,
                                 fZ10_tr_t, fZ10_te_t, pfx, pfz10,
                                 epochs, batch_size,
                                 ascent_steps, descent_steps, lr_ascent, lr_descent,
                                 Lambda_G, Lambda_H)
  res01 <- train_delta_phi_block("R=(0,1)", R2_t, Rb_t, fX_tr_t, fX_te_t,
                                 fZ01_tr_t, fZ01_te_t, pfx, pfz01,
                                 epochs, batch_size,
                                 ascent_steps, descent_steps, lr_ascent, lr_descent,
                                 Lambda_G, Lambda_H)
  res00 <- train_delta_phi_block("R=(0,0)", R0_t, Rb_t, fX_tr_t, fX_te_t,
                                 fZ00_tr_t, fZ00_te_t, pfx, pfz00,
                                 epochs, batch_size,
                                 ascent_steps, descent_steps, lr_ascent, lr_descent,
                                 Lambda_G, Lambda_H)
  list(
    train = list(d10 = res10$train, d01 = res01$train, d00 = res00$train),
    test  = list(d10 = res10$test,  d01 = res01$test,  d00 = res00$test)
  )
}

estimate_gamma_eta_lasso <- function(MX_tr, MX_te, X_tr, X_te,
                                     Y_tr, TT_tr, Rb_tr,
                                     delta_hat_tr,
                                     alpha = 1,   # LASSO penalty (1=LASSO, 0.5=elastic net)
                                     nfolds = 5,  # cross-validation folds
                                     seed = 1234) {
  
  set.seed(seed)
  
  
  Y_tr <- as.vector(Y_tr)
  MX_tr <- as.matrix(MX_tr)
  Rb_tr <- as.vector(Rb_tr)
  TT_tr <- as.vector(TT_tr)
  delta_hat_tr <- as.vector(delta_hat_tr)
  X_tr <- as.matrix(X_tr)
  MX_te <- as.matrix(MX_te)
  X_te <- as.matrix(X_te)
  
  # -------------------------
  # Step 1. γ₀(m, x)
  # γ₀(m,x) = E[Y*ω | R=1,T=1,M,X] / E[ω | R=1,T=1,M,X]
  # -------------------------
  idx_11_tr <- which(Rb_tr == 1 & TT_tr == 1)
  X_gamma <- MX_tr[idx_11_tr, , drop = FALSE]
  Y_gamma_num <- Y_tr[idx_11_tr] * delta_hat_tr[idx_11_tr]
  Y_gamma_den <- delta_hat_tr[idx_11_tr]
  
  # Numerator model: lasso 
  fit_num <- cv.glmnet(X_gamma, Y_gamma_num, alpha = alpha, nfolds = nfolds)
  fit_den <- cv.glmnet(X_gamma, Y_gamma_den, alpha = alpha, nfolds = nfolds)
  
  # predict
  pred_num_tr <- as.vector(predict(fit_num, MX_tr, s = "lambda.min"))
  pred_den_tr <- as.vector(predict(fit_den, MX_tr, s = "lambda.min"))
  pred_num_te <- as.vector(predict(fit_num, MX_te, s = "lambda.min"))
  pred_den_te <- as.vector(predict(fit_den, MX_te, s = "lambda.min"))
  
  gamma_train <- pred_num_tr / (pred_den_tr + 1e-8)
  gamma_test  <- pred_num_te / (pred_den_te + 1e-8)
  
  # -------------------------
  # Step 2. η₀(x)
  # η₀(x) = E[γ(M,X)*ω | R=1,T=0,X] / E[ω | R=1,T=0,X]
  # -------------------------
  idx_10_tr <- which(Rb_tr == 1 & TT_tr == 0)
  X_eta <- X_tr[idx_10_tr, , drop = FALSE]
  Y_eta_num <- gamma_train[idx_10_tr] * delta_hat_tr[idx_10_tr]
  Y_eta_den <- delta_hat_tr[idx_10_tr]
  
  fit_num2 <- cv.glmnet(X_eta, Y_eta_num, alpha = alpha, nfolds = nfolds)
  fit_den2 <- cv.glmnet(X_eta, Y_eta_den, alpha = alpha, nfolds = nfolds)
  
  # 预测
  pred_num2_tr <- as.vector(predict(fit_num2, X_tr, s = "lambda.min"))
  pred_den2_tr <- as.vector(predict(fit_den2, X_tr, s = "lambda.min"))
  pred_num2_te <- as.vector(predict(fit_num2, X_te, s = "lambda.min"))
  pred_den2_te <- as.vector(predict(fit_den2, X_te, s = "lambda.min"))
  
  eta_train <- pred_num2_tr / (pred_den2_tr + 1e-8)
  eta_test  <- pred_num2_te / (pred_den2_te + 1e-8)
  
  return(list(
    gamma_train = gamma_train,
    gamma_test  = gamma_test,
    eta_train   = eta_train,
    eta_test    = eta_test
  ))
}

estimate_b1_b2_lasso <- function(MX_tr, MX_te, X_tr, X_te,
                                 T_tr, Rb_tr,
                                 delta_hat_tr,
                                 alpha = 1, nfolds = 5, seed = 1234) {
  set.seed(seed)
  
  
  T_tr <- as.vector(T_tr)
  Rb_tr <- as.vector(Rb_tr)
  delta_hat_tr <- as.vector(delta_hat_tr)
  X_tr <- as.matrix(X_tr)
  MX_tr <- as.matrix(MX_tr)
  X_te <- as.matrix(X_te)
  MX_te <- as.matrix(MX_te)
  
  # -------------------------
  # Step 1. b_1(x) = E[omega0 | R=1,x] / E[(1-T)*omega0 | R=1,x]
  # -------------------------
  idx_r1 <- which(Rb_tr == 1)
  X_b1 <- X_tr[idx_r1, , drop = FALSE]
  
  Y_num <- delta_hat_tr[idx_r1]
  Y_den <- (1 - T_tr[idx_r1]) * delta_hat_tr[idx_r1]
  
  fit_num <- cv.glmnet(X_b1, Y_num, alpha = alpha, nfolds = nfolds)
  fit_den <- cv.glmnet(X_b1, Y_den, alpha = alpha, nfolds = nfolds)
  
  pred_num_tr <- as.vector(predict(fit_num, X_tr, s = "lambda.min"))
  pred_den_tr <- as.vector(predict(fit_den, X_tr, s = "lambda.min"))
  pred_num_te <- as.vector(predict(fit_num, X_te, s = "lambda.min"))
  pred_den_te <- as.vector(predict(fit_den, X_te, s = "lambda.min"))
  
  b1_train <- pred_num_tr / (pred_den_tr + 1e-8)
  b1_test  <- pred_num_te / (pred_den_te + 1e-8)
  
  # -------------------------
  # Step 2. b_2(m,x) = E[(1-T)omega0*b1 | R=1, m,x] / E[T*omega0 | R=1, m,x]
  # -------------------------
  idx_r1 <- which(Rb_tr == 1)
  MX_b2 <- MX_tr[idx_r1, , drop = FALSE]
  Y_num2 <- (1 - T_tr[idx_r1]) * delta_hat_tr[idx_r1] * b1_train[idx_r1]
  Y_den2 <- T_tr[idx_r1] * delta_hat_tr[idx_r1]
  
  fit_num2 <- cv.glmnet(MX_b2, Y_num2, alpha = alpha, nfolds = nfolds)
  fit_den2 <- cv.glmnet(MX_b2, Y_den2, alpha = alpha, nfolds = nfolds)
  
  pred_num2_tr <- as.vector(predict(fit_num2, MX_tr, s = "lambda.min"))
  pred_den2_tr <- as.vector(predict(fit_den2, MX_tr, s = "lambda.min"))
  pred_num2_te <- as.vector(predict(fit_num2, MX_te, s = "lambda.min"))
  pred_den2_te <- as.vector(predict(fit_den2, MX_te, s = "lambda.min"))
  
  b2_train <- pred_num2_tr / (pred_den2_tr + 1e-8)
  b2_test  <- pred_num2_te / (pred_den2_te + 1e-8)
  
  # b2_train <- as.numeric(winsorize_tensor(b2_train, lower_q = 0.01))
  # b2_test <- as.numeric(winsorize_tensor(b2_test, lower_q = 0.01))
  
  return(list(
    b1_train = b1_train,
    b1_test  = b1_test,
    b2_train = b2_train,
    b2_test  = b2_test
  ))
}

estimate_nus <- function(fZ10_tr_t, fZ10_te_t,
                         fZ01_tr_t, fZ01_te_t,
                         fZ00_tr_t, fZ00_te_t,
                         kappa_tr_t,
                         fX_tr_t, Rb_tr_t, pfx,
                         pfz10, pfz01, pfz00,
                         epochs, batch_size,
                         ascent_steps, descent_steps,
                         lr_ascent, lr_descent,
                         Lambda_G, Lambda_H) {
  nu10 <- train_nu_block("nu10", pfz10, fZ10_tr_t, fZ10_te_t, kappa_tr_t,
                         fX_tr_t, Rb_tr_t, pfx,
                         epochs, batch_size,
                         ascent_steps, descent_steps,
                         lr_ascent, lr_descent, Lambda_G, Lambda_H)
  nu01 <- train_nu_block("nu01", pfz01, fZ01_tr_t, fZ01_te_t, kappa_tr_t,
                         fX_tr_t, Rb_tr_t, pfx,
                         epochs, batch_size,
                         ascent_steps, descent_steps,
                         lr_ascent, lr_descent, Lambda_G, Lambda_H)
  nu00 <- train_nu_block("nu00", pfz00, fZ00_tr_t, fZ00_te_t, kappa_tr_t,
                         fX_tr_t, Rb_tr_t, pfx,
                         epochs, batch_size,
                         ascent_steps, descent_steps,
                         lr_ascent, lr_descent, Lambda_G, Lambda_H)
  list(nu10 = nu10, nu01 = nu01, nu00 = nu00)
}

estimate_mu_pair <- function(X_tr_t, X_te_t, Y_tr_t,
                             TT_tr_t, Rb_tr_t, deltahat_tr_t,
                             px, epochs, batch_size) {
  w_mu1 <- TT_tr_t * Rb_tr_t * deltahat_tr_t
  w_mu0 <- (1 - TT_tr_t) * Rb_tr_t * deltahat_tr_t
  mu1 <- train_weighted_regressor("mu1", eta_model_fn, px, X_tr_t, Y_tr_t, w_mu1, X_te_t, epochs, batch_size, lr = 1e-3)
  mu0 <- train_weighted_regressor("mu0", eta_model_fn, px, X_tr_t, Y_tr_t, w_mu0, X_te_t, epochs, batch_size, lr = 1e-3)
  list(mu1 = mu1, mu0 = mu0)
}

# helper to build composite signals
build_kappa_train <- function(eta_tr_t, gamma_tr_t, b1_tr_t, b2_tr_t, TT_tr_t, Y_tr_t) {
  eta_tr_t + (1 - TT_tr_t) * b1_tr_t * (gamma_tr_t - eta_tr_t) +
    TT_tr_t * b2_tr_t * (Y_tr_t - gamma_tr_t)
}

build_chi_train <- function(mu_tr_t, TT_tr_t, ps_inv_tr_t, Y_tr_t, treated = TRUE) {
  if (treated) {
    mu_tr_t + TT_tr_t * ps_inv_tr_t * (Y_tr_t - mu_tr_t)
  } else {
    mu_tr_t + (1 - TT_tr_t) * ps_inv_tr_t * (Y_tr_t - mu_tr_t)
  }
}

# =========================
# main
# =========================
est_dml <- function(data.input,
                    seed = NA,
                    Lambda_G = 0.1,
                    Lambda_H = 0.00001,
                    batch_size = 128,
                    epochs = 50,
                    ascent_steps = 1,
                    descent_steps = 3,
                    learning_rate_ascent = 5e-4,
                    learning_rate_descent = 1e-3,
                    K_fold = 5,
                    eps = 1e-6) {
  
  # data.input <- datagen(1000)
  R1 <- data.input$R1
  R2 <- data.input$R2
  Rb <- R1 * R2
  R0 <- (R1 == R2)
  X  <- data.input$X.obs
  Y  <- data.input$Y
  M  <- data.input$M
  Z  <- data.input$Z
  TT <- data.input$TT
  
  Y_norm <- as.numeric(normalize(Y))
  M_norm <- as.numeric(normalize(M))
  
  n <- length(TT)
  if(is.na(seed)){seed=sample(9999999,1)}
  set.seed(seed)
  folds <- createFolds(1:n, k = K_fold, list = TRUE, returnTrain = TRUE)
  
  # storage
  delta_hat10 <- delta_hat01 <- delta_hat00 <- numeric(n)
  gamma_hat   <- eta_hat     <- numeric(n)
  b1_hat      <- b2_hat      <- numeric(n)
  nu10_hat    <- nu01_hat    <- nu00_hat <- numeric(n)
  mu1_hat     <- mu0_hat     <- numeric(n)
  Y11nu10_hat <- Y11nu01_hat <- Y11nu00_hat <- numeric(n)
  Y00nu10_hat <- Y00nu01_hat <- Y00nu00_hat <- numeric(n)
  
  for (fold in 1:K_fold) {
    cat(sprintf("=====Fold %d/%d=======\n",fold,K_fold))
    
    train_idx <- folds[[fold]]
    test_idx  <- setdiff(1:n, train_idx)
    
    # features
    fX   <- cbind(X, TT, M_norm, Y_norm)
    # fX   <- cbind(X)
    fZ10 <- cbind(Z, X[,-2],      TT, M_norm, Y_norm)
    fZ01 <- cbind(Z, X[,-1],      TT, M_norm, Y_norm)
    fZ00 <- cbind(Z, X[,-c(1,2)], TT, M_norm, Y_norm)
    MX   <- cbind(M_norm, X)
    MX1 <- cbind(M, X, X^2)
    X1 <- cbind(X, X^2)
    
    pfx <- ncol(fX); pfz10 <- ncol(fZ10); pfz01 <- ncol(fZ01); pfz00 <- ncol(fZ00)
    pmx <- ncol(MX); px <- ncol(X)
    
    # split
    fX_tr   <- fX[train_idx, , drop = FALSE];   fX_te   <- fX[test_idx, , drop = FALSE]
    fZ10_tr <- fZ10[train_idx, , drop = FALSE]; fZ10_te <- fZ10[test_idx, , drop = FALSE]
    fZ01_tr <- fZ01[train_idx, , drop = FALSE]; fZ01_te <- fZ01[test_idx, , drop = FALSE]
    fZ00_tr <- fZ00[train_idx, , drop = FALSE]; fZ00_te <- fZ00[test_idx, , drop = FALSE]
    MX_tr   <- MX[train_idx, , drop = FALSE];   MX_te   <- MX[test_idx, , drop = FALSE]
    MX1_tr   <- MX1[train_idx, , drop = FALSE];   MX1_te   <- MX1[test_idx, , drop = FALSE]
    X_tr    <- X[train_idx, , drop = FALSE];    X_te    <- X[test_idx, , drop = FALSE]
    X1_tr    <- X1[train_idx, , drop = FALSE];    X1_te    <- X1[test_idx, , drop = FALSE]
    Y_tr    <- Y[train_idx];                    Y_te    <- Y[test_idx]
    TT_tr   <- TT[train_idx];                   TT_te   <- TT[test_idx]
    R0_tr   <- R0[train_idx];                   R0_te   <- R0[test_idx]
    R1_tr   <- R1[train_idx];                   R1_te   <- R1[test_idx]
    R2_tr   <- R2[train_idx];                   R2_te   <- R2[test_idx]
    Rb_tr   <- Rb[train_idx];                   Rb_te   <- Rb[test_idx]
    
    # tensors
    fX_tr_t   <- torch_tensor(fX_tr,   dtype = torch_float()); fX_te_t   <- torch_tensor(fX_te,   dtype = torch_float())
    fZ10_tr_t <- torch_tensor(fZ10_tr, dtype = torch_float()); fZ10_te_t <- torch_tensor(fZ10_te, dtype = torch_float())
    fZ01_tr_t <- torch_tensor(fZ01_tr, dtype = torch_float()); fZ01_te_t <- torch_tensor(fZ01_te, dtype = torch_float())
    fZ00_tr_t <- torch_tensor(fZ00_tr, dtype = torch_float()); fZ00_te_t <- torch_tensor(fZ00_te, dtype = torch_float())
    MX_tr_t   <- torch_tensor(MX_tr,   dtype = torch_float()); MX_te_t   <- torch_tensor(MX_te,   dtype = torch_float())
    X_tr_t    <- torch_tensor(X_tr,    dtype = torch_float()); X_te_t    <- torch_tensor(X_te,    dtype = torch_float())
    
    Y_tr_t  <- torch_tensor(Y_tr,  dtype = torch_float())$reshape(c(-1,1))
    TT_tr_t <- torch_tensor(TT_tr, dtype = torch_float())$reshape(c(-1,1))
    R0_tr_t <- torch_tensor(R0_tr, dtype = torch_float())$reshape(c(-1,1))
    R1_tr_t <- torch_tensor(R1_tr, dtype = torch_float())$reshape(c(-1,1))
    R2_tr_t <- torch_tensor(R2_tr, dtype = torch_float())$reshape(c(-1,1))
    Rb_tr_t <- torch_tensor(Rb_tr, dtype = torch_float())$reshape(c(-1,1))
    
    # ---------- delta blocks ----------
    dres <- estimate_deltas(R1_tr_t, R2_tr_t, R0_tr_t, Rb_tr_t,
                            fX_tr_t, fX_te_t,
                            fZ10_tr_t, fZ10_te_t,
                            fZ01_tr_t, fZ01_te_t,
                            fZ00_tr_t, fZ00_te_t,
                            pfx, pfz10, pfz01, pfz00,
                            epochs, batch_size,
                            ascent_steps, descent_steps,
                            learning_rate_ascent, learning_rate_descent,
                            Lambda_G, Lambda_H)
    delta_hat10_train <- dres$train$d10
    delta_hat01_train <- dres$train$d01
    delta_hat00_train <- dres$train$d00
    delta_hat10[test_idx] <- dres$test$d10
    delta_hat01[test_idx] <- dres$test$d01
    delta_hat00[test_idx] <- dres$test$d00
    
    deltahat_train <- delta_hat00_train + delta_hat01_train + delta_hat10_train - 2
    deltahat_tr_t  <- torch_tensor(deltahat_train, dtype = torch_float())$reshape(c(-1,1))
    
    # ---------- gamma & eta ----------
    gres <- estimate_gamma_eta_lasso(MX1_tr, MX1_te, X1_tr, X1_te,
                                     Y_tr, TT_tr, Rb_tr, deltahat_train)
    gamma_tr <- gres$gamma_train
    eta_tr   <- gres$eta_train
    gamma_hat[test_idx] <- gres$gamma_test
    eta_hat[test_idx]   <- gres$eta_test
    
    # ---------- b1 & b2 ----------
    bres <- estimate_b1_b2_lasso(MX1_tr, MX1_te, X1_tr, X1_te, 
                                 TT_tr, Rb_tr, deltahat_train)
    b1_tr <- bres$b1_train; b2_tr <- bres$b2_train
    b1_hat[test_idx] <- bres$b1_test
    b2_hat[test_idx] <- bres$b2_test
    
    # ---------- kappa (train) & nu ----------
    eta_tr_t <- torch_tensor(eta_tr, dtype = torch_float())$reshape(c(-1,1))
    gamma_tr_t <- torch_tensor(gamma_tr, dtype = torch_float())$reshape(c(-1,1))
    b1_tr_t <- torch_tensor(b1_tr, dtype = torch_float())$reshape(c(-1,1))
    b2_tr_t <- torch_tensor(b2_tr, dtype = torch_float())$reshape(c(-1,1))
    kappa_tr_t <- build_kappa_train(eta_tr_t, gamma_tr_t, b1_tr_t, b2_tr_t, TT_tr_t, Y_tr_t)
    nres <- estimate_nus(fZ10_tr_t, fZ10_te_t,
                         fZ01_tr_t, fZ01_te_t,
                         fZ00_tr_t, fZ00_te_t,
                         kappa_tr_t,
                         fX_tr_t, Rb_tr_t, pfx,
                         pfz10, pfz01, pfz00,
                         epochs, batch_size,
                         ascent_steps, descent_steps,
                         learning_rate_ascent, learning_rate_descent,
                         Lambda_G, Lambda_H)
    nu10_hat[test_idx] <- nres$nu10
    nu01_hat[test_idx] <- nres$nu01
    nu00_hat[test_idx] <- nres$nu00
    
    # ---------- mu1/mu0 & chi (train) ----------
    mres <- estimate_mu_pair(X_tr_t, X_te_t, Y_tr_t,
                             TT_tr_t, Rb_tr_t, deltahat_tr_t,
                             px, epochs, batch_size)
    mu1_tr_t <- mres$mu1$train_tensor; mu0_tr_t <- mres$mu0$train_tensor
    mu1_hat[test_idx] <- mres$mu1$test
    mu0_hat[test_idx] <- mres$mu0$test
    
    ps0inv_tr_t <- b1_tr_t
    ps1inv_tr_t <- b1_tr_t / (b1_tr_t - 1 + eps)
    chi1_tr_t <- build_chi_train(mu1_tr_t, TT_tr_t, ps1inv_tr_t, Y_tr_t, treated = TRUE)
    chi0_tr_t <- build_chi_train(mu0_tr_t, TT_tr_t, ps0inv_tr_t, Y_tr_t, treated = FALSE)
    
    # ---------- representers for Y11 & Y00 ----------
    y11res <- estimate_nus(fZ10_tr_t, fZ10_te_t,
                           fZ01_tr_t, fZ01_te_t,
                           fZ00_tr_t, fZ00_te_t,
                           chi1_tr_t,
                           fX_tr_t, Rb_tr_t, pfx,
                           pfz10, pfz01, pfz00,
                           epochs, batch_size,
                           ascent_steps, descent_steps,
                           learning_rate_ascent, learning_rate_descent,
                           Lambda_G, Lambda_H)
    Y11nu10_hat[test_idx] <- y11res$nu10
    Y11nu01_hat[test_idx] <- y11res$nu01
    Y11nu00_hat[test_idx] <- y11res$nu00
    
    y00res <- estimate_nus(fZ10_tr_t, fZ10_te_t,
                           fZ01_tr_t, fZ01_te_t,
                           fZ00_tr_t, fZ00_te_t,
                           chi0_tr_t,
                           fX_tr_t, Rb_tr_t, pfx,
                           pfz10, pfz01, pfz00,
                           epochs, batch_size,
                           ascent_steps, descent_steps,
                           learning_rate_ascent, learning_rate_descent,
                           Lambda_G, Lambda_H)
    Y00nu10_hat[test_idx] <- y00res$nu10
    Y00nu01_hat[test_idx] <- y00res$nu01
    Y00nu00_hat[test_idx] <- y00res$nu00
  } # end fold
  
  # ---------- global aggregates ----------
  delta_hat  <- delta_hat00 + delta_hat01 + delta_hat10 - 2
  kappa_hat  <- eta_hat + (1 - TT) * b1_hat * (gamma_hat - eta_hat) + TT * b2_hat * (Y - gamma_hat)
  
  ps0inv <- b1_hat
  ps1inv <- b1_hat / (b1_hat - 1 + eps)
  chi1_hat <- mu1_hat + TT * ps1inv * (Y - mu1_hat)
  chi0_hat <- mu0_hat + (1 - TT) * ps0inv * (Y - mu0_hat)
  
  # ---------- influence functions (pseudo-outcomes) ----------
  # psi for theta, alpha1, alpha0
  psi_theta  <- Rb * delta_hat * kappa_hat +
    nu10_hat * (R1 - Rb * delta_hat10) +
    nu01_hat * (R2 - Rb * delta_hat01) +
    nu00_hat * (R0 - Rb * delta_hat00)
  
  psi_alpha1 <- Rb * delta_hat * chi1_hat +
    Y11nu10_hat * (R1 - Rb * delta_hat10) +
    Y11nu01_hat * (R2 - Rb * delta_hat01) +
    Y11nu00_hat * (R0 - Rb * delta_hat00)
  
  psi_alpha0 <- Rb * delta_hat * chi0_hat +
    Y00nu10_hat * (R1 - Rb * delta_hat10) +
    Y00nu01_hat * (R2 - Rb * delta_hat01) +
    Y00nu00_hat * (R0 - Rb * delta_hat00)
  
  theta_hat  <- mean(psi_theta)
  alpha1_hat <- mean(psi_alpha1)
  alpha0_hat <- mean(psi_alpha0)
  
  
  psi_ate <- psi_alpha1 - psi_alpha0
  psi_nie <- psi_alpha1 - psi_theta
  psi_nde <- psi_theta  - psi_alpha0
  
  ate_hat <- mean(psi_ate)
  nde_hat <- mean(psi_nde)
  nie_hat <- mean(psi_nie)
  
  # ---------- SD (of IF) & SE ----------
  n <- length(TT)
  
  theta_sd  <- sd(psi_theta)
  alpha1_sd <- sd(psi_alpha1)
  alpha0_sd <- sd(psi_alpha0)
  ate_sd    <- sd(psi_ate)
  nde_sd    <- sd(psi_nde)
  nie_sd    <- sd(psi_nie)
  
  theta_se  <- theta_sd  / sqrt(n)
  alpha1_se <- alpha1_sd / sqrt(n)
  alpha0_se <- alpha0_sd / sqrt(n)
  ate_se    <- ate_sd    / sqrt(n)
  nde_se    <- nde_sd    / sqrt(n)
  nie_se    <- nie_sd    / sqrt(n)
  
  
  list(theta = theta_hat,
       alpha1 = alpha1_hat,
       alpha0 = alpha0_hat,
       nde = nde_hat,
       nie = nie_hat,
       ate = ate_hat,
       theta_sd = theta_sd,
       alpha1_sd = alpha1_sd,
       alpha0_sd = alpha0_sd,
       nde_sd = nde_sd,
       nie_sd = nie_sd,
       ate_sd = ate_sd,
       aux_sd = list(delta_hat10 = delta_hat10,
                     delta_hat01 = delta_hat01,
                     delta_hat00 = delta_hat00,
                     b1_hat = b1_hat, b2_hat = b2_hat,
                     nu10_hat = nu10_hat, nu01_hat = nu01_hat, nu00_hat = nu00_hat,
                     mu1_hat = mu1_hat, mu0_hat = mu0_hat,
                     gamma_hat = gamma_hat, eta_hat = eta_hat)
  )
}

winsorize_tensor <- function(t, lower_q = NULL, upper_q = 0.99, hard_cap = NULL) {
  if (!is.null(hard_cap)) {
    return(torch_clamp(t, min = 0, max = hard_cap))
  }
  v <- as.numeric(t)
  uq <- as.numeric(quantile(v, probs = upper_q, na.rm = TRUE))
  lq <- if (is.null(lower_q)) min(v, na.rm = TRUE) else as.numeric(quantile(v, probs = lower_q, na.rm = TRUE))
  torch_clamp(t, min = lq, max = uq)
}

