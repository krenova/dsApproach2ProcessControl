
#################################################################
# If one chooses not to run the simulation, one can load
# the simulation results below.
#################################################################
# load("./results/theory_experiment.RData")





#################################################################
# Initialization
#################################################################
library(tidyverse)
theme_set(theme_bw())

plot_scale_factor <- 2



N <- 40000                    # number of observations
threshold <- 35               # threshold for binary outcomes


a <- 30                       # DGP intercept
b1 <- 1                       # DGP slope
b2 <- 3                       # 2nd DGP slope to simulate missing variable endogeneity
b <- c(b1,b2)                 # DGP vector of slopes
eps_sig_sq <- 1               # DGP error variance


mu <- c(1,2)                  # VAR process mean
rho <- c(0.37,0.2,0.3)        # VAR error correlations for simulations 1, 2 and 3 respectively 
ups <- list(                  
  c(1.65,2.5),
  c(1,1.5),
  c(0.27,0.21)
)                             # VAR error variance terms for simulations 1, 2 and 3 respectively 

PI <- list(
  matrix(c(0.4,0.4,0.4,0.4),2,2),
  matrix(c(0.7,0.2,0.2,0.7),2,2),
  matrix(c(0.95,0.03,0.02,0.96),2,2)
)                             # Phi matrix used to model the process variance terms for simulations 1, 2 and 3 respectively



#################################################################
# Simulations with varying batch sizes
#################################################################
for (batch in c(1,2,4,8,10,16,20,40,50,100)) {

  dgp_params <- vector(mode="list",3)
  dat <- vector(mode="list",3)
  names(dgp_params) <- paste0("Simulation_",1:3)
  names(dat) <- paste0("Simulation_",1:3)
  
  for (i in 1:3) {
    dgp_params[[i]]$Rho <- matrix(c(1,rho[[i]],rho[[i]],1),nrow=2)
    dgp_params[[i]]$Sigma <- diag(sqrt(ups[[i]])) %*% dgp_params[[i]]$Rho %*% diag(sqrt(ups[[i]]))
    dgp_params[[i]]$var_mat <- matrix(solve(diag(2*2) - kronecker(PI[[i]],PI[[i]])) %*% as.vector( dgp_params[[i]]$Sigma),ncol=2)
    dgp_params[[i]]$cor_mat <- diag(1/sqrt(diag(dgp_params[[i]]$var_mat))) %*% dgp_params[[i]]$var_mat %*% diag(1/sqrt(diag(dgp_params[[i]]$var_mat)))
    dgp_params[[i]]$const <- as.vector((diag(2)-PI[[i]])%*%mu)
    
    
    set.seed(123)
    dat[[i]]$innov <- mvtnorm::rmvnorm(
      n=N,
      mean = c(0,0),
      sigma = dgp_params[[i]]$Sigma
    )
    dat[[i]]$X <- matrix(nrow=N,ncol=2)
    colnames(dat[[i]]$X) <- c("x1","x2")
    dat[[i]]$X[1,] <- dgp_params[[i]]$const + dat[[i]]$innov[1,]
    for (t in 1:(N-1)) {
      dat[[i]]$X[t+1,] <- dgp_params[[i]]$const + PI[[i]] %*% dat[[i]]$X[t,] + dat[[i]]$innov[t+1,]
    }
    
    dat[[i]]$y <- a + dat[[i]]$X %*% b + rnorm(n=N,sd=sqrt(eps_sig_sq))
  }
  
  
  
  
  #################################################################
  # Process Simulation Data
  #################################################################
  dat_prep <- vector(mode='list',3)
  dat_train <- vector(mode='list',3)
  dat_test <- vector(mode='list',3)
  results <- vector(mode='list',3)
  fit <- vector(mode='list',3)
  names(dat_prep) <- paste0("Simulation_",1:3)
  names(dat_train) <- paste0("Simulation_",1:3)
  names(dat_test) <- paste0("Simulation_",1:3)
  
  
  
  for (i in 1:3) {
    
    dat_prep[[i]] <- data.frame(
      y = dat[[i]]$y,
      z = 1*(dat[[i]]$y>threshold),
      dat[[i]]$X 
    )
    
    
    dat_train[[i]] <- dat_prep[[i]][1:(N/2),] %>%
      mutate(
        t = 1:(N/2),
        subgroup = rep(1:((N/2)/batch),each=batch)
      ) %>%
      group_by(subgroup) %>%
      mutate(
        p=mean(z),
        ybar=mean(y),
      ) %>%
      ungroup()
    
    dat_test[[i]] <- dat_prep[[i]][(N/2)+1:(N/2),] %>%
      mutate(
        t = (N/2)+1:(N/2),
        subgroup = rep(((N/2)/batch+1):(N/batch),each=batch)
      ) %>%
      group_by(subgroup) %>%
      mutate(
        p=mean(z),
        ybar=mean(y),
      ) %>%
      ungroup()
    
    
    
    
    # LM Model
    # fit[[i]] <- lm(y~x1,dat_train[[i]])
    fit[[i]] <- lm(y~x1+x2,dat_train[[i]])
    fit[[i]]$coefficients
    
    
    
    # LM Agg Optimization Model
    rmse_bar <- function(par){
      
      a_hat <- par[1]
      b1_hat <- par[2]
      b2_hat <- par[3]
      
      dat_train[[i]] %>%
        mutate(
          # ybar_hat = a_hat + b1_hat*x1,
          ybar_hat = a_hat + b1_hat*x1 + b2_hat*x1,
        ) %>%
        group_by(subgroup) %>%
        summarize(
          ybar_hat = mean(ybar_hat),
          ybar = first(ybar)
        ) %>%
        mutate(sqErr = (ybar_hat-ybar)^2) %>%
        summarize(
          RMSE_bar = mean(sqErr)
        ) %>%
        .$RMSE_bar
    } 
    
    
    results[[i]] <- optim(
      fn = rmse_bar,
      # par = c(a_hat = 5, b1_hat = 1)
      par = c(a_hat = 5, b1_hat = 1, b2_hat = 1)
    )
    
  
  }
  
  
  # Agg model to LM model on Agg performance
  # NOTE: no diff using y or y_bar
  dat_perf_metrics <- sapply(
    1:3,
    function(i)
    dat_test[[i]] %>%
      mutate(
        yhat_lm = predict(fit[[i]],dat_test[[i]]),
        # yhat_agg = results[[i]]$par[1] + results[[i]]$par[2]*x1,
        yhat_agg = results[[i]]$par[1] + results[[i]]$par[2]*x1 + results[[i]]$par[3]*x1,
        sqErr_lm_indiv = (yhat_lm  - y)^2,
        sqErr_agg_indiv = (yhat_agg  - y)^2
      ) %>%
      group_by(subgroup) %>%
      summarize(
        ybar_hat_lm = mean(yhat_lm),
        ybar_hat_agg = mean(yhat_agg),
        ybar = first(ybar),
        sqErr_lm_indiv = sum(sqErr_lm_indiv),
        sqErr_agg_indiv = sum(sqErr_agg_indiv)
      ) %>%
      mutate(
        sqErr_lm = (ybar_hat_lm-ybar)^2,
        sqErr_agg = (ybar_hat_agg-ybar)^2
      ) %>%
      summarize(
        RMSE_lm = sqrt(mean(sqErr_lm)),
        RMSE_agg = sqrt(mean(sqErr_agg)),
        RMSE_lm_indiv = sqrt(sum(sqErr_lm_indiv)/(n()*batch)),
        RMSE_agg_indiv = sqrt(sum(sqErr_agg_indiv)/(n()*batch))
      )
    ) %>%
    t() %>%
    data.frame() %>%
    mutate_at(
      .vars=vars(contains("RMSE_")),
      .funs=~unlist(.)
    ) %>%
    mutate(
      Simulation = paste0("Sim ", 1:3),
      "Relative Performance" = RMSE_agg/RMSE_lm
    ) %>%
    tibble()



  if (exists('dat_perf_metrics_all')) {
    dat_perf_metrics_all <- dat_perf_metrics_all %>%
      bind_rows(
        dat_perf_metrics %>%
          mutate(
            batch_size = batch
          )
      )
  } else {
    dat_perf_metrics_all <- dat_perf_metrics %>%
      mutate(
        batch_size = batch
      )
  }
  
  cat("current batch: ", batch," of (1,2,4,8,10,16,20,40,50,100)\n",sep="")
} #batch






#################################################################
# Overview of Results and Generate Data Frame for Analysis and
# Reporting
#################################################################
lapply(dgp_params,function(x) x$cor_mat)
cor(dat[[1]]$X);cor(dat[[2]]$X);cor(dat[[3]]$X)
lapply(dgp_params,function(x) x$var_mat)
var(dat[[1]]$X);var(dat[[2]]$X);var(dat[[3]]$X)

dat_train[[3]] %>%
  group_by(subgroup) %>%
  summarize(subgroup=first(subgroup),ybar=first(ybar)) %>%
  ggplot(aes(x=subgroup,y=ybar)) +
  geom_line()

latex2exp::TeX('x\\_1')
#-------------------------------------------
# Auto-Covariance and Auto-Correlation
#-------------------------------------------
dat_auto <- sapply(
  dat,
  function(x) acf(x$X[,1],type='covariance',plot=F,lag.max=100)$acf
) %>%
  data.frame() %>% tibble() %>%
  mutate(Variable = 'x_1') %>%
  bind_rows(
    sapply(
      dat,
      function(x) acf(x$X[,2],type='covariance',plot=F,lag.max=100)$acf
    ) %>%
      data.frame() %>% tibble() %>%
      mutate(Variable = 'x_2')
  ) %>%
  rename_at(
    .vars=vars(matches("^Sim")),
    .funs=~paste0("Sim ",gsub("(?:.*)(\\d)(?:.*)","\\1",x=set_names(.),perl=T))
  ) %>%
  mutate(
    Lag=rep(0:100,2),
    Measure = "Covariance",
    Variable = factor(Variable)
  ) %>%
  pivot_longer(
    cols=`Sim 1`:`Sim 3`,
    names_to = "Simulation",
    values_to = "Values",
    names_ptypes = factor()
  ) %>%
  bind_rows(
    sapply(
      dat,
      function(x) acf(x$X[,1],type='correlation',plot=F,lag.max=100)$acf
    ) %>%
      data.frame() %>% tibble() %>%
      mutate(Variable = 'x_1') %>%
      bind_rows(
        sapply(
          dat,
          function(x) acf(x$X[,2],type='correlation',plot=F,lag.max=100)$acf
        ) %>%
          data.frame() %>% tibble() %>%
          mutate(Variable = 'x_2')
      ) %>%
      rename_at(
        .vars=vars(matches("^Sim")),
        .funs=~paste0("Sim ",gsub("(?:.*)(\\d)(?:.*)","\\1",x=set_names(.),perl=T))
      ) %>%
      mutate(
        Lag=rep(0:100,2),
        Measure = "Correlation",
        Variable = factor(Variable)
      ) %>%
      pivot_longer(
        cols=`Sim 1`:`Sim 3`,
        names_to = "Simulation",
        values_to = "Values",
        names_ptypes = factor()
      ) 
  ) 

dat_auto %>%
  ggplot(aes(x=Lag,y=Values,color=Simulation)) +
  # geom_segment(aes(xend=Lag,yend=Covariance-Covariance)) +
  geom_point() +
  facet_grid(
    Measure~Variable,
    scales='free', 
    labeller=labeller(Variable=c(x_1=latex2exp::TeX(r"("$x_1$")"),x_2=latex2exp::TeX("$x_2$")))
  ) +
  ylab("Covariance / Correlation") +
  theme(text = element_text(size = plot_scale_factor*10))
  

ggsave(
  filename = "s3_autocor.svg",
  plot = last_plot(),
  device = "svg",
  path = "./results/charts",
  scale = plot_scale_factor,
  width = 16,
  height = 8,
  units = "cm",#c("in", "cm", "mm", "px"),
  dpi = 1000
)


################################################################
# Simulations' Performance Comparisons
################################################################
#-------------------------------------------
# Absolute Performance
#-------------------------------------------
p1 <- dat_perf_metrics_all %>%
  mutate(Simulation=factor(Simulation)) %>%
  pivot_longer(
    cols=c(RMSE_lm,RMSE_agg),
    names_to="Model",
    values_to="RMSE",
    names_ptypes=factor()
  ) %>%
  ggplot(aes(x=batch_size, y=RMSE,color=Simulation)) +
  geom_line(aes(linetype=Model),linewidth=0.7) +
  theme(text = element_text(size = plot_scale_factor*10)) +
  scale_linetype_manual( labels = c(latex2exp::TeX("$g^u$"), latex2exp::TeX("$g^a$")), values = c('solid', 'dotted') )

p1_legend <- cowplot::get_legend(p1)

#-------------------------------------------
# Performance Ratio
#-------------------------------------------
p2 <- dat_perf_metrics_all %>%
  mutate(Simulation=factor(Simulation)) %>%
  ggplot(aes(x=batch_size, y=`Relative Performance`,color=Simulation)) +
  geom_line(linewidth=0.7) +
  theme(text = element_text(size = plot_scale_factor*10))



#-------------------------------------------
# Combined Plot
#-------------------------------------------
cowplot::plot_grid(
  p1 + theme(legend.position='none'),
  p2 + theme(legend.position='none'),
  p1_legend,
  ncol=3,
  rel_widths=c(3.5, 3.5, 1)
)


ggsave(
  filename = "s3_sim_performance.svg",
  plot = last_plot(),
  device = "svg",
  path = "./results/charts",
  scale = plot_scale_factor,
  width = 20,
  height = 6.75,
  units = "cm",#c("in", "cm", "mm", "px"),
  dpi = 1000
)



# save.image("./results/theory_experiment.RData")























