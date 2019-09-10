source("/Users/minguyen/Documents/Documents/College/MSc/Thesis/R Code/GPD WQE/wqe.R")

library('fitdistrplus')
library('evd')


############################## GG plot ######################################

th <- theme(legend.position = 'bottom',
            axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
            axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
            axis.title.x = element_text(color = "grey20", size = 14,face = "plain"),
            axis.title.y = element_text(color = "grey20", size = 14,face = "plain"),
            legend.text=element_text(size=14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())


############################## Source Data ###################################

data <- data.frame(read.csv("/Users/minguyen/Documents/Documents/College/MSc/Thesis/R Code/GPD WQE/SampleData.CSV",header = TRUE))
data <- data$x

# data <- as.numeric(as.character(data$price))
# data <- data - lag(data,1)
# data <- data[data  > 0 & !is.na(data)]
# data <- sort(data)
# data <-data*100

# data <- sort( rweibull(8456, 0.2,300) )
# true_par <- c(0.2,300)
# true_q <- qweibull

############################# Preliminary Information ############################

n <- length(data)
period <- 6.5
stoptype <- 'grad'
tol <- 1e-15
rad0 <- 500
eta <- 0.02
gamma = 1e-2
theta0 = c(0,150,1.1)


############################# Log-Normal ##################################

muln <- sum(log(data))/n
sigmaln <- sqrt(sum( (log(data) - muln)^2)/n)
fnvar(data,plnorm,rlnorm,period,0,100000,c(muln,sigmaln),0.999)



################################# Pickands ###################################

pickands <- xi_pickands(data)

median(pickands$xi)

scale_gof <- 5

p <- ggplot(pickands, aes(x = k))
p <- p + geom_line(aes(y = xi*scale_gof, colour = "xi"),size = 1.5) 
 geom_line(aes(y = gof*scale_gof, colour = "gof"),size = 1.5) +
 geom_hline(yintercept=0.068*scale_gof,linetype="dotted", color = "springgreen4",size = 1.5) +
 annotate("text", x = 4000, y = 0.05*scale_gof, label = "0.068",parse = TRUE,color = 'springgreen4', size = 6) +
 annotate("text", x = 3000, y = 0.15*scale_gof, label = "GoF",parse = TRUE,color = 'springgreen4', size = 6)

p <- p + geom_line(aes(y = sigma, colour = "sigma"),size = 1.5) 

p <- p + scale_y_continuous( 
                             limits = c(-2,5),
                              sec.axis = sec_axis(~./scale_gof,name = expression(xi))
                             )

p <- p + labs(y = expression(sigma),
              x = "k",
              colour = "Parameter")
p <- p + scale_colour_manual(values = c(
                                        "red", 
                                        "blue"))

p <- p + th + theme(legend.position = 'none',
                    axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.title.x = element_text(color = "grey20", size = 20,face = "plain"),
                    axis.title.y.left = element_text(color = "red", size = 20,face = "plain"),
                    axis.title.y.right = element_text(color = "blue", size = 20,face = "plain"),
                    legend.text=element_text(size=20),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())
p



################################### Hill ###################################

hill <- xi_hill(data)

median(hill$xi)

# 95% Confidence Interval
ggplot(hill, aes(x = k)) + geom_line(aes( y = xi,colour = 'hill'), size = 1.5) + 
  geom_line( aes(x = k,y = xi - conf,colour = 'hill_l'),linetype = 'dotted', size = 1.5) + 
  geom_line( aes(x = k,y = xi + conf,colour = 'hill_u'),linetype = 'dotted', size = 1.5) + 
  geom_line( pickands, mapping = aes(x = k, y = xi, colour = 'pick'), size = 1.5) +
  geom_line( pickands, mapping = aes(x = k, y = xi - conf, colour = 'pick_l'),linetype = 'dotted', size = 1.5) +
  geom_line( pickands, mapping = aes(x = k, y = xi + conf, colour = 'pick_u'),linetype = 'dotted', size = 1.5) +
  scale_y_continuous( limits = c(-1,4) ) +
  labs(y = expression(xi), x = "k") +
  scale_colour_manual(values = c("red","red","red","blue","blue","blue")) +
  th + theme(legend.position = 'none',
               axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
               axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
               axis.title.x = element_text(color = "grey20", size = 20,face = "plain"),
               axis.title.y.left = element_text(color = "grey20", size = 20,face = "plain"),
               axis.title.y.right = element_text(color = "grey20", size = 20,face = "plain"),
               legend.text=element_text(size=20),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank())

# Compare with ODR
ggplot(hill, aes(x = k)) + geom_line(aes( y = xi,colour = 'Hill'), size = 1.5) + 
  geom_line( pickands, mapping = aes(x = k, y = xi, colour = 'Pickands'), size = 1.5) +
  geom_hline(aes(yintercept=0.228,color = "ODR"),linetype="dashed",  size = 1.5) +
  geom_hline(aes(yintercept=0.15,color = NA),linetype="dashed",  size = 1.5) +
  scale_y_continuous( limits = c(-0.5,2) ) +
  labs(y = expression(xi), x = "k") +
  scale_colour_manual(values = c("red","springgreen4","blue")) +
  th + theme(legend.position = 'bottom',
             axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
             axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
             axis.title.x = element_text(color = "grey20", size = 20,face = "plain"),
             axis.title.y.left = element_text(color = "grey20", size = 20,face = "plain"),
             axis.title.y.right = element_text(color = "grey20", size = 20,face = "plain"),
             legend.text=element_text(size=20),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank())






################################# Fit A ###################################


# Stage 1 

fittype <- 'A'

cut <- 0

s1_A_results_GN <- stage1(data,theta0,algo = 'GN',period = period, mquants = NULL, cut_pct = cut, rad0 = rad0, 
                          fittype = fittype, stoptype = stoptype, eta = eta, maxIter = 20000, 
                          tol = tol,printIter = TRUE)

pop_bias(  cut, s1_A_results_GN , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s1_A_results_GN , 0.999 , qgpd ,  true_q , true_par )   

s1_A_results_LM <- stage1(data,theta0,algo = 'LM',period = period, mquants = NULL, cut_pct = cut, rad0 = rad0, 
                          fittype = fittype, stoptype = stoptype, eta = eta, maxIter = 100000, 
                          tol = tol,printIter = TRUE)

s1_A_rhohat <- s1_A_results_LM$fullresults$rhohat

pop_bias(  cut, s1_A_results_LM , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s1_A_results_LM , 0.999 , qgpd ,  true_q , true_par )   


ggplot( s1_A_results_LM$F_data_frame, aes(x = X))  + geom_line(aes(y = F_n,colour = 'Empirical'), size = 1.5 ) +
  geom_line(aes(x = X,y = F_x,colour = 'LM Fitted'), size = 1.5 ) +
  geom_line(data = s1_A_results_GN$F_data_frame,  aes(x = X,y = F_x,colour = 'GN Fitted'), size = 1.5 ) +
  xlab('X') +
  ylab('F(X)') +
  scale_x_continuous(trans = 'log10')  +
  scale_colour_manual('', values = c('Empirical' = 'red','LM Fitted' = 'blue', 'GN Fitted' = 'springgreen4')) +
  th


# Stage 2

v_A <- 15

# NULL w2 returns equal weights (vector of 1's is default)
w2 <- NULL

s2_A_results <- stage2(data,theta0, v = v_A, epsilon0 = NULL, gamma = NULL, rhohat = NULL, w2 = w2, 
                       cut_pct = 0,period = period, fittype = fittype, mquants = NULL,
                       rad0 = 500, stoptype = stoptype,tol = tol, maxIter = 30000, 
                       eta = eta, printIter = TRUE)

pop_bias(  cut, s2_A_results , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s2_A_results , 0.999 , qgpd ,  true_q , true_par )   

ggplot( s2_A_results$F_data_frame, aes(x = X,y = F_n,colour = 'Empirical')) + geom_line( size = 1.5 ) +
  geom_line(aes(x = X,y = F_x,colour = 'Fitted'), size = 1.5 ) +
  xlab('X') +
  ylab('F(X)') +
  scale_x_continuous(trans = 'log10')  +
  scale_colour_manual('', values = c('Empirical' = 'red','Fitted' = 'blue')) +
  th +
  theme( legend.position = 'bottom')



# Truncate the tail

tails <- seq(0,0.95,0.05)

nIter_vector <- seq(10000,30000,length.out = length(tails))

truncated_s2_A_results <- truncate_s2_tail(data,theta0,tail_cutoffs = tails, 
                                           mquants = NULL, v = v_A,
                                           fittype = fittype, period = period,rad0 = rad0, w2 = w2,
                                           stoptype = stoptype, tol = tol, eta = eta, 
                                           hold_weight_ratio = TRUE,
                                           nIter_vector = nIter_vector,                                            
                                           print_cut_Iter = TRUE)

scale_gof <- 10

p_B <- ggplot(truncated_s2_A_results, aes(x = cuts))

p_B <- p_B + geom_line(aes(y = sambias_t, colour = "mcvar_t"),linetype="dashed",size = 1.5) + 
  geom_line(aes(y = sambias_f, colour = "mcvar_f"),size = 1.5)

p_B <- p_B + 
  geom_line(aes(y = gof_t*scale_gof, colour = "gof_t"),size = 1.5,linetype="dashed") +
  geom_line(aes(y = gof_f*scale_gof, colour = "gof_f"),size = 1.5) +
  geom_hline(yintercept=0.068*scale_gof,linetype="dotted", color = "blue",size = 1.5) +
  annotate("text", x = 0.5, y = 0.075*scale_gof, label = "0.068",parse = TRUE,color = 'blue', size = 6) 

p_B <- p_B + scale_y_continuous(limits = c(0,2),
                            sec.axis = sec_axis(~./scale_gof, name = "GoF Stat."))

p_B <- p_B + scale_colour_manual(values = c("blue","blue","red", "red"))

p_B <- p_B + labs(y = "MC VaR (mln GBP)",
              x = "Cutoff (% of data set)",
              colour = "Parameter")

p_B <- p_B + th + theme(legend.position = 'none',
                    axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.title.x = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.title.y.left = element_text(color = "red", size = 14,face = "plain"),
                    axis.title.y.right = element_text(color = "blue", size = 14,face = "plain"),
                    legend.text=element_text(size=14),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())

p_B












################################# Fit B ###################################

# Stage 1 

fittype <- 'B'

cut_B <- 0

s1_B_results_GN <- stage1(data,theta0,algo = 'GN',period = period, mquants = NULL, cut_pct = cut_B, rad0 = rad0, 
                          fittype = fittype, stoptype = stoptype, eta = eta, maxIter = 20000, 
                          tol = tol,printIter = TRUE)

pop_bias(  cut, s1_B_results_GN , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s1_B_results_GN , 0.999 , qgpd ,  true_q , true_par )   


s1_B_results_LM <- stage1(data,theta0,algo = 'LM',period = period, mquants = NULL, cut_pct = cut_B, rad0 = rad0, 
                          fittype = fittype, stoptype = stoptype, eta = eta, maxIter = 30000, 
                          tol = tol,printIter = TRUE)

pop_bias(  cut, s1_B_results_LM , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s1_B_results_LM , 0.999 , qgpd ,  true_q , true_par )   

ggplot( s1_B_results_GN$F_data_frame, aes(x = X,y = F_n,colour = 'Empirical') ) + geom_line() +
  geom_line(aes(x = X,y = F_x,colour = 'GN Fitted') ) +
  geom_line(data = s1_B_results_LM$F_data_frame,  aes(x = X,y = F_x,colour = 'LM Fitted') ) +
  geom_vline(xintercept=s1_B_results_GN$par[1],linetype="dashed", color = "blue") +
  annotate("text", x = s1_B_results_GN$par[1] , y = 0.5, label = "mu_GN",parse = TRUE,color = 'blue') +
  geom_vline(xintercept=s1_B_results_LM$par[1],linetype="dashed", color = "springgreen4") +
  annotate("text", x = s1_B_results_LM$par[1] , y = 0.5, label = "mu_LM",parse = TRUE,color = 'springgreen4') +
  xlab('X') +
  ylab('F(X)') +
  scale_x_continuous(trans = 'log10')  +
  scale_colour_manual('', values = c('Empirical' = 'red','GN Fitted' = 'blue', 'LM Fitted' = 'springgreen4')) +
  th


s1_B_rhohat <- s1_B_results_LM$fullresults$rhohat



# Stage 2

v_B <- 70

w2 <- NULL

s2_B_results <- stage2(data,s1_B_results_LM$par, v = v_B, epsilon0 = NULL, gamma = NULL, rhohat = NULL, w2 = w2, 
                       cut_pct = 0,period = period, fittype = fittype, mquants = NULL,
                       rad0 = 500, stoptype = stoptype,tol = tol, maxIter = 100000, 
                       eta = eta, printIter = TRUE)

pop_bias(  cut, s2_B_results , 0.95 , qgpd ,  true_q , true_par )   
pop_bias(  cut, s2_B_results , 0.999 , qgpd ,  true_q , true_par )  

ggplot( s2_B_results$F_data_frame, aes(x = X,y = F_n,colour = 'Empirical') ) + geom_line() +
  geom_line(aes(x = X,y = F_x,colour = 'Fitted') ) +
  geom_vline(xintercept=s2_B_results$par[1],linetype="dashed", color = "springgreen4") +
  annotate("text", x = s2_B_results$par[1] - 30, y = 0.5, label = "mu",parse = TRUE,color = 'springgreen4') +
  xlab('X') +
  ylab('F(X)') +
  scale_x_continuous(trans = 'log10')  +
  scale_colour_manual('', values = c('Empirical' = 'red','Fitted' = 'blue')) +
  th


# Truncate the tail

tails <- seq(0,0.98,0.02)

nIter_vector <- seq(50000,60000,length.out = length(tails))

truncated_s2_B_results <- truncate_s2_tail(data,theta0,tail_cutoffs = tails, 
                                           mquants = mquants, v = v_B,gamma = NULL, rhohat = NULL,
                                           fittype = fittype, period = period,rad0 = rad0, w2 = NULL,
                                           stoptype = stoptype, tol = tol, eta = eta, 
                                           hold_weight_ratio = TRUE,nIter_vector = nIter_vector, 
                                           print_cut_Iter = TRUE)

scale_gof <- 500

p_B <- ggplot(truncated_s2_B_results, aes(x = cuts))

p_B <- p_B + geom_line(aes(y = mcvar_t, colour = "mcvar_t"),linetype="dashed",size = 1.5) + 
  geom_line(aes(y = mcvar_f, colour = "mcvar_f"),size = 1.5)

p_B <- p_B + 
  geom_line(aes(y = gof_t*scale_gof, colour = "gof_t"),size = 1.5,linetype="dashed") +
  geom_line(aes(y = gof_f*scale_gof, colour = "gof_f"),size = 1.5) +
  geom_hline(yintercept=0.068*scale_gof,linetype="dotted", color = "blue",size = 1.5) +
  annotate("text", x = 0.5, y = 0.075*scale_gof, label = "0.068",parse = TRUE,color = 'blue', size = 6) 

p_B <- p_B + scale_y_continuous(limits = c(0,250),
                            sec.axis = sec_axis(~./scale_gof, name = "GoF Stat."))

p_B <- p_B + scale_colour_manual(values = c("blue","blue","red", "red"))

p_B <- p_B + labs(y = "MC VaR (mln GBP)",
              x = "Cutoff (% of data set)",
              colour = "Parameter")

p_B <- p_B + th + theme(legend.position = 'none',
                    axis.text.x = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.text.y = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.title.x = element_text(color = "grey20", size = 14,face = "plain"),
                    axis.title.y.left = element_text(color = "red", size = 14,face = "plain"),
                    axis.title.y.right = element_text(color = "blue", size = 14,face = "plain"),
                    legend.text=element_text(size=14),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())

p_B






################################# Tail Index ###############################

xi_s2_A_low <- get_ODR_xi_range(data,theta0 = theta0,epsilon0 = NULL,tail_cutoffs = tails, 
                             mquants = NULL,v = 1.5,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype = 'A', period = 6.5,rad0 = 500,
                             stoptype = 'grad',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = nIter_vector, print_cut_Iter = TRUE)


xi_s2_A_up <- get_ODR_xi_range(data,theta0 = theta0,epsilon0 = NULL,tail_cutoffs = tails, 
                             mquants = NULL,v = 300,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype =  'A', period = 6.5,rad0 = 500,
                             stoptype = 'grad',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = nIter_vector, print_cut_Iter = TRUE)

xi_s2_B_low <- get_ODR_xi_range(data,theta0 = theta0,epsilon0 = NULL,tail_cutoffs = tails, 
                             mquants = NULL,v = 30,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype =  'B', period = 6.5,rad0 = 500,
                             stoptype = 'grad',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = nIter_vector, print_cut_Iter = TRUE)

xi_s2_B_up <- get_ODR_xi_range(data,theta0 = theta0,epsilon0 = NULL,tail_cutoffs = tails, 
                             mquants = NULL,v = 20,gamma = NULL, rhohat = NULL,
                             w2 = NULL,fittype =  'B', period = 6.5,rad0 = 500,
                             stoptype = 'grad',tol = 1e-10, eta = 0.2,hold_weight_ratio = TRUE,
                             nIter_vector = nIter_vector, print_cut_Iter = TRUE)

pickands <- xi_pickands(data)

hill <- xi_hill(data)

ggplot( xi_s2_A_low, aes(x = rev(cuts),y = xi,colour = 'A Lower')) + geom_line( size = 1.5, linetype = 'twodash' ) +
#  geom_line(xi_s2_A_up, mapping = aes(x = rev(cuts), y = xi,colour = 'A Upper'), size = 1.5 ) +
  geom_line(xi_s2_B_low,mapping =  aes(x = rev(cuts), y = xi,colour = 'B Lower'), size = 1.5,linetype = 'dotted' ) +
#  geom_line(xi_s2_B_up,mapping =  aes(x = rev(cuts), y = xi,colour = 'B Upper'), size = 1.5,linetype = 'longdash' ) +
  geom_line(pickands,mapping =  aes(x = k/(length(pickands$xi)*4), y = xi,colour = 'Pickands'), size = 1.5 ) +
  geom_line(hill,mapping =  aes(x = k/(length(hill$xi)), y = xi,colour = 'Hill'), size = 1.5 ) +
  scale_y_continuous(limits = c(0,2.3) ) +
#  scale_colour_manual( labels = c('A Lower','A Upper','B Lower','B Upper','Hill','Pickands'),
 #                      values = c('blue','red','springgreen4','#FFFF00','purple','#00CCCC')) +
  xlab("Cutoff (% of data set)") +
  ylab(expression(xi) ) +
  th 










