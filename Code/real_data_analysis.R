################################################################################

# Goal: Run CLRM model using Gibbs sampler

rm(list=ls())

library("dplyr")
library("stringr")
library("ggpubr")
library("ggplot2")
library("tidyr")
library("tidytext")

# ----------------Set the working directory to current path -------------------#

library("rstudioapi")
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

###------ Data Read In

test.data = read.csv("./Datasets/Real Data/global_data_72216obs_38mut.csv", header= T)  
str(test.data)

# Data cleaning
test.data = test.data[, -1] 
test.data = test.data[, c(6:ncol(test.data), 1:5)] # Move mutation matrix to the front
test.data$sex <- as.factor(test.data$sex)
test.data$Patient_status <- as.factor(test.data$Patient_status)
test.data$Continent <- as.factor(test.data$Continent) 
table(test.data$sex)
table(test.data$Patient_status)
table(test.data$Continent)

# Define the custom order of patient status
status_order <- c("asymptomatic", "mild", "moderate", "recovered", "severe", "deceased")

# Convert the Patient_status variable to a factor with custom order
test.data$Patient_status <- factor(test.data$Patient_status, levels = status_order)

# Number of mutations
numMutation = ncol(test.data)-5

### Check the number of occurrences of each mutation
mutation_result <- sort(colSums(test.data[ , 1:(numMutation)]), decreasing = TRUE)
mutation_result <- t(as.data.frame(t(mutation_result)))
colnames(mutation_result) <- "Count"

### Mutation counts by each patient status
counts.by.ps <- matrix(0, nrow=numMutation, ncol=length(levels(test.data$Patient_status)))
rownames(counts.by.ps) <- colnames(test.data)[1:numMutation]
colnames(counts.by.ps) <- levels(test.data$Patient_status)
for (i in 1:numMutation){
  counts.by.ps[i,] <- tapply(test.data[,i], test.data$Patient_status, FUN = sum)
}
counts.by.ps


# Female with the presence of A23403G C3037T C14408T C241T
library(dplyr)

# Filter data for females with all four mutations present
female_with_mutations <- test.data %>%
  filter(sex == "Female") %>%
  filter(A23403G == 1, C3037T == 1, C14408T == 1, C241T == 1)

# Filter data for male with all four mutations present
male_with_mutations <- test.data %>%
  filter(sex == "Male") %>%
  filter(A23403G == 1, C3037T == 1, C14408T == 1, C241T == 1)

# Calculate ratio
ratio <- nrow(female_with_mutations) / nrow(male_with_mutations)
ratio



######################### Data Visualization #################################
# Plot patient age distribution
png("../Figures/AgeDistribution.png", width=2800, height=2500, res=500)
par(mfrow = c(1, 1))
ggplot(test.data, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(x = "Age", y = "Count") +
  theme_minimal() + 
  theme_pubclean() +
  theme(axis.text = element_text(size=15), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size=20))
dev.off()

## Plot count by patient status
png("../Figures/PatientStatusCount.png", width=2800, height=2800, res=500)
par(mfrow = c(1, 1))
# Convert the Patient_status variable to a factor with custom order
test.data$Patient_status <- factor(test.data$Patient_status, levels = status_order)
# Visualize by patient status
test.data %>% group_by(Patient_status) %>% summarise(Count=n()) %>%
  ggplot(aes(x = Patient_status, y = Count))+
  geom_bar(fill="#0073C2FF", stat = "identity", width=0.7)+
  geom_text(aes(label=Count), vjust=-0.3, size=5)+
  labs(x="Patient Status", y= "Count") +
  # coord_cartesian(ylim = c(0, 18000)) +
  theme_pubclean() +
  theme(axis.text = element_text(size=15), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size=20)) 
dev.off()

## Plot count by continents
png("../Figures/ContinentCount.png", width=2000, height=1400, res=500)
par(mfrow = c(1, 1))
# Visualize by continents
test.data %>% group_by(Continent) %>% summarise(Count=n()) %>%
  ggplot(aes(x = Continent, y = Count))+
  geom_bar(fill="#0073C2FF", stat = "identity", width=0.7)+
  geom_text(aes(label=Count), vjust=-0.3, size=6)+
  labs(x="Continent", y= "Count") +
  # coord_cartesian(ylim = c(0, 6000)) +
  theme_pubclean() +
  theme(axis.text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_text(size=20)) 
dev.off()


## Plot count by continents, patient status, and sex
png("../Figures/Continent_PatientStatus_Sex.png", width=5500, height=4500, res=500)
# Count the occurrences of patient status by sex and continent
count_data <- with(test.data, table(Patient_status, sex, Continent))
# Convert the count data to a data frame
count_df <- as.data.frame(count_data)
# Plot the count data using a stacked bar plot
ggplot(count_df, aes(x = Patient_status, y = Freq, fill = Continent)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ sex) +
  labs(x = "Patient Status", y = "Count", fill = "Continent") +
  theme_minimal() +
  theme_pubclean() +
  # coord_cartesian(ylim = c(0, 6000)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))  
dev.off()


### Draw number of top mutations by continent
# Mutation List
Mutations = names(test.data[,1:numMutation])
str(test.data)
CountsByContinent <- matrix(0, nrow = nlevels(as.factor(test.data$Continent)), ncol = length(Mutations))
rownames(CountsByContinent) <- aggregate(test.data[,Mutations[1]], by=list(Continent=test.data$Continent), FUN=sum)$Continent
colnames(CountsByContinent) <- Mutations
for (i in 1:(length(Mutations))){
  CountsByContinent[,i] <- aggregate(test.data[,Mutations[i]], by=list(Continent=test.data$Continent), FUN=sum)$x
}
CountsByContinent <- as.data.frame(CountsByContinent)
CountsByContinent 


CountList <- c()
for (j in 1:nrow(CountsByContinent)){
  CountList <- c(CountList, as.numeric(CountsByContinent[j, ]))
}

MutationCounts.df <- data.frame(
  Continent = c(rep("Africa", numMutation),
                rep("Asia", numMutation),
                rep("Europe", numMutation),
                rep("North America", numMutation),
                rep("Oceania", numMutation),
                rep("South America", numMutation)), 
  Mutation = rep(Mutations, 6),
  Count = CountList
)

png("../Figures/Mutation_Counts_by_Continents.png", width=5500, height=8500, res=500)
MutationCounts.df %>% mutate(Continent = as.factor(Continent), Mutation = reorder_within(Mutation, Count, Continent)) %>%
  ggplot(aes(x=Mutation, y=Count, fill = Continent)) +
  # geom_bar(stat='identity')+
  geom_col(show.legend = FALSE) +
  facet_wrap(~Continent,  ncol=2, scales = "free_y") +
  coord_flip() +
  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Count", y="Mutation", title = "Counts of Top 38 mutations in Each Continent")+
  # scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_pubclean() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

##############################################################################



######################### Extract Y, X, M, N, D ###############################
levels(test.data$Patient_status)
levels(test.data$Continent)

# Create binary columns for sex, patient status and continent
library(fastDummies)
test.data.dummy <- dummy_cols(test.data, 
                              select_columns = c("sex", "Patient_status"), 
                              remove_first_dummy = TRUE) # Use Asymptomatic status as first dummy variable
head(test.data.dummy)
str(test.data.dummy)

# Define the binary matrix Y_{N×M}
M <- numMutation
Y <- test.data.dummy[,1:M]

# Define the design matrix
N <- nrow(Y)
X <- as.matrix(cbind("intercept" = rep(1, N),
                     "age" = test.data.dummy$age,
                     test.data.dummy[, (ncol(test.data.dummy)-5):ncol(test.data.dummy)]))
D <- ncol(X)-1


### ---- Run Gibbs sampler

# Load required packages
library(foreach)
library(doParallel)

# Set up parallel processing
num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl=3)


# Find optimal K
source("mixture_adaptive_gibbs_sampler_single_lambda.R") # With mixture adaptive proposal on beta_k
run_Gibbs_optimal_K <- function(K) {
  num_iter <- 30000
  const_var <- 1
  scale_param <- 0.01 # Scales covariance in the fixed proposal
  R <- 100  # Initial steps for updating beta_k with fixed proposal
  p <- 1
  psi <- 0.05
  ptm <- proc.time()
  Gibbs_result <- Gibbs_CLRM(Y, X, K, num_iter, const_var, scale_param, R, p, psi)
  proc.time()-ptm
  saveRDS(Gibbs_result, file = paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, 
                                      "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var,".RDS"))
}

# Run Gibbs sampler in parallel for different values of K
K.list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)

## With parallel processing
ptm <- proc.time()
foreach(K = K.list, .combine = 'c', .verbose = TRUE) %dopar% {
  run_Gibbs_optimal_K(K)
}
proc.time()-ptm

# Stop parallel processing
stopCluster(cl)

# Run multiple chains with optimal K
source("mixture_adaptive_gibbs_sampler_single_lambda.R") # With mixture adaptive proposal on beta_k
run_Gibbs_multiple_chain <- function(num_chain) {
  K <- 4
  num_iter <- 30000
  const_var <- 1
  scale_param <- 0.01 # Scales covariance in the fixed proposal
  R <- 100  # Initial steps for updating beta_k with fixed proposal
  p <- 1
  psi <- 0.05
  ptm <- proc.time()
  Gibbs_result <- Gibbs_CLRM(Y, X, K, num_iter, const_var, scale_param, R, p, psi)
  proc.time()-ptm
  saveRDS(Gibbs_result, file = paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", 
                                      p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var, "_Chain", 2, ".RDS"))
}


# Run three chains
chain.number.list <- c(1,2,3)
foreach(chain.number = chain.number.list, .combine = 'c', .verbose = TRUE) %dopar% {
  run_Gibbs_multiple_chain(chain.number)
}
stopCluster(cl)



########################### Summarize Gibbs result #############################

##### Find optimal K

K.list <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
BIC.list <- c()
loglk.list <- c()
<<<<<<< Updated upstream
ICL.list <- c()
silhouette_scores <- numeric(length(K.list))
burnin <- 1000
=======
silhouette_scores <- numeric(length(k_range))
burnin <- 10000
>>>>>>> Stashed changes
num_iter <- 30000


for (K in K.list) {
  
  Gibbs_result <- readRDS(paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, 
                                 "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var, ".RDS"))
 
  # Drop burn-in iterations
  # selected.iterations <- burnin:num_iter
  selected.iterations <- 1000:30000

  # Draw trace plot of beta
  output_file <- paste0("./Figures/SingleLambdaRealDataBeta_TracePlot_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var, ".pdf")
  pdf(output_file, width = 8, height = 12)
  for (k in 1:K) {
    par(mfrow = c(5, 2))
    for (d in 2:(D+1)) {
      beta_values <- Gibbs_result$beta_list[k, d, selected.iterations]
      plot(
        beta_values,
        type = "l",
        xlab = "Iteration",
        ylab = bquote(beta[.(k) * "," * .(d-1)]),
        ylim = c(min(beta_values), max(beta_values))
      )
      # lines(Gibbs_result_2$beta_list[k, d, selected.iterations], col="blue")
      # lines(Gibbs_result_3$beta_list[k, d, selected.iterations], col="green")
      # abline(h = 0, col = "red", lty=2)

      density_values <- density(beta_values)
      plot(density_values,
           xlab="",
           ylab = bquote(beta[.(k) * "," * .(d-1)]),
           main = "",
           ylim = c(0, max(density_values$y)+1))
      # lines(density(Gibbs_result_2$beta_list[k, d, selected.iterations]), col="blue") # standard mcmc
      # lines(density(Gibbs_result_3$beta_list[k, d, selected.iterations]), col="green")
      # abline(v = 0, col = "red", lty=2)
    }
  }
  dev.off()
  print(paste("Saved:", output_file))
  
  # 
  # # Draw histogram of beta
  # for (k in 1:K) {
  #   # png(paste0("../Figures/RealDataBeta_", k, "_TracePlot_K", K, "_N", N, "_M", M, "_D", D, ".png"), width=3500, height=3500, res=500)
  #   par(mfrow = c(4, 2))
  #   for (d in 1:(D+1)) {
  #     # p1 <- hist(Gibbs_result_1$beta_list[k, d, selected.iterations], plot=F)
  #     # p2 <- hist(Gibbs_result_2$beta_list[k, d, selected.iterations], plot=F)
  #     # p3 <- hist(Gibbs_result_3$beta_list[k, d, selected.iterations], plot=F)
  #     # plot(p2, col=rgb(0,0,1,1/4), xlab=bquote(beta[.(k=1)*","*.(d-1)]), main="")
  #     # plot(p3, col=rgb(1,0,0,1/4), add=T)
  #     # plot(p3, col=rgb(0,1,0,1/4), add=T)
  #   }
  #   # dev.off()
  # }
  
  # Get posterior lambda
  posterior.lambda <- matrix(NA, nrow = 1, ncol = K)
  output_file <- paste0("./Figures/SingleLambdaRealDataLambda_TracePlot_K", K, "_N", N, "_M", M, "_D", D, ".pdf")
  pdf(output_file, width = 8, height = 12)
  par(mfrow = c(5, 2))
  for (k in 1:K) {
    lambda_values <- Gibbs_result$lambda_list[1, k, selected.iterations]
    plot(
      lambda_values,
      type = "l",
      xlab = "Iteration",
      ylab = bquote(lambda[.(k)])
    )
    
    density_values <- density(lambda_values)
    plot(density_values,
         xlab="",
         ylab = bquote(lambda[.(k)]),
         main = "",
         ylim = c(0, max(density_values$y)+1))
    
    posterior.lambda[1,k] <- mean(lambda_values)
  }
  dev.off()
  print(paste("Saved:", output_file))
  
  
  # # Draw trace plot of lambda
  # posterior.lambda <- matrix(NA, nrow = 1, ncol = K)
  # par(mfrow = c(3, 2))
  # # for (m in 1:M) {
  #   for (k in 1:K) {
  #     lambda_values <- Gibbs_result$lambda_list[1, k, selected.iterations]
  #     posterior.lambda[m,k] <- mean(lambda_values)
  #     plot(
  #       lambda_values,
  #       type = "l",
  #       xlab = "Iteration",
  #       ylab = bquote(lambda[.(m) * "," * .(k)])
  #     )
  #     # lines(Gibbs_result_2$lambda_list[m, k, selected.iterations], col="blue")
  #     # lines(Gibbs_result_3$lambda_list[m, k, selected.iterations], col="green")
  #     
  #     plot(density(Gibbs_result$lambda_list[1, k,selected.iterations]), 
  #          xlab="", 
  #          ylab = bquote(lambda[.(m) * "," * .(k)]), 
  #          main = "")
  #     # lines(density(Gibbs_result_2$lambda_list[m, k, selected.iterations]), col="blue") # standard mcmc
  #     # lines(density(Gibbs_result_3$lambda_list[m, k, selected.iterations]), col="green")
  #   }
  # #}
  
  
  # Look at acceptance rate of beta_k
  accept.matrix <- apply(Gibbs_result$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K){
    print(paste0("beta_", k))
    accept.rate <- sum(accept.matrix[burnin:(num_iter-1),k])/length(accept.matrix[burnin:(num_iter-1),k])
    print(accept.rate)
  }
  
  
  ### Print estimated results from Gibbs sampler
  
  # Get posterior Z and the cluster membership for each mutation
  clusterIndex <- rep(0, M)
  MutationName <- rep("", M)
  posterior.Z <- matrix(0, nrow = M, ncol = K)
  for (m in 1:M) {
    MutationName[m] <- colnames(Y)[m]
    clusterIndex[m] <- which.max(rowSums(Gibbs_result$Z_list[m, , selected.iterations]))
    posterior.Z[m, clusterIndex[m]] <- 1
  }
  
  # Get posterior beta
  posterior.beta <- matrix(NA, nrow = K, ncol = D+1)
  for (k in 1:K) {
    for (d in 1:(D+1)) {
      posterior.beta[k, d] <- mean(Gibbs_result$beta_list[k, d, selected.iterations])
    }
  }
  
  
  # Create text file
  file_name <- paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, 
                      "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var, ".txt")

  # open a file connection for writing
  file_conn <- file(file_name, open = "w")
  
  # Print cluster specific mutations, odds ratio, and 95% CI
  accept.matrix <- apply(Gibbs_result$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K) {
    writeLines(paste("Cluster ", k), file_conn)
    writeLines(paste(MutationName[which(clusterIndex==k)]), file_conn)
    summaryTable = t(apply(Gibbs_result$beta_list[k, , selected.iterations], 1, function(x)
      c(exp(mean(x)), exp(quantile(x, c(0.025, 0.975))))))
    rownames(summaryTable)=c("Intercept",
                             "age",
                             "sex_Male",
                             "Patient_status_mild",
                             "Patient_status_moderate",
                             "Patient_status_recovered ",
                             "Patient_status_severe",
                             "Patient_status_deceased")
    print(summaryTable)
    write.table(summaryTable, file_conn, row.names=TRUE, col.names=TRUE)
    
    # Write acceptance rate of beta_k
    accept.rate <- sum(accept.matrix[selected.iterations-1,k])/length(accept.matrix[selected.iterations-1,k])
    print(paste0("Acceptance rate of beta_", k, ":", accept.rate))
    writeLines(paste0("Acceptance rate of beta_", k, ":", accept.rate), file_conn)
  }
  
  # Compute BIC value
  bic <- BIC(Y, X, posterior.beta, posterior.lambda, posterior.Z)
  writeLines(paste("BIC: ", bic), file_conn)
  

  # Store BIC value
  BIC.list <- c(BIC.list, bic)
  
  # Compute loglikelihood
  loglk <- log_likelihood(Y, X, posterior.beta, posterior.lambda, posterior.Z)
  
  # # Store loglikelihood values
  # loglk.list <- c(loglk.list, loglk)
  
  # # Compute ICL value
  # icl <- ICL(Y, X, posterior.beta, posterior.lambda, posterior.Z)
  # writeLines(paste("ICL: ", icl), file_conn)
  # 
  # # Store ICL value
  # ICL.list <- c(ICL.list, icl)
  # 
  # Compute Silhouette Scores
  silhouette_scores[k - 1] <- compute_clrm_silhouette(X, posterior.Z, posterior.beta)
  
  # Compute Silhouette Scores
  silhouette_scores[k - 1] <- compute_clrm_silhouette(X, posterior.Z, posterior.beta, distance_type = "prediction")
  
  # close the file connection
  close(file_conn)

}


# # BIC for different K
# BIC.list <- c(3126641.67, 3119996.09, 3084936.30, 3085179.93, 3085581.02)

BIC.list <- c(3126676, 3085177, 3075396, 3072373, 3067441, 3067375, 3067808, 3068045, 3068174)
silhouette_scores <- c(1.0000000, 1.0000000, 1.0000000, 0.9736842, 0.9736842, 0.9736842, 0.9736842,
                       0.9736842, 0.9210526)

# Plot BIC
plot(K.list, BIC.list, ylab = "BIC", xlab="K", type="b", xaxt="n" )
axis(1, at = K.list, labels = K.list)

delta_bic <- diff(BIC.list)
plot(K.list[2:length(K.list)], delta_bic, type = "b",
     xlab = "K", ylab = "ΔBIC",
     main = "Diminishing ΔBIC for K")
abline(h = 0, col = "red")


# Plot log-likelihood
plot(K.list, loglk.list, ylab = "Log-likelihood", xlab="K", type="b", xaxt="n" )
axis(1, at = K.list, labels = K.list)

<<<<<<< Updated upstream
# # Plot ICL
# plot(K.list, ICL.list, ylab = "ICL", xlab="K", type="b", xaxt="n" )
# axis(1, at = K.list, labels = K.list)


# Visualize Silhouette Scores
plot(K.list, silhouette_scores, type = "b", xlab = "Number of Clusters (K)",
     ylab = "Average Silhouette Width", main = "")
axis(1, at = K.list, labels = K.list)

=======
# Visualize Silhouette Scores
plot(K.list, silhouette_scores[1:6], type = "b", xlab = "Number of Clusters (K)",
     ylab = "Average Silhouette Width", main = "Silhouette Method")
>>>>>>> Stashed changes

##### For multiple MCMC chains with optimal K
K = 4
burnin <- 1000
num_iter <- 30000
selected.iterations <- burnin:num_iter

Gibbs_result_1 <- readRDS(paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, 
                                 "_psi", psi, "_scaleparam", scale_param,  "_sigma^2_", const_var, "_Chain1.RDS"))
Gibbs_result_2 <- readRDS(paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, 
                                 "_psi", psi, "_scaleparam", scale_param,  "_sigma^2_", const_var, "_Chain2.RDS"))

# Draw trace plot of beta from two chains
draw_beta_two_chains <- function(Gibbs_result_1, Gibbs1_cluster_index, Gibbs_result_2, Gibbs2_cluster_index, display_index, D, selected.iterations){
  old_par <- par(mar = c(5, 5, 2, 2))  # Increase left and bottom margins
  for (d in 1:(D+1)) {
    beta_values <- Gibbs_result_1$beta_list[Gibbs1_cluster_index, d, selected.iterations]
    plot(
      beta_values,
      type = "l",
      xlab = "Iteration",
      ylab = bquote(beta[.(k = display_index) * "," * .(d-1)]),
      ylim = c(min(beta_values) - min(beta_values)*0.01, max(beta_values) + max(beta_values)*0.01),
      cex.lab = 1.2)  # Increase label size
    lines(Gibbs_result_2$beta_list[k = Gibbs2_cluster_index, d, selected.iterations], col="blue")
    
    
    plot(density(beta_values),
         xlab = "",
         ylab = bquote(beta[.(k = display_index) * "," * .(d-1)]),
         main = "",
         cex.lab = 1.2)  # Increase label size
    lines(density(Gibbs_result_2$beta_list[k = Gibbs2_cluster_index, d, selected.iterations]), col="blue")
  }
}


output_file <- paste0("./Figures/SingleLambdaRealDataBeta_TracePlot_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_scaleparam", scale_param, "_sigma^2_", const_var, ".pdf")
pdf(output_file, width = 12, height = 18)
par(mfrow = c(8, 2))
draw_beta_two_chains(Gibbs_result_1, 3, Gibbs_result_2, 3, 1, D, selected.iterations)
draw_beta_two_chains(Gibbs_result_1, 2, Gibbs_result_2, 4, 2, D, selected.iterations)
draw_beta_two_chains(Gibbs_result_1, 1, Gibbs_result_2, 2, 3, D, selected.iterations)
draw_beta_two_chains(Gibbs_result_1, 4, Gibbs_result_2, 1, 4, D, selected.iterations)
dev.off()


for (chain.number in chain.number.list) {
  
<<<<<<< Updated upstream
  Gibbs_result_1 <- readRDS(paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, 
                                   "_psi", psi, "_scaleparam", scale_param,  "_sigma^2_", const_var, "_Chain1.RDS"))
  Gibbs_result_2 <- readRDS(paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, 
                                   "_psi", psi, "_scaleparam", scale_param,  "_sigma^2_", const_var, "_Chain2.RDS"))
=======
  Gibbs_result_1 <- readRDS(paste0("../Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, ".RDS"))
  Gibbs_result_2 <- readRDS(paste0("../Results/Real Data Results/MixtureAdaptiveGibbsSamplerSingleLambdaRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_Chain1.RDS"))
>>>>>>> Stashed changes
 
  # Drop burn-in iterations
  selected.iterations <- burnin:num_iter
  
  # Draw trace plot of beta
  for (k in 1:K) {
    # png(paste0("../Figures/RealDataBeta_", k=1, "_TracePlot_K", K, "_N", N, "_M", M, "_D", D, ".png"), width=3500, height=6000, res=500)
    par(mfrow = c(3, 2))
    for (d in 2:(D+1)) {
      beta_values <- Gibbs_result_1$beta_list[k, d, selected.iterations]
      plot(
        beta_values,
        type = "l",
        xlab = "Iteration",
        ylab = bquote(beta[.(k) * "," * .(d-1)]),
        ylim = c(min(beta_values)-0.002, max(beta_values)+0.002)
      )
      lines(Gibbs_result_2$beta_list[k, d, selected.iterations], col="blue")
      # lines(Gibbs_result_3$beta_list[k=2, d, selected.iterations], col="green")
      # abline(h = 0, col = "red", lty=2)
      
      density_values <- density(Gibbs_result_1$beta_list[k, d, selected.iterations])
      plot(density_values, 
           xlab="", 
           ylab = bquote(beta[.(k) * "," * .(d-1)]), 
           main = "",
           ylim = c(0, max(density_values$y)+0.005))
      # lines(density(Gibbs_result_2$beta_list[k, d, selected.iterations]), col="blue") # standard mcmc
      # lines(density(Gibbs_result_3$beta_list[k=2, d, selected.iterations]), col="green")
      # abline(v = 0, col = "red", lty=2)
    }
    # dev.off()
  }
  
  # 
  # # Draw histogram of beta
  # for (k in 1:K) {
  #   # png(paste0("../Figures/RealDataBeta_", k, "_TracePlot_K", K, "_N", N, "_M", M, "_D", D, ".png"), width=3500, height=3500, res=500)
  #   par(mfrow = c(4, 2))
  #   for (d in 1:(D+1)) {
  #     # p1 <- hist(Gibbs_result_1$beta_list[k, d, selected.iterations], plot=F)
  #     # p2 <- hist(Gibbs_result_2$beta_list[k, d, selected.iterations], plot=F)
  #     # p3 <- hist(Gibbs_result_3$beta_list[k, d, selected.iterations], plot=F)
  #     # plot(p2, col=rgb(0,0,1,1/4), xlab=bquote(beta[.(k=1)*","*.(d-1)]), main="")
  #     # plot(p3, col=rgb(1,0,0,1/4), add=T)
  #     # plot(p3, col=rgb(0,1,0,1/4), add=T)
  #   }
  #   # dev.off()
  # }
  
  
  # Draw trace plot of lambda
  posterior.lambda <- matrix(NA, nrow = 1, ncol = K)
  par(mfrow = c(3, 2))
  # for (m in 1:M) {
    for (k in 1:K) {
      lambda_values <- Gibbs_result_2$lambda_list[1, k, selected.iterations]
      posterior.lambda[1,k] <- mean(lambda_values)
      plot(
        lambda_values,
        type = "l",
        xlab = "Iteration",
        ylab = bquote(lambda[.(k)])
      )
      # lines(Gibbs_result_2$lambda_list[m, k, selected.iterations], col="blue")
      # lines(Gibbs_result_3$lambda_list[m, k, selected.iterations], col="green")
      
      plot(density(Gibbs_result_2$lambda_list[1, k,selected.iterations]), 
           xlab="", 
           ylab = bquote(lambda[.(k)]), 
           main = "")
      # lines(density(Gibbs_result_2$lambda_list[m, k, selected.iterations]), col="blue") # standard mcmc
      # lines(density(Gibbs_result_3$lambda_list[m, k, selected.iterations]), col="green")
    }
  # }
  
  
  # Look at acceptance rate of beta_k
  accept.matrix <- apply(Gibbs_result_2$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K){
    print(paste0("beta_", k))
    accept.rate <- sum(accept.matrix[burnin:(num_iter-1),k])/length(accept.matrix[burnin:(num_iter-1),k])
    print(accept.rate)
  }
  
  
  ### Print estimated results from Gibbs sampler
  
  # Get posterior Z and the cluster membership for each mutation
  clusterIndex <- rep(0, M)
  MutationName <- rep("", M)
  posterior.Z <- matrix(0, nrow = M, ncol = K)
  for (m in 1:M) {
    MutationName[m] <- colnames(Y)[m]
    clusterIndex[m] <- which.max(rowSums(Gibbs_result_2$Z_list[m, , selected.iterations]))
    posterior.Z[m, clusterIndex[m]] <- 1
  }
  
  # Get posterior beta
  posterior.beta <- matrix(NA, nrow = K, ncol = D+1)
  for (k in 1:K) {
    for (d in 1:(D+1)) {
      posterior.beta[k, d] <- mean(Gibbs_result_2$beta_list[k, d, selected.iterations])
    }
  }
  
  
  # Create text file
  file_name <- paste0("./Results/Real Data Results/MixtureAdaptiveGibbsSamplerDataDrivenVarianceRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_psi", psi, "_Chain2.txt")
  # file_name <- paste0("../Results/Real Data Results/AdaptiveGibbsSamplerRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter, "_p", p, "_R", R, "_q", q, ".txt")
  # file_name <- paste0("../Results/Real Data Results/GibbsSamplerRealData_K", K, "_N", N, "_M", M, "_D", D, "_iter", num_iter,".txt")
  
  # open a file connection for writing
  file_conn <- file(file_name, open = "w")
  
  # Print cluster specific mutations, odds ratio, and 95% CI
  accept.matrix <- apply(Gibbs_result_2$accept_beta_list, 2, c) # convert array to matrix
  for (k in 1:K) {
    writeLines(paste("Cluster ", k), file_conn)
    writeLines(paste(MutationName[which(clusterIndex==k)]), file_conn)
    summaryTable = t(apply(Gibbs_result_2$beta_list[k, , selected.iterations], 1, function(x)
      c(exp(mean(x)), exp(quantile(x, c(0.025, 0.975))))))
    rownames(summaryTable)=c("Intercept",
                             "age",
                             "sex_Male",
                             "Patient_status_mild",
                             "Patient_status_moderate",
                             "Patient_status_recovered ",
                             "Patient_status_severe",
                             "Patient_status_deceased")
    print(summaryTable)
    write.table(summaryTable, file_conn, row.names=TRUE, col.names=TRUE)
    
    # Write acceptance rate of beta_k
    accept.rate <- sum(accept.matrix[selected.iterations-1,k])/length(accept.matrix[selected.iterations-1,k])
    print(paste0("Acceptance rate of beta_", k, ":", accept.rate))
    writeLines(paste0("Acceptance rate of beta_", k, ":", accept.rate), file_conn)
  }
  
  # close the file connection
  close(file_conn)
  
}
