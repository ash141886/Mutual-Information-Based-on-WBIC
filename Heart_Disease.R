library(rstan)
library(ggplot2)
library(plotly)
library(igraph)
library(bnlearn)
library(BNSL)
df <- read.csv("Heart_Disease.csv")
selected_normal <- c("age", "chol", "trestbps", "thalach", "oldpeak")
selected_discrete <- setdiff(colnames(df), selected_normal)
discrete_data <- df[, selected_discrete]
normal_data <- df[, selected_normal]


# Min-Max normalize the data
normal_data <- as.data.frame(lapply(normal_data, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))

#normal_data <- scale(normal_data)
# Converting data frames to matrices
discrete_matrix <- as.matrix(discrete_data)
normal_matrix <- as.matrix(normal_data)
# Initializing an empty matrix for the binomial data
binomial_matrix <- matrix(0, nrow = nrow(discrete_matrix), ncol = ncol(discrete_matrix))
# Converting to binomial
for (i in 1:ncol(discrete_matrix)) {
  threshold <- median(unique(discrete_matrix[, i]))
  binomial_matrix[, i] <- ifelse(discrete_matrix[, i] > threshold, 1, 0)
}
N=nrow(df)
n1=ncol(binomial_matrix)
n2=ncol(normal_matrix)
# Combine data into a list
data_list <- list(N = N, num_binomial = n1, num_normal = n2, y1 = binomial_matrix, y2 = t(normal_matrix))

# Compile the Stan model
fit <-stan(file="Clayton_copula.stan", chains = 4,seed = 1234, data=data_list)


# Function to calculate free energy
free_energy <- function(log_likelihood) {
  return(-mean(rowSums(log_likelihood)))
}


log_lik_pairwise <- extract(fit, 'log_lik_pairwise')$log_lik_pairwise
log_lik_bernoulli <- extract(fit, 'log_lik_bernoulli')$log_lik_bernoulli
log_lik_normal <- extract(fit, 'log_lik_normal')$log_lik_normal

# Initialize lists to store free energy values
free_energy_bernoulli <- list()
free_energy_normal <- list()
free_energy_pairwise <- list()


# Calculate free energy for each Bernoulli variable
for (i in 1:dim(log_lik_bernoulli)[2]) {
  free_energy_bernoulli[paste0("x", i)] <- free_energy(log_lik_bernoulli[, i, ])
}

# Calculate free energy for each Normal variable
for (i in 1:dim(log_lik_normal)[2]) {
  free_energy_normal[paste0("y", i)] <- free_energy(log_lik_normal[, i, ])
}

# Calculate free energy for each pair
for (i in 1:dim(log_lik_pairwise)[2]) {
  pair_key <- paste0("Pair", i)
  free_energy_pairwise[pair_key] <- free_energy(log_lik_pairwise[, i, ])
}


# Initialize the J matrix for free energy
mi_bernoulli_normal <- matrix(0, nrow = length(free_energy_bernoulli), ncol = length(free_energy_normal))

# Populate the J matrix with free energy
for (i in 1:length(free_energy_bernoulli)) {
  for (j in 1:length(free_energy_normal)) {
    # Extract free energy values for the i-th Bernoulli variable, j-th Normal variable, and (i, j)-th pair
    fe_b <- free_energy_bernoulli[[paste0("x", i)]]
    fe_n <- free_energy_normal[[paste0("y", j)]]
    fe_p <- free_energy_pairwise[[paste0("Pair", (i - 1) * length(free_energy_normal) + j)]]
    
    # Calculate J_ij using the formula
    mi_bernoulli_normal[i, j] <- (1 / N) * (fe_b + fe_n - fe_p)
  }
}

# MI Matrix for Bernoulli
#mi_bernoulli_bernoulli=mi_matrix(binomial_matrix)


f_1=function(x,m){
  n=length(x); cc=array(0,dim=m)
  S=0
  for(i in 1:n){
    S=S-(1/n)*log((cc[x[i]]+0.5)/(i-1+0.5*m)); cc[x[i]]=cc[x[i]]+1
  }
  return(S) }
f_2=function(x,y,m.1,m.2){
  n=length(x); cc=array(0,dim=c(m.1,m.2))
  S=0
  for(i in 1:n){
    S=S-(1/n)*log((cc[x[i],y[i]]+0.5)/(i-1+0.5*m.1*m.2)); cc[x[i],y[i]]=cc[x[i],y[i]]+1
  }
  return(S)
}
multiplication=
  function(n,prob){
    x=runif(n); y=array(dim=n); m=length(prob)
    for(i in 1:n)for(j in 1:m){if(x[i] < prob[j]){y[i]=j; break}; x[i]=x[i]-prob[j]}
    return(y)
  }
r=100; n=1000
S=T=NULL
for(k in 1:r){
  x=multiplication(n,c(3/4,1/4)); y=multiplication(n,c(1/2,1/2))
  S=c(S, f_1(x,2)+f_1(y,2)-f_2(x,y,2,2))
}
T=NULL
for(i in 1:r){
  x=multiplication(n,c(3/4,1/4)); z=(x+multiplication(n,c(1/10,9/10)))%%2+1
  T=c(T,f_1(x,2)+f_1(z,2)-f_2(x,z,2,2))
}
MI3=function(x,y) (1/n)*(max(f_1(x,2)+f_1(y,2)-f_2(x,y,2,2),0))
n=nrow(binomial_matrix)
p=ncol(binomial_matrix)
x=matrix(nrow=n,ncol=p)
for(i in 1:p)x[,i]=binomial_matrix[[i]]
mi_bernoulli_bernoulli=matrix(nrow=p,ncol=p)
for(i in 1:(p-1))for(j in (i + 1):p) mi_bernoulli_bernoulli[i,j] = MI3(x[,i], x[,j])



# MI matrix for normal 

library(lg) 

l_lik=function(x){
  if(ncol(as.matrix(x))==1) dlg_marginal(x, eval_points=x)$f_est else dlg(lg_main(x),grid=x)$f_est
}

MI2=function(x,y){
  z=cbind(x,y)
  n=nrow(z)
  value=sum(log(l_lik(z)))/n-sum(log(l_lik(x)))/n-sum(log(l_lik(y)))/n
  return(value)
}


mat=as.matrix(normal_matrix)
p=ncol(mat); mi_normal_normal=matrix(nrow=p,ncol=p)
for(i in 1:(p-1))for(j in (i+1):p) mi_normal_normal[i,j]=MI2(mat[,i], mat[,j])

# Final MI for the mixture of Bernoulli and normal mixture

# Get dimensions
n_bernoulli <- dim(mi_bernoulli_bernoulli)[1]
n_normal <- dim(mi_normal_normal)[1]

# Initialize the combined MI matrix
n_total <- n_bernoulli + n_normal
mi_combined <- matrix(0, nrow=n_total, ncol=n_total)

# Fill in the blocks
mi_combined[1:n_bernoulli, 1:n_bernoulli] <- mi_bernoulli_bernoulli
mi_combined[(n_bernoulli + 1):n_total, (n_bernoulli + 1):n_total] <- mi_normal_normal
mi_combined[1:n_bernoulli, (n_bernoulli + 1):n_total] <- mi_bernoulli_normal
mi_combined[(n_bernoulli + 1):n_total, 1:n_bernoulli] <- t(mi_bernoulli_normal) # transpose



# Generate new row and column names
bernoulli_indices <- as.character(1:n_bernoulli)
normal_indices <- as.character((n_bernoulli + 1):(n_bernoulli + n_normal))
combined_indices <- c(bernoulli_indices, normal_indices)

# Assign the new names to the matrix
colnames(mi_combined) <- combined_indices
rownames(mi_combined) <- combined_indices

# Set lower triangular part to NA
mi_combined[lower.tri(mi_combined, diag = FALSE)] <- NA

# Set negative terms to 0
mi_combined[mi_combined <=0] <- 0

# Print or return the combined MI matrix
print(mi_combined)


#Plot with Krushkal

kruskal_mi=function(w)
{
  p=ncol(w)
  p=nrow(w)
  for(i in 1:(p-1))for(j in (i+1):p) w[j,i]=w[i,j];
  for(i in 1:p)w[i,i]=0;
  parent=array(dim=p); u=rep(0,p); v=rep(-Inf,p)
  j=1; parent[1]=-1;
  while(j<=p){
    u[j]=1;
    for(i in 1:p)if(u[i]==0&&w[j,i]>v[i]){v[i]=w[j,i]; parent[i]=j;}
    max=0; for(i in 1:p) if(u[i]==0&&v[i]>max){j=i; max=v[i];}
    if(max==0){j=1; while(j<=p&&u[j]==1)j=j+1; if(j<=p)parent[j]=-1;}
  }
  pair.1=pair.2=NULL
  for(j in 1:p)if(parent[j]!=-1){
    pair.1=c(pair.1,parent[j]);
    pair.2=c(pair.2,j)
  }
  return(cbind(pair.1,pair.2));
}




# Assuming mi_combined is your MI matrix with NAs for the lower triangle and diagonal
isolated_nodes <- rownames(mi_combined)[apply(mi_combined, 1, function(row) {
  all(row <= 0, na.rm = TRUE)
})]
edge_list <- data.frame(Source = character(), Target = character(), Weight = numeric())
for(i in 1:(nrow(mi_combined) - 1)) {
  for(j in (i + 1):ncol(mi_combined)) {
    weight <- mi_combined[i, j]
    if(!is.na(weight) && weight > 0) {
      edge_list <- rbind(edge_list, data.frame(Source = rownames(mi_combined)[i], Target = rownames(mi_combined)[j], Weight = weight))
    }
  }
}
# Create dummy edges linking each isolated node to a common dummy node
dummy_edges <- data.frame(Source = isolated_nodes, Target = "dummy_node", Weight = rep(0, length(isolated_nodes)))
# Combine real and dummy edge lists
combined_edge_list <- rbind(edge_list, dummy_edges)
# Save combined edge list
#write.csv(combined_edge_list, "mi_final_heart_disease_1.csv", row.names = FALSE)


par(mfrow=c(2,2))
stan_diag(fit)
stan_rhat(fit, bins=100)
stan_ess(fit, bins=100)