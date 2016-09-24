# Data S1

# Optimizing species detection and sampling schemes

# Define function to simulate data 

data.fn <- function(R = h, T = n,
                    psi = g,
                    p = f){
  
  # Create empty data frames
  y <- array(NA, dim = c(R, T))		
  
  z <- rbinom(R, 1, psi)
  
  prob <- z * p		
  
  for(j in 1:T){
    
    y[,j] <- rbinom(R, 1, prob)
    
  }
  
  return(list(R = R, T = T,
              psi = psi, p = p,
              z = z, y = y))
}

# Define model in BUGS language

sink("model.txt")
cat("
    model{
    
    # Priors
    
    psi ~ dunif(0,1)
    p ~ dunif(0,1)
    
    # Ecological model for the true abundance	
    
    for(i in 1:R){					
    
    z[i] ~ dbern(psi)											
    
    } 
    
    # Observational model
    
    for(i in 1:R){
    
    for(j in 1:T){		
    
    p.eff[i,j] <- z[i] * p 					
    
    y[i,j] ~ dbern(p.eff[i,j]) 			# Detection
    
    y.new[i,j] ~ dbern(p.eff[i,j])		# Simulated data			
    
    } #js									
    
    } #is
    
    }
    ", fill = TRUE)
sink()

# Monitor Parameters
params <- c("psi", "p")

# MCMC settings
ni <- 10000
nt <- 10
nb <- 100
nc <- 3

# Create vectors with parameters for simulated data
Tpsi <- seq(from = 0.1, to = 1, by = 0.1)

Tp <-  seq(from = 0.1, to = 1, by = 0.1)

sit <- seq(from = 5, to = 200, by = 20)

surv <- seq(from = 1, to = 10, by = 2)


nsheet <- length(Tpsi) * length(Tp) * length(sit) * length(surv)
nrow <- 25
ncol <- 8

# Create empty matrix to store data
store <- array(NA, dim = c(nrow, ncol, nsheet))

colnames(store) <- c("True_psi", "True_p", "Nsurveys", "Nsites", "Psi_Mean", "SD", "ylo", "yhi")

Nsheet <- 1


for(g in Tpsi){  			# Psi   
  
  for(f in Tp){				# p     
    
    for(h in sit){			# site 
      
      for(n in surv){			# Number of surveys  (4) 
        
        for(q in 1:25){		# number of simulated data sets
          
          # Simulate the data
          sodata <- data.fn(R = h, T = n,
                            psi = g,
                            p = f)	
          
          # Bundle the data
          win.data <- list(y = sodata$y, R = nrow(sodata$y), T = ncol(sodata$y))
          
          # Create initial values
          zst <- apply(win.data$y, 1, max, na.rm = T) 
          
          inits <- function() {list(z = zst)}
          
          # Run the model
          output2 <- run.jags(data = win.data,  inits = inits, monitor = params, burnin = nb, model = "model.txt", sample = ni, n.chains = nc, method = "parallel")
          
          # Save the output
          
          store[q, 1, Nsheet] <- sodata$psi
          
          store[q, 2, Nsheet] <- sodata$p
          
          store[q, 3, Nsheet] <- sodata$T
          
          store[q, 4, Nsheet] <- sodata$R
          
          store[q, c(5,6,7,8), Nsheet] <- c(output2$summaries[1,4], output2$summaries[1,5], output2$summaries[1,1], output2$summaries[1,3])
          
        }				
        Nsheet <- Nsheet + 1	
      }
    }
  }
}

