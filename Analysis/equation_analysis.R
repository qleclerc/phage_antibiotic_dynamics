
phage_q = function(P, N, beta=1) (1-exp(-beta*N))*P

phage_alt = function(P, N, beta) beta*P*N
plot(seq(0,100,1),phage_alt(10,seq(0,100,1),1),
     type= "p", xlab = "Bacteria", ylab="Phage binding / bacteria infected", main = "10 phage in system")
plot(seq(0,100,1),phage_alt(seq(0,100,1),10,1),
     type= "p", xlab = "Phage", ylab="Phage binding / bacteria infected", main = "10 bacteria in system")


#the one below doesn't work, because implicitly assumes that 1 phage at most can bind to a bacteria
#so, number of phage binding saturates quickly
#phage_alt = function(P, N, beta) (1-exp(-beta*N/P))*P

bac_q = function(P, N, beta2=1) (1-exp(-beta2*P/N))*N

#the one below doesn't work because implicitly assumes phage can bind to more than one bacteria over 1 time point
#so number of bac infected can go over number of bac in system
#bac_alt = function(P, N) (1-exp(-N/P))*N

plot(seq(0,100,1), phage_q(10,seq(0,100,1),1),
     type= "p", xlab = "Bacteria", ylab="Phage binding", main = "10 phage in system")
plot(seq(0,100,1), phage_q(seq(0,100,1),10,1),
     type= "p", xlab = "Phage", ylab="Phage binding", main = "10 bacteria in system")

plot(seq(0,100,1), bac_q(10,seq(0,100,1)),
     type= "p", xlab = "Bacteria", ylab="Bacteria infected", main = "10 phage binding in system")
plot(seq(0,100,1), bac_q(seq(0,100,1),10),
     type= "p", xlab = "Phage binding", ylab="Bacteria infected", main = "10 bacteria in system")

plot(phage_q(seq(0,10,1),10,1), type= "p", xlab = "Bacteria", ylab="Phage binding", main = "10 bacteria in system")


scenario = function(P, N, beta=1, beta2=1){
  P_inf = phage_q(P,N,beta)
  B_inf = bac_q(P_inf, N,beta2)
  
  cat(N, "bacteria and", P, "phage, binding at a rate of", beta, ", and infecting at a rate of", beta2, ":\n",
      P_inf, "phage bind and", B_inf, "bacteria are infected")
  
}

scenario(10,100,0.5,2)

#














chosen_beta = 1
tt = matrix(0, nrow=20, ncol=20)
rownames(tt) = paste0("P", as.character(seq(1,nrow(tt),1)))
colnames(tt) = paste0("B", as.character(seq(1,ncol(tt),1)))

for(i in 1:nrow(tt)){
tt[i,] = phage_q(i, seq(1,ncol(tt),1), chosen_beta)/i
}

pheatmap(tt, cluster_cols = F, cluster_rows = F)


tt = matrix(0, nrow=30, ncol=30)
rownames(tt) = paste0("P", as.character(seq(1,nrow(tt),1)))
colnames(tt) = paste0("B", as.character(seq(1,ncol(tt),1)))

for(i in 1:nrow(tt)){
  tt[i,] = bac_q(phage_q(i, seq(1,ncol(tt),1), chosen_beta),
                 seq(1,ncol(tt),1))/seq(1,ncol(tt),1)
}

pheatmap(tt, cluster_cols = F, cluster_rows = F)


tt = matrix(0, nrow=30, ncol=30)
rownames(tt) = paste0("P", as.character(seq(1,nrow(tt),1)))
colnames(tt) = paste0("B", as.character(seq(1,ncol(tt),1)))

for(i in 1:nrow(tt)){
  tt[i,] = i*seq(1,ncol(tt),1)*chosen_beta/seq(1,ncol(tt),1)
}

pheatmap(tt, cluster_cols = F, cluster_rows = F)
