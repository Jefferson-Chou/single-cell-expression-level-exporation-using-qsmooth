#scale normalization
object <- scale(object,center = T,scale = T)

#compute sample quantile 
Q = apply(object, 2, sort) 

#quantile normalization begin by calculate a reference distribution
# compute reference quantile (Qref)
Qref = apply(Q, 1, mean)


#compute the estimated regression coefficients(Qbeta)
if(is.factor(group_factor)){
  X = model.matrix(~ 0 + group_factor)
} else {
  X = model.matrix(~ group_factor)
}
Qbeta = t(solve(t(X) %*% X) %*% t(X) %*% t(Q))

#Qhat
Qhat = Qbeta %*% t(X)

#compute SST¡¢SSB
SST = rowSums((Q - Qref)^2)
SSB = rowSums((Qhat - Qref)^2)

#compute weight
Weights = 1 - SSB / SST

#compute k (to smooth the weights) (window=0.05)
window = 0.05
k = floor(nrow(Q) * window)
if (k %% 2 == 0) 
  k = k + 1
w = runmed(Weights, k = k, endrule="constant")

#compute the smooth quantile normalized data
Qsmooth = w * Qref + (1 - w) * Qhat


