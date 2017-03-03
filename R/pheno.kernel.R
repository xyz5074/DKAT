pheno.kernel <-
function(Y, rho=0.1){
s=var(Y)
N=glasso(s,rho=0.1)$wi
Ky=Y%*%N%*%t(Y)
return(Ky)  
}
