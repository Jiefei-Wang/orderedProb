library(KSgeneral)
l=c(0.1,0.79)
h=c(0.5,0.8)

stat <- 0.1
n=1000

l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
h=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))

total <- sort(c(0,l,h,1))


get_g <- function(value,func_data){
    n<- length(value)
    j<- 0
    result <- rep(0,n)
    for(i in seq_len(n)){
        if(value[i]==func_data[j+1]){
            j = j+1
            result[i] = j - 1
        }else{
            result[i] = j - 1
        }
        if(j == length(func_data)&&i!=n){
            result[(i+1):n] = length(func_data)-1
            break
        }
    }
    result
}
get_h <- function(value,func_data){
    n<- length(value)
    j<- 0
    result <- rep(0,n)
    for(i in seq_len(n)){
        if(value[i]==func_data[j+1]){
            j = j+1
            result[i] = j
        }else{
            result[i] = j + 1 
        }
        if(j == length(func_data)&&i!=n){
            result[(i+1):n] = length(func_data)+1
            break
        }
    }
    result
}
#g(t_i)
g_value <- get_g(total,h)
#h(t_i)
h_value <- get_h(total,l)

system.time({
    n_t <- length(total)
    diff_t <- diff(total)
    m <- length(l)
    # record <-c()
    Q <- rep(0, m+1)
    Q[1] <- 1
    for(i in seq_len(n_t-1)-1){
        # record <- rbind(record,Q)
        gt_i <- g_value[i+1]
        gt_i_plus<- g_value[i+2]
        ht_i_plus<- h_value[i+2]
        
        index <- (gt_i+1) : (ht_i_plus-1) +1
        tmpQ <- Q[index]
        tmpPi<-dpois((ht_i_plus- gt_i - 2):0,m*diff_t[i+1])
        res <- convolve(tmpQ, tmpPi, type = "o")
        
        # Q[] <- 0
        Q[index] <- res[seq_len(length(index))]
    }
    Q[m+1]/dpois(m,m)
})
i=0
q<- c(Q[i+1,],rep(0,m+1))
mypi<- c(dpois(0:m,m*diff_t[i+1]),rep(0,m+1))

fq<- fft(q)
fpi<- fft(mypi)

Re(fft(fq*fpi,TRUE))/length(fq)
Q[i+2,]

system.time({
    Q <- rep(0, m+1)
    Q[1] <- 1
    # Q <- matrix(rep(0, m+1),1)
    # Q[1,1] <- 1
    for(i in seq_len(n_t-1)-1){
        gt_i <- g_value[i+1]
        gt_i_plus <- g_value[i+2]
        ht_i_plus<- h_value[i+2]
        newQ <- rep(0, m+1)
        for(curM in (gt_i_plus+1) : (ht_i_plus-1)){
            l_list <- seq.int(gt_i+1,curM)
            # Q_list <- Q[i+1,l_list+1]
            Q_list <- Q[l_list+1]
            pi_list <- dpois(curM-l_list,m*diff_t[i+1])
            newQ[curM+1] <- sum(Q_list*pi_list)
        }
        # Q <- rbind(Q,newQ)
        # Q <- newQ
    }
    Q[m+1]/dpois(m,m)
})
# Q

# Q[n_t,m+1]/dpois(m,m)


df <- data.frame(rbind(h, l))
write.table(df, "Boundary_Crossing_Time.txt", 
            sep = ", ", row.names = FALSE, col.names = FALSE)
1 - KSgeneral::ks_c_cdf_Rcpp(ncol(df))
file.remove("Boundary_Crossing_Time.txt")


system.time({1 - KSgeneral::ks_c_cdf_Rcpp(ncol(df))})
