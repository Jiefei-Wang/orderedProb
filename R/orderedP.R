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


orderedProb <- function(l,h){
    total <- sort(c(0,l,h,1))
    #g(t_i)
    g_value <- get_g(total,h)
    #h(t_i)
    h_value <- get_h(total,l)
    
    n_t <- length(total)
    diff_t <- diff(total)
    m <- length(l)
    # record <-c()
    Q <- rep(0, m+1)
    Q[1] <- 1
    for(i in seq_len(n_t-1)-1){
        # record <- rbind(record,Q)
        gt_i <- g_value[i+1]
        ht_i_plus<- h_value[i+2]
        
        index <- (gt_i+1) : (ht_i_plus-1) +1
        tmpQ <- Q[index]
        tmpPi<-dpois((ht_i_plus- gt_i - 2):0,m*diff_t[i+1])
        res <- convolve(tmpQ, tmpPi, type = "o")
        
        # Q[] <- 0
        Q[index] <- res[seq_len(length(index))]
    }
    Q[m+1]/dpois(m,m)
}

