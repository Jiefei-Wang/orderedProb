library(KSgeneral)
l=c(0.1,0.79)
h=c(0.5,0.8)

stat <- 0.1
n=40000

l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
h=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))

total <- sort(c(0,l,h,1))


#g(t_i)
g_value <- get_g(total,h)
#h(t_i)
h_value <- get_h(total,l)
n_t <- length(total)
diff_t <- diff(total)
m <- length(l)
compute_prob_fft(m,g_value,h_value,n_t,diff_t)


compute_prob(m,g_value,h_value,n_t,diff_t)



Q <- matrix(rep(0, m+1),1)
Q[1,1] <- 1
for(i in seq_len(n_t-1)-1){
    gt_i <- g_value[i+1]
    gt_i_plus <- g_value[i+2]
    ht_i_plus<- h_value[i+2]
    newQ <- rep(0, m+1)
    for(curM in (gt_i_plus+1) : (ht_i_plus-1)){
        l_list <- seq.int(gt_i+1,curM)
        Q_list <- Q[i+1,l_list+1]
        pi_list <- dpois(curM-l_list,m*diff_t[i+1])
        newQ[curM+1] <- sum(Q_list*pi_list)
    }
    Q <- rbind(Q,newQ)
}
Q[n_t,m+1]/dpois(m,m)


library(fftw)
# record <-c()
Q <- rep(0, m+1)
Q[1] <- 1
for(i in seq_len(n_t-1)-1){
    # record <- rbind(record,Q)
    gt_i <- g_value[i+1]
    ht_i_plus<- h_value[i+2]
    
    index <- (gt_i+1) : (ht_i_plus-1) +1
    tmpQ <- Q[index]
    tmpPi<-dpois(0:(ht_i_plus- gt_i - 2),m*diff_t[i+1])
    Q_f= FFT(c(tmpQ,rep(0,length(index))))
    Pi_f = FFT(c(tmpPi,rep(0,length(index))))
    
    Q[index] <- Re(IFFT(Q_f*Pi_f)[seq_len(length(index))])
    
}
Q[m+1]/dpois(m,m)





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
            Q_list <- Q[l_list+1]
            pi_list <- dpois(curM-l_list,m*diff_t[i+1])
            newQ[curM+1] <- sum(Q_list*pi_list)
        }
        Q <- newQ
    }
    Q[m+1]/dpois(m,m)
})
# Q


system.time({
    a= compute_prob(m,g_value,h_value,n_t,diff_t)
})
system.time({
    b= compute_prob_fft(m,g_value,h_value,n_t,diff_t)
})


stat <- 0.1
n=10000

l=sapply(1:n,function(x)qbeta(stat,x,n-x+1))
h=sapply(1:n,function(x)qbeta(1 - stat,x,n-x+1))

orderedProb(l,h)


df <- data.frame(rbind(h, l))
write.table(df, "Boundary_Crossing_Time.txt", 
            sep = ", ", row.names = FALSE, col.names = FALSE)
1 - KSgeneral::ks_c_cdf_Rcpp(ncol(df))
file.remove("Boundary_Crossing_Time.txt")



system.time({a=orderedProb(l,h)})
system.time({b=1 - KSgeneral::ks_c_cdf_Rcpp(ncol(df))})







devtools::load_all()

x_real <- as.double(x)
x_img<-as.double(rep(0,length(x_real)))
y_real <- as.double(y)
y_img <- as.double(rep(0,length(x_real)))

round(convolve(x_real, rev(y_real), type = "o"),3)
test(x_real,x_img,rev(y_real),y_img)

x_real

testfft(x)
x<-c(1,2,3,4)
fft(x)

main1()

x<-c(1,2,3,4,0,0,0,0)
y<-c(4,5,6,7,0,0,0,0)

testConv(x,y)
x

x_f= fft(x)
y_f = fft(y)

round(Re(fft(x_f*y_f,inverse=T)),3)/length(x)

x<-c(1,2,3,4,0,0)
testfft(x)


fft(x)
