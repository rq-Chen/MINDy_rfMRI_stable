# Modified from source code to fix a bug
my_wrskGap <- function (data, K, S, npermute = 10, cores)
{
    rownames(data) <- 1:nrow(data)
    colnames(data) <- 1:ncol(data)
    gap_s <- se_s <- logBsa <- logBs <- c()
    B_s <- nonZerW <- c()
    resFinal <- list()
    for (j in 1:length(S)) {
        cat(S[j])
        ress <- wrsk(dat = data, k = K, s = S[j])  # The whole point is to suppress warnings here
        resFinal[[j]] <- ress
        B_s <- ress$WBCSS
        nonZerW <- c(nonZerW, length(which(ress$varweights != 
            0)))
        B_sa <- numeric(npermute)
        res.per <- mclapply(1:npermute, function(x, dat, k, s, 
            weights_par) {
            set.seed(x)
            dat.per <- apply(dat, 2, function(z) {
                sample(z, length(z))
            })
            r <- wrsk(dat.per, k, s)
        }, dat = data, k = K, s = S[j], mc.cores = cores)
        B_sa <- numeric(npermute)
        for (a in 1:npermute) {
            B_sa[a] <- res.per[[a]]$WBCSS
        }
        gap_s <- c(gap_s, log(B_s) - mean(log(B_sa)))
        se_s <- c(se_s, 1 * (sqrt(1 + 1/length(B_sa)) * sqrt(var(log(B_sa))) * 
            sqrt((length(B_sa) - 1)/length(B_sa))))
        if (nonZerW[j] == ncol(data) | S[j] == tail(S, 1)) {
            break
            return(list(gap = gap_s, se = se_s, nonZerW = nonZerW, 
                s = S[1:j], resFinal = resFinal))
        }
    }
    return(list(gap = gap_s, se = se_s, nonZerW = nonZerW, s = S[1:j], 
        resFinal = resFinal))
}
