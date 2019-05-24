# Additional functions for confiance intervals and tests, from M.A. Remiche's lecture.
# By: Pierre Beaujean (April 2019).

# --- CI = confiance interval

qz <- function(q, df=NULL) {
    return(qnorm(q, 0, 1))
}

pz <- function(p, df=NULL, lower.tail=TRUE) {
    return(pnorm(p, 0, 1, lower.tail=lower.tail))
}

CI.mean <- function(x, sd=NULL, conf.level=.95)
{
    mu <- mean(x)
    n <- length(x)
    empirical_sd <- FALSE
    
    if (is.null(sd)) {
        sd <- sd(x)
        empirical_sd <- TRUE
    }

    return(CI.mean.base(n, mu, sd, empirical_sd, conf.level))
}

CI.mean.base <- function(n, mu, sd, empirical_sd, conf.level=.95) {
    if (empirical_sd) {
        e = qt((1+conf.level)/2, n-1) * sd / sqrt(n)
    } else {
        e = qz((1+conf.level)/2) * sd / sqrt(n)
    }
    return(c(mu-e, mu+e))
}

CI.var.base <- function(n, v, conf.level=.95) {
    return(c((n-1)*v/qchisq((1+conf.level)/2, n-1), (n-1)*v/qchisq((1-conf.level)/2, n-1)))
}

CI.var <- function(x, conf.level=.95)
{
    v <- var(x)
    n <- length(x)

    return(CI.var.base(n, v, conf.level))
}

CI.prop.base <- function(n, prop.char, conf.level=.95) {
    e = qnorm((1+conf.level)/2, 0, 1) * sqrt((prop.char * (1 - prop.char)) / n)
    return(c(prop.char-e, prop.char+e))
}

# --- hypothesis tests

test.base <- function(t, df, pfunc, qfunc, vars, risk.level=0.05, alternative=c("two.sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    
    e1 <- -Inf
    e2 <- Inf
    
    p.value <- Inf
    
    if (alternative == "two.sided") {
        e1 = qfunc(risk.level/2, df)
        e2 = qfunc(1-risk.level/2, df)
        
        p.value = 2*min(c(pfunc(abs(t), df), pfunc(abs(t), df, lower.tail=FALSE)))
        # p.value = pfunc(-abs(t), df) + pfunc(abs(t), df, lower.tail=FALSE)
        
        cat(sprintf("| H0: %s == %s\n", vars[1], vars[2]))
        cat(sprintf("| H1: %s != %s\n", vars[1], vars[2]))

    } else {
        if (alternative == "less") {
            e1 = qfunc(risk.level, df)
        
            p.value = pfunc(t, df)
            
            cat(sprintf("| H0: %s >= %s\n", vars[1], vars[2]))
            cat(sprintf("| H1: %s <  %s\n", vars[1], vars[2]))
        } else {
            e2 = qfunc(1-risk.level, df)
        
            p.value = pfunc(t, df, lower.tail=FALSE)
            
            cat(sprintf("| H0: %s <= %s\n", vars[1], vars[2]))
            cat(sprintf("| H1: %s >  %s\n", vars[1], vars[2]))
        }
    }
    
    if(length(df) != 1) {
        cat(sprintf("df = (%s)", paste(df,collapse=",")))
    } else {
        cat(sprintf("df = %d", df))
    }
    cat(", ")
    if (log(p.value) > -4)  {
        cat(sprintf("p-value = %.5f", p.value))
    } else {
        cat(sprintf("p-value = %.3e", p.value))
    }
    cat(sprintf("\nIC = [%.5f; %.5f]", e1, e2))
    
    if (t < e2 && t > e1) {
        cat(" -> H0 cannot be rejected\n")
        return(FALSE)
    } else {
        cat(" -> H0 can be rejected, alternative (H1) is true\n")
        return(TRUE)
    }
    
    return(NULL)
}

test.mean.normal <- function(x, mu_0, sd=NULL, risk.level=0.05, alternative="two.sided") {
    mu_n <- mean(x)
    n <- length(x)
    cat(sprintf("mu_n = %.5f\n", mu_n))
    empirical_sd <- FALSE
    
    if (is.null(sd)) {
        sd <- sd(x)
        empirical_sd <- TRUE
    }
    
    return(test.mean.normal.base(mu_0, n, mu_n, sd, empirical_sd, risk.level, alternative))
}

test.mean.normal.base <- function(mu_0, n, mu_n, sd, empirical_sd, risk.level=0.05, alternative="two.sided") {
    t <- (mu_n - mu_0) * sqrt(n) / sd
    cat(sprintf("T(x) = %.5f\n", t))
    
    if (empirical_sd) {
        cat("NOTE: use UNKNOWN sd version!\n")
        return(test.base(t, n-1, pt, qt, c("mu_n", "mu_0"), risk.level, alternative))
    } else {
        cat("NOTE: use KNOWN sd version!\n")
        return(test.base(t, n, pz, qz, c("mu_n", "mu_0"), risk.level, alternative))
    }
}

test.mean <- function(x, mu_0, sd=NULL, risk.level=0.05, alternative="two.sided") {
    alternative <- match.arg(alternative)
    mu_n <- mean(x)
    n <- length(x)
    cat(sprintf("mu_n = %.5f\n", mu_n))
    empirical_sd <- FALSE
    
    if (is.null(sd)) {
        sd <- sd(x)
        empirical_sd <- TRUE
    }
    
    return(test.mean.base(mu_0, n, mu_n, sd, empirical_sd, risk.level, alternative))
}

test.mean.base <- function(mu_0, n, mu_n, sd, empirical_sd, risk.level=0.05, alternative="two.sided") {
    alternative <- match.arg(alternative)
    if (empirical_sd) {
        t <- (mu_n - mu_0) * sqrt(n - 1) / sd
        cat(sprintf("T(x) = %.5f\n", t))
        cat("NOTE: use UNKNOWN sd version!\n")
        return(test.base(t, n, pz, qz, c("mu_n", "mu_0"), risk.level, alternative))
    } else {
        t <- (mu_n - mu_0) * sqrt(n) / sd
        cat(sprintf("T(x) = %.5f\n", t))
        cat("NOTE: use KNOWN sd version!\n")
        return(test.base(t, n, pz, qz, c("mu_n", "mu_0"), risk.level, alternative))
    }
}

test.var <- function(x, var_0, mu=NULL, risk.level=0.05, alternative="two.sided") {
    var_n <- var(x)
    n <- length(x)
    cat(sprintf("var_n = %.5f\n", var_n))
    
    if (is.null(mu)) {
        mu <- mean(x)
        return(test.var.base(var_0, n, var_n, mu, risk.level, alternative))
    } else {  # known mu
        t <- sum((x - mu)^2 / var_0)
        cat(sprintf("T(x) = %.5f\n", t))
        cat("NOTE: use KNOWN mean version!\n")
        return(test.base(t, n, pchisq, qchisq, c("var_n", "var_0"), risk.level, alternative))
    }
}

test.var.base <- function(var_0, n, var_n, mu, risk.level=0.05, alternative="two.sided") {
    
    t <- (n-1)*var_n / var_0
    cat(sprintf("T(x) = %.5f\n", t))
    cat("NOTE: use UNKNOWN mean version!\n")
    
    return(test.base(t, n-1, pchisq, qchisq, c("var_n", "var_0"), risk.level, alternative))
}

test.compare_means <- function(x, y, risk.level=0.05, alternative="two.sided") {
    mean_x <- mean(x)
    var_x <- var(x)
    n_x <- length(x)
    mean_y <- mean(y)
    var_y <- var(y)
    n_y <- length(y)
    
    return(test.compare_means.base(n_x, mean_x, var_x, n_y, mean_y, var_y, risk.level, alternative))
}


test.compare_means.base<- function(n_x, mean_x, var_x, n_y, mean_y, var_y, risk.level=0.05, alternative="two.sided") {
    cat(sprintf("x :: n = %d, mean = %.5f, var = %.5f\n", n_x, mean_x, var_x))
    cat(sprintf("y :: n = %d, mean = %.5f, var = %.5f\n", n_y, mean_y, var_y))
    
    df <- (var_x/(n_x-1)+var_y/(n_y-1))^2 / (var_x^2/((n_x-1)*n_x^2)+var_y^2/((n_y-1)*n_y^2))
    
    dfx <- df - as.integer(df)
    if (dfx < 0.5) {
        ndf <- as.integer(df)
    } else {
        ndf <- as.integer(df + 1)
    }
    
    cat(sprintf("nu = %.5f\n", df)) 
    t <- (mean_x - mean_y) * (var_x / (n_x-1) + var_y / (n_y-1))^(-1/2)
    cat(sprintf("T(x,y) = %.5f\n", t))
    
    return(test.base(t, ndf, pt, qt, c("mean_x", "mean_y"), risk.level, alternative))
}

test.compare_vars <- function(x, y, risk.level=0.05, alternative="two.sided") {
    var_x <- var(x)
    n_x <- length(x)
    var_y <- var(y)
    n_y <- length(y)
    
    return(test.compare_vars.base(n_x, var_x, n_y, var_y, risk.level, alternative))
}

qfmod <- function(q, df) {
    return(qf(q,  df1=df[1], df2=df[2]))
}

pfmod <- function(p, df, lower.tail=TRUE) {
    return(pf(p, df1=df[1], df2=df[2], lower.tail=lower.tail))
}


test.compare_vars.base<- function(n_x, var_x, n_y, var_y, risk.level=0.05, alternative="two.sided") {
    cat(sprintf("x :: n = %d, var = %.5f\n", n_x, var_x))
    cat(sprintf("y :: n = %d, var = %.5f\n", n_y, var_y))
    
    t <- var_x / var_y
    cat(sprintf("T(x,y) = %.5f\n", t))
    
    return(test.base(t, c(n_x-1, n_y-1), pfmod, qfmod, c("var_x", "var_y"), risk.level, alternative))
}

test.equal_distributions <- function(x, y, risk.level=0.05) {
    n = sum(x)
    
    if (length(x) != length(y)) {
        message("x and y must contain the same number of elements")
        return(NULL)
    }
    
    py = y / sum(y)
    npi = n*py
    
    # check limit
    if(n < 30) {
        message("WARNING: statistic may be inapropriate: n<30")
    }
    
    if (length(npi[npi < 1]) > 0) {
        message("WARNING: statistic may be inapropriate, at least one n*p(i) < 1 -> group some categories")
    }
    
    if (length(npi[npi < 5]) > .2 * length(x)) {
        message("WARNING: statistic may be inapropriate, more than 20% of n*p(i) < 5 -> group some categories")
    }
    
    # go ahead !
    t <- sum((x - npi)^2 / npi)
    cat(sprintf("T(x,y) = %.5f\n", t))
    df = length(x)-1
    
    e = qchisq(1-risk.level, df)
    p.value = pchisq(t, df, lower.tail=FALSE)
    cat(sprintf("df = %d, p-value = %.5f\n", df, p.value))
    
    cat("| H0: x follows y distribution\n| H1: x doesn't follow y distribution\n")
    
    cat(sprintf("IC = [-Inf;%.5f]", e))
    
    if (t < e) {
        cat(" -> H0 cannot be rejected\n")
        return(FALSE)
    } else {
        cat(" -> H0 can be rejected, alternative (H1) is true\n")
        return(TRUE)
    }
    
    return(NULL)
}

test.independance <- function(x, risk.level=0.05){
    if(!is.table(x)){
        message("ERROR: x is not a table")
        return(NULL)
    }
    
    if(length(dim(x)) != 2) {
        message("ERROR: x is not a 2 dimension table")
        return(NULL)
    }
    
    nr = nrow(x)
    nc = ncol(x)
    n = sum(x)
    
    cat("data\n----\n")
    print(addmargins(x))
    
    ft = prop.table(x)
    extended_ft = addmargins(ft)
    
    p_dotj = extended_ft[nr + 1, 1:nc]
    p_idot = extended_ft[1:nr, nc + 1]
    
    cat("\n")
    cat("p_i*\n----\n")
    print(p_idot)
    cat("\n")
    cat("p_*j\n----\n")
    print(p_dotj)
    
    cat('\nResults\n-------\n')
    
    t <- 0
    
    for(i in c(1:nr)) {
        for(j in c(1:nc)) {
            e <- n*p_idot[i]*p_dotj[j]
            t <- t + (x[i,j]-e)^2 / e
        }
    }
    
    cat(sprintf("T(x,y) = %.5f\n", t))
    df = (nr-1)*(nc-1)
    
    e = qchisq(1-risk.level, df)
    p.value = pchisq(t, df, lower.tail=FALSE)
    
    cat(sprintf("df = %d, p-value = %.5f\n", df, p.value))
    
    cat("| H0: x and y are independants\n| H1: x and y are dependant\n")
    
    cat(sprintf("IC = [-Inf;%.5f]", e))
    
    if (t < e) {
        cat(" -> H0 cannot be rejected\n")
        return(FALSE)
    } else {
        cat(" -> H0 can be rejected, alternative (H1) is true\n")
        return(TRUE)
    }
    
    return(NULL)
    
    
}


