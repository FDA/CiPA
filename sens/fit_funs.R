try_polr<-function(datadf, predname, penalty=0, tol=NA, maxit=1e6){
    if(length(penalty)>1 || penalty!=0)
        print("Ignoring penalty parameter for polr function!")
    if(!is.na(tol))
        print("Ignoring tolerance parameter for polr function!")

    # get starting guess (polr's default often fails)
    u <- as.integer(table(datadf$class))
    u <- (cumsum(u)/sum(u))[1:(length(u)-1)]
    zetas <- qlogis(u)
    s0 <- c(0, zetas[1], log(diff(zetas)))

    try({
        polr(as.formula(paste0("class~",predname)), data=datadf, start=s0, Hess=TRUE, control=list(maxit=maxit)) # note: large maxit is very slow but seems necessary for convergence
    })
}

try_glm<-function(datadf, predname, penalty=0, tol=NA, maxit=1e6){
    if(length(penalty)>1 || penalty!=0)
        print("Ignoring penalty parameter for polr function!")
    if(!is.na(tol))
        print("Ignoring tolerance parameter for polr function!")

    try({
        glm(as.formula(paste0("class~",predname)), data=datadf, family=binomial(link="logit"), control=list(maxit=maxit))
    })
}

try_repolr<-function(datadf,predname,penalty=0,tol=1e-10,maxit=1e6){
    
   try({
        repolr(as.formula(paste0("class~",predname)), subjects="drug",data=datadf, categories=3,times=1:maxsamp,corr.mod="independence",alpha=0.9,diffmeth="analytic", fit.opt=c(cmaxit = maxit, omaxit = maxit, ctol = tol, otol = tol, h = 0.01))
    })

}#repolr

try_ordLORgee<-function(datadf, predname, penalty=0,tol=1e-10,maxit=1e6){
     options(warn=1)
   try({
        ordLORgee(as.formula(paste0("class~",predname)), id=datadf$drug,data=datadf, link="logit", repeated=datadf$sample,LORstr="independence",control=LORgee.control(verbose=FALSE))
    })

}#ordLORgee


try_lrm<-function(datadf, predname, penalty, tol=1e-10, maxit=1e6){
    if(class(penalty)=="list" || length(penalty)==1){ # specified penalty
        fout<-try({
            lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=penalty, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
        })
        return(fout)
    }

    pen_vec<-penalty # candidate penalties to try if penalty=0 doesn't work
    fout<-try({
        lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=0, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
    })
    if(inherits(fout, "try-error")){ # try all the penalties
        print("Unpenalized regression didn't work, using first penalty that works!")
        i<-1
        while(inherits(fout, "try-error") && i<=length(pen_vec)){
            pen<-pen_vec[i]
            i<-i+1
            fout<-try({
                lrm(as.formula(paste0("class~",predname)), data=datadf, penalty=pen, x=TRUE, y=TRUE, tol=tol, maxit=maxit)
            })
        }
        #print(sprintf("Using penalty = %g",pen))
    }
    if(inherits(fout, "try-error"))
        return(fout)

    p <- pentrace(fout, penalty=pen_vec, tol=tol, maxit=maxit, noaddzero=TRUE)
    if(inherits(p, "try-error"))
        return(p)
    #print(sprintf("Using penalty = %g",p$penalty))
    update(fout, penalty=p$penalty)
}
