# Rarefaction Analysis


# Load GT.TABLE from combined 012 files
data.files <- paste0(
    '../data/rarefaction/Spatial_144_ind_rarefaction_012_format_',
    c('LRV.012.GT.TABLE', 'OHN.012.GT.TABLE', 'OM3.012.GT.TABLE', 'TER1.012.GT.TABLE', 
      'ZW2.012.GT.TABLE', 'MO.012.GT.TABLE',  'OM2.012.GT.TABLE', 'TER2.012.GT.TABLE', 'ZW4.012.GT.TABLE')
)

data.sets <- lapply(data.files, read.delim)


# solving NA issues
clean.na <- function(mtx){
    cat('# lines =', dim(mtx)[1], '\n')
    cat('# lines without NA =', dim(na.omit((mtx)))[1], '\n')
    cat('% lines without NA =', nrow(na.omit((mtx)))/nrow(mtx), '\n')
    cat('# NAs =', length(mtx[is.na(mtx)]), '\n\n')
    mtx[is.na(mtx)] <- '0'
    mtx[,3:ncol(mtx)] <- apply(mtx[,3:ncol(mtx)], 2, as.numeric)
    stopifnot(length(mtx[is.na(mtx)])==0)
    return(mtx)
}

data.sets <- lapply(data.sets, clean.na)


# selecting the input file ####
final.selection.1 <- function(mtx){
    input0 <- na.omit(mtx)#[sample(nrow(na.omit(mtx)),2000),] 
    input <- na.omit(input0[,3:ncol(mtx)])
    ngtype <- ncol(input)
    input$bigA <- rowSums(input)
    input$polymorphic <- ifelse(input$bigA %in% c(1:2*ngtype-1), 1, 0)
    return(input)
}

final.selection.2 <- function(mtx){
    input0 <- na.omit(mtx)#[sample(nrow(na.omit(mtx)),2000),] 
    input <- input0[,3:ncol(mtx)]
    ngtype <- ncol(input)
    LOCI <- function(x) {sum(!is.na(x))}
    GENOTYPES <- function(x) {length(levels(factor(unlist(x))))}
    input$loci <- apply(input[,1:ngtype], 1, LOCI)
    input$genotypes <- apply(input[,1:ngtype], 1, GENOTYPES)
    input$bigA <- rowSums(input[,1:ngtype], na.rm=TRUE)
    input$bigAfreq <- input$bigA / (2*input$loci)
    input$polymorphic <- ifelse(!input$bigAfreq %in% c(0,1), 1, 0)
    input$polymorphic2 <- ifelse(input$genotypes > 1, 1, 0)
    return(input)
}

inputs <- lapply(data.sets, final.selection.2)


# calculate the PPLs
PPL.1 <- function(x,mtx){
    X <- subset(mtx, select=x)
    genes <- ncol(X)*2
    X$bigA <- rowSums(X)
    table(X$bigA) # values between 0 and genes
    X$polymorphic <- ifelse(X$bigA %in% c(1:(genes-1)), 1, 0)
    ppl <- sum(X$polymorphic) / nrow(X) # percentage polymorphic loci
    return(ppl)
}

PPL.2 <- function(x,mtx){
    X <- subset(mtx, select=x)
    genes <- ncol(X)*2
    LOCI <- function(x) {rowSums(x, na.rm=TRUE)}
    Y <- data.frame(X[,1:(genes/2)])
    Z <- Y; Z[!is.na(Z)] <- 1
    X$loci <- LOCI(Z)
    X$bigA <- rowSums(Y, na.rm=TRUE)
    X$bigAfreq <- X$bigA / (2*X$loci)
    X$polymorphic <- ifelse(is.na(X$bigAfreq), NA, ifelse(!X$bigAfreq %in% c(0,1), 1, 0))
    ppl <- sum(X$polymorphic,na.rm=TRUE) / sum(!is.na(X$polymorphic))
    return(ppl)
}

apply.PPL <- function(input){
    ngtype <- ncol(input)-6

    if (ngtype > 15) {
        # Random down sampling
        LIST0 <- lapply(1:3,      function(x) combn(ngtype,x)) # MUST ADJUST BASED ON THE VALUE OF ngtype
        LIST1 <- lapply(4:ngtype, function(x) sapply(1:1000, function(i) sample(ngtype,x)))
        LIST <- c(LIST0, LIST1)
    } else {
        # OR full prmutation
        LIST <- lapply(1:ngtype, function(x) combn(ngtype,x))
    }
                       
    ppl <- function(lst) lapply(as.list(data.frame(lst)), function(x) PPL.2(x,input))
    FINAL0 <- data.frame(do.call(rbind, lapply(1:length(LIST), function(i) cbind(ppl(LIST[[i]]),i))))
    FINAL <- data.frame(unlist(FINAL0[,1]), unlist(FINAL0[,2]))
    colnames(FINAL) <- c('ppl','samplesize')
    return(FINAL)
}

finals <- lapply(inputs, apply.PPL)

                         
# save results
save(finals, file = '../data/rarefaction/rarefaction_Spatial_144_ind_rarefaction_012_format.RData')

                         
# plot comparisons                         
cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm,as.matrix)
    n <- max(sapply(nm,nrow)) 
    do.call(cbind, lapply(nm, function(x) rbind(x, matrix(,n-nrow(x),ncol(x))))) 
}

plot.all <- function(finals, maxngt, labs){
    finals <- lapply(finals, function (x) x[x$samplesize <= maxngt,])
    nplots <- length(finals)
    
    NEW <-data.frame(cbind(
        1:maxngt, 
        do.call(cbind.fill,
                lapply(finals, function (x) do.call(rbind,as.list(by(x$ppl,x$samplesize,mean)))))
    ))

    plot(NEW$X1, seq(0, 1, length.out=length(NEW$X1)),
         xlab='Number of genotypes', ylab='Proportion polymorphic loci', 
         ylim=c(0,1), xaxt='n', type='n', pch=16) # ,main="random loci in T"
    axis(side=1, at=1:maxngt, labels=1:maxngt, cex.axis=0.8)

    sf <- seq(0, 0.2, length.out=nplots) / nplots
    cl <- rainbow(nplots)
    for (i in 1:nplots){
        x <- NEW$X1
        y <- NEW[,i+1]
        x <- x[!is.na(y)]
        y <- y[!is.na(y)]
        points(finals[[i]]$samplesize+i*sf[i]-0.1, finals[[i]]$ppl, col=cl[i], pch=16, cex=0.4, lwd=2)
        lines(x, y, col=cl[i], pch=16, cex=0.4, lty='dashed')
    }

    legend('bottomright', legend=labs, col=cl[1:nplots], bty='n', pch=16, cex=0.8, lty='dashed')
}
              
pdf(file='../data/rarefaction/rarefaction_Spatial_144_ind_rarefaction_012_format.pdf', width=9, height=7)
plot.all(finals, 12, c('LRV', 'OHN', 'OM3', 'TER1', 'ZW2', 'MO', 'OM2', 'TER2', 'ZW4'))
dev.off()
