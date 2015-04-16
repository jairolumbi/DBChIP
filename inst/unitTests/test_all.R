test_all <- function(){
    data("PHA4")
    conds <- factor(c("emb","emb","L1", "L1"), levels=c("emb", "L1"))
    
    bs.list <- read.binding.site.list(binding.site.list)
    
    checkTrue(is.list(bs.list))
    #each chr in each dataset is a data.frame
    checkTrue(is.data.frame(bs.list[[1]][[1]]))
    #the example has two fields
    checkTrue(identical(colnames(bs.list[[1]][[1]]), c("pos", "weight")))

    ## compute consensus site
    consensus.site <- site.merge(bs.list, in.distance=100, out.distance=250)
    
    checkTrue(is.data.frame(consensus.site[[1]]))
    checkTrue(identical(colnames(consensus.site[[1]]), c("pos", "nsig", "origin", "ori.pos")))
    
    dat <- load.data(chip.data.list=chip.data.list, conds=conds, consensus.site=consensus.site, input.data.list=input.data.list, data.type="MCS")
    
    checkTrue(is.list(dat))
    checkTrue(identical(names(dat$norm.factor.vec), names(dat$chip.list)))
    
    ## count ChIP reads around each binding site
    dat <- get.site.count(dat, window.size=250)
    
    checkTrue(!is.null(dat$site.count))
    checkTrue(identical(colnames(dat$site.count), names(dat$chip.list)))
    checkTrue(identical(nrow(dat$site.count), sum(sapply(dat$consensus.site, nrow))))
}
