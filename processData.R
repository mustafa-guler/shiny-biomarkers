filterData <- function(tcga, gores) {
    filtered <- tcga[.removeVers(tcga$gene_id) %in% gores$ensembl_gene_id, ]
    
    names <- unlist(filtered$gene_id)
    filtered <- as.data.frame(t(filtered[,-1]))
    colnames(filtered) <- names
    return(filtered)
}

addSurvivalData <- function(tcga, clinical) {
    sub_ids <- row.names(tcga)
    sub_ids <- sub_ids[!grepl("^go", sub_ids)]
    sub_ids <- sub("-.{3}$", "", sub_ids)
    tcga <- tcga[sub_ids %in% clinical$submitter_id,]
    sub_ids <- sub_ids[sub_ids %in% clinical$submitter_id]
    tcga$times <- vapply(
        sub_ids,
        function(x)
            return(clinical[clinical$submitter_id == sub("-.{3}$", "", x),]$years_to_death),
        FUN.VALUE = numeric(1)
    )
    tcga$is_dead <- vapply(
        sub_ids,
        function(x)
            return(clinical[clinical$submitter_id == sub("-.{3}$", "", x),]$vital_status == "dead"),
        FUN.VALUE = numeric(1)
    )
    tcga
}

trainingData <- function(dat, fraction) {
    sample(1:nrow(dat), round(nrow(dat) * fraction))
}

train <- function(dat, training) {
    rfsrc(Surv(times, is_dead) ~ ., dat[training, ], splitrule = "logrankscore")
}

test <- function(data, model, training) {
    predict.rfsrc(model, data[-training, ])
}

plotVIMP <- function(model) {
    plot(gg_vimp(model))
}

calcVIMP <- function(model) {
    vimps <- vimp(model)$importance
    vimps <- sort(vimps, decreasing = TRUE)
    vimps[1:8]
}

plotVariables <- function(model, vimps) {
    plot.variable(model, xvar.names = names(vimps), plots.per.page = 4, surv.type = "rel.freq")
}

plotPie <- function(clinical, target, titlename) {
    tab <- table(clinical[target])
    tab <- as.data.frame(round(tab / sum(tab) * 100, 2))
    colnames(tab)[1] <- target
    ggplot(tab, aes_string(x = NA, y = "Freq", fill=target)) + 
        geom_bar(width = 1, stat = "identity") + 
        coord_polar("y", start = 0) + 
        ggtitle(titlename) + 
        theme(
               axis.title = element_blank(),
               # plot.title = element_text(size = 32),
               axis.text.y = element_blank()
        )
}

plotKM <- function(dat, vimps) {
    times <- dat$times
    is_dead <- dat$is_dead
    dat <- dat[, names(vimps)]
    dat[] <- lapply(dat, function(x) ifelse(x < median(x), "LOW", "HIGH"))
    dat$surv_obj <- Surv(times, is_dead)
    
    par(mfrow = c(2,4))
    plots <- list()
    for(i in 1:6) {
        fit <- survfit(as.formula(paste("surv_obj~", names(vimps)[i])), data = dat)
        curr <- autoplot(fit, conf.int = FALSE) + ggtitle(names(vimps)[i])
        plots[[i]] <- curr
    }
    
    grid.arrange(grobs = as.list(plots), nrow = 2)
    

}

.removeVers <- function(ids) {
    return(
        vapply(
            unlist(ids),
            function(x)
                return(sub("\\..*", "", x)),
            FUN.VALUE = character(1)
        )
    )
}
