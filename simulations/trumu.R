#function (attributes, population.data, data.method, time.interval.method = "linear", 
#    ...)
    
attributes
population.data = m.plot.data
data.method = "list"
time.interval.method = "linear"
#function (attributes, population.data, data.method, time.interval.method = "linear", 
#    ...) 
#{
    if (!(time.interval.method %in% c("linear", "equal"))) {
        base::stop("time.interval.method must be one of \"linear\" or \"equal\".")
    }

    if (base::ncol(attributes) == 3) {
        parent.list <- attributes[, c(1, 2)]
        color.list <- attributes[, 3]
    }
    
num.ssp = 10
    
    
####     else{ if (base::ncol(attributes) == 2) {
####        parent.list <- attributes[, c(1, 2)]
####        color.list <- NA
####    }else {
####        base::stop("attributes must have 2 or 3 columns containing OTU names, parents and colors (optional), respectively.")
####    }}
####    num.ssp <- base::nrow(parent.list)

####    if (base::all(base::is.na(color.list))) {
####        if (num.ssp > 12) {
####            color.list <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, 
####                "Paired")))(num.ssp)
####        }
####        else {
####            color.list <- RColorBrewer::brewer.pal(num.ssp, "Paired")
####        }
####    }
####    else {
####        if (base::length(color.list) != num.ssp) {
####            base::stop("Length of color.list doesn't match the number of OTUs.")
####        }
####    }
####    if (base::length(unique(color.list)) < length(color.list)) {
####        base::warning("Some OTUs have the same color in attributes.")
####    }
####    if (data.method == "list") {
####        if (base::ncol(population.data) != 3) {
####            base::stop("population.data must have 3 colomns containing OTU names, time/generation and abundances, respectively.")
####        }
####        base::colnames(population.data) <- c("names", "times", 
####            "abundances")
####        species.name.pop <- base::unique(base::as.character(population.data$names))
####        if (base::nrow(base::unique(population.data)) != base::nrow(population.data)) {
####            base::stop("There are some repeated rows in population.data.")
####        }
####        if (base::any(base::is.na(base::as.numeric(population.data$times)))) {
####            base::stop("Second column of population.data must be numbers (time steps or generations).")
####        }
####    }
####    else if (data.method == "table") {
####        if (base::length(base::rownames(population.data)) < base::nrow(population.data)) {
####            base::stop("rownames of population.data must be OUT names.")
####        }
####        else {
####            species.name.pop <- base::rownames(population.data)
####        }
####        if (base::length(base::which(base::table(species.name.pop) != 
####            1)) > 0) {
####            base::stop("There are some repeated entry in the rownames of population.data.")
####        }
####        if (base::length(base::colnames(population.data)) < base::ncol(population.data)) {
####            base::stop("colnames of population.data must be time steps or generations.")
####        }
####        else {
####            time.points.input <- base::as.numeric(base::colnames(population.data))
####        }
####        if (base::any(base::is.na(time.points.input))) {
####            base::stop("colnames of population.data must be numbers (time steps or generations).")
####        }
####        if (base::length(base::which(base::table(base::as.numeric(base::colnames(population.data))) != 
####            1)) > 0) {
####            base::stop("There are some repeated entry in the colnames of population.data.")
####        }


####        num.step <- base::length(time.points.input)
####        pop.data.name <- c()
####        pop.data.time <- c()
####        pop.data.abundace <- c()
####        for (t in 1:num.step) {
####            time <- base::as.numeric(time.points.input[t])
####            live.id <- base::which(population.data[, t] > 0)
####            pop.data.name <- base::c(pop.data.name, base::rownames(population.data)[live.id])
####            pop.data.time <- base::c(pop.data.time, base::rep(time, 
####                base::length(live.id)))
####            pop.data.abundace <- base::c(pop.data.abundace, base::as.numeric(population.data[live.id, 
####                t]))
####            if (base::any(pop.data.abundace < 0) || base::any(base::is.na(pop.data.abundace))) {
####                base::stop("Abundance of species must be numbers greater than or equal to zero.")
####            }
####        }
####        population.data <- base::data.frame(names = pop.data.name, 
####            times = pop.data.time, abundances = pop.data.abundace)
####    }


####    else {
####        base::stop("data.method must be one of \"list\" or \"table\".")
####    }

#####    f <- FALSE
#####    for (i in 1:num.ssp) {
#####        if (!base::any(population.data$abundances[which(population.data$names == 
#####            species.name.pop[i])] != 0)) {
#####            f <- TRUE
#####            population.data <- population.data[base::which(population.data$names != 
#####                species.name.pop[i]), ]
#####        }
#####    }
#####    if (f) {
#####        warning("Abundances of some OTUs are always zero and they will be ignored.", 
#####            immediate. = TRUE)
#####        species.name.pop <- base::unique(base::as.character(population.data$names))
#####    }

    parent.list <- base::as.data.frame(parent.list)
#    base::colnames(parent.list) <- c("names", "parent")
    species.name.par <- base::unique(base::as.character(parent.list$names))
#####    if (base::length(base::intersect(species.name.pop, species.name.par)) == 
#####        0) {
#####        base::stop("OTUs names in attributes (first column) and population.data (rownames) do not match.")
#####    }
#####    if (!base::all(species.name.pop %in% species.name.par)) {
#####        base::stop("Parents of some OTUs are not specified. Please Specify all the parents and put NA for those OTUs which have not a parent in the present data.")
#####    }
#####    if (base::length(base::which(base::table(base::as.character(parent.list$names)) != 
#####        1)) > 0) {
#####        base::stop("There are some repeated entry in the OTUs of parent list.")
#####    }
#####    if (!base::all(base::unique(base::as.character(parent.list$parent)) %in% 
#####        base::c(NA, base::unique(base::as.character(parent.list$names))))) {
#####        base::stop("Parents of some OTUs are not from the OTU list. Please put NA as parent of those OTUs which their parents are not present in the population dynamics. Please check attributes.")
#####    }
#####    if (base::all(!base::is.na(parent.list$parent))) {
#####        base::stop("Parent of at least one OUTs (ancesstor of other) must be NA. Please check attributes.")
#####    }
#####    if (!base::all(species.name.par %in% species.name.pop)) {

#####        not.present <- base::as.character(parent.list$names[base::which(!(species.name.par %in% 
#####            species.name.pop))])
#####        base::warning("These OTUs are not present in the dynamics and are ignored:", 
#####            immediate. = TRUE)
#####        for (sp in not.present) {
#####            base::message(paste(sp, " "), appendLF = FALSE)
#####        }
#####        base::message("")
#####        for (sp in not.present) {
#####            parent.list$parent[base::which(parent.list$parent == 
#####                sp)] <- parent.list$parent[base::which(parent.list$names == 
#####                sp)]
#####        }
#####        parent.list <- parent.list[-base::which(!(species.name.par %in% 
#####            species.name.pop)), ]
#####        color.list <- color.list[-base::which(!(species.name.par %in% 
#####            species.name.pop))]
#####    }

    num.ssp <- base::length(species.name.pop)
    for (i in 1:num.ssp) {
        if (!(parent.list$parent[i] %in% parent.list$names)) {
            parent.list$parent[i] <- NA
        }
    }
#    name.color <- base::data.frame(name = base::as.character(parent.list$names), 
#        color = color.list)
#    tempid.name <- base::data.frame(temp.ids = 1:num.ssp, name = species.name.pop)
#    parent.tempid.list <- base::data.frame(temp.id = base::as.numeric(tempid.name$temp.ids[base::match(x = base::as.character(parent.list$names), 
#        base::as.character(tempid.name$name))]), parent.id = base::as.numeric(tempid.name$temp.ids[base::match(x = base::as.character(parent.list$parent), 
#        base::as.character(tempid.name$name))]))
#    orphan.list <- base::as.numeric(parent.tempid.list$temp.id[base::which(base::is.na(base::as.numeric(parent.tempid.list$parent.id)))])
#    tempid.order <- orphan.list
#    for (ssp in orphan.list) {
#        base::environment(FillOrderList) <- environment()
#        tempid.order <- FillOrderList(ssp, tempid.order)
#    }
#    temp.ids <- base::match(x = population.data$names, species.name.pop)
#    pop.dyn <- base::as.data.frame(cbind(temp.ids, population.data[, 
#        -1]))
#    time.points <- base::sort(base::unique(pop.dyn$times), decreasing = FALSE)
#    start <- base::min(time.points)
#    num.step <- base::length(time.points)
#    emergence.time <- c()
#    extinction.time <- c()
#    for (i in 1:num.ssp) {
#        time.point.sp <- pop.dyn$times[base::intersect(base::which(pop.dyn$temp.ids == 
#            i), base::which(pop.dyn$abundances > 0))]
#        if (base::length(base::unique(time.point.sp)) < base::length(time.point.sp)) {
#            base::stop("There is at least one repeated time point for at least one OTU. Please check population.data.")
#        }
#        emergence.time <- base::c(emergence.time, min(time.point.sp))
#        extinction.time <- base::c(extinction.time, max(time.point.sp))
#    }
#    temp.ids <- 1:num.ssp
#    species.info <- base::as.data.frame(base::cbind(temp.ids, 
#        emergence.time, extinction.time))
#    sp.t.lost <- c()
#    t.lost.total <- c()
#    for (i in 1:num.ssp) {
#        time.seri <- time.points[base::which(time.points == species.info$emergence.time[i]):base::which(time.points == 
#            species.info$extinction.time[i])]
#        time.seri.sp <- base::sort(pop.dyn$times[base::which(pop.dyn$temp.ids == 
#            i)], decreasing = FALSE)
#        t.lost <- time.seri[base::which(!(time.seri %in% time.seri.sp))]
#        if (base::length(t.lost) > 0) {
#            t.lost.total <- base::c(t.lost.total, list(t.lost))
#            sp.t.lost <- base::c(sp.t.lost, i)
#        }
#    }
#    if (base::length(sp.t.lost) > 0) {
#        base::message("Please provide these missing data points:")
#        for (i in 1:base::length(sp.t.lost)) {
#            base::message(paste("OTU", sp.t.lost[i], "lacks these time points:", 
#                t.lost.total[i], sep = " "))
#        }
#        base::stop("Some time points are missing for at least one OTUs or at least one zero abundance occures between two other non-zeros.")
#    }
#    for (i in 1:num.ssp) {
#        name <- base::as.character(tempid.name$name[i])
#        par <- base::which(base::as.character(tempid.name$name) == 
#            base::as.character(parent.list$parent[base::which(parent.list$names == 
#                name)]))
#        if (base::length(par) > 0) {
#            if (species.info$emergence.time[i] < species.info$emergence.time[par]) {
#                base::stop("Emergence time of at least one OTU is smaller than its parent. Please check attributes or population.data.")
#            }
#        }
#        if (!(base::is.na(as.character(parent.list$parent[i])))) {
#            if (base::as.character(parent.list$names[i]) == base::as.character(parent.list$parent[i])) {
#                base::stop("Name of at least one OTU is equal to its parent. Please check attributes.")
#            }
#        }
#    }
#    sp.id <- 1:num.ssp
#    species.info <- base::as.data.frame(base::cbind(sp.id, species.info[base::match(x = tempid.order, 
#        species.info$temp.ids), ]))
#    temp.to.real.id <- base::as.data.frame(base::cbind(species.info$temp.ids, 
#        species.info$sp.id))
#    base::colnames(temp.to.real.id) <- c("temp.ids", "sp.id")
#    temp.to.real.id <- temp.to.real.id[base::order(temp.to.real.id$temp.ids), 
#        ]
#    color.list.sorted <- base::as.character(name.color$color[base::match(x = base::as.character(tempid.name$name[base::match(x = tempid.order, 
#        tempid.name$temp.ids)]), base::as.character(name.color$name))])
#    sp.id <- base::rep(NA, base::dim(population.data)[1])
#    pop.dyn <- base::cbind(sp.id, pop.dyn)
#    for (id in 1:num.ssp) {
#        pop.dyn$sp.id[base::which(pop.dyn$temp.ids == species.info$temp.ids[id])] <- species.info$sp.id[id]
#    }
#    pop.dyn <- pop.dyn[base::order(pop.dyn$times), ]
#    pop.dyn <- pop.dyn[base::order(pop.dyn$sp.id), ]
#    species.info <- species.info[-2]
#    pop.dyn <- pop.dyn[-2]
#    for (i in 1:num.step) {
#        times <- base::which(pop.dyn$times == time.points[i])
#        tot.pop <- base::sum(pop.dyn$abundances[times])
#        pop.dyn$abundances[times] <- pop.dyn$abundances[times]/tot.pop
#    }
#    p.list <- base::rep(NA, num.ssp)
#    s.list <- base::rep(NA, num.ssp)
#    for (i in 1:num.ssp) {
#        s.list[i] <- temp.to.real.id$sp.id[base::which(species.name.pop == 
#            base::as.character(parent.list$names)[i])]
#        p <- temp.to.real.id$sp.id[base::which(species.name.pop == 
#            base::as.character(parent.list$parent)[i])]
#        if (base::length(p) != 0) 
#            p.list[i] = p
#    }
#    parent.id.list <- base::as.data.frame(base::cbind(s.list, 
#        p.list))
#    parent.id.list <- parent.id.list[base::order(parent.id.list$s.list), 
#        ]
#    ssp.children.list = base::rep(0, num.ssp)
#    for (i in 1:num.ssp) {
#        ssp.children.list[i] <- base::list(base::which(parent.id.list$p.list == 
#            i))
#    }
#    ssp.parent <- parent.id.list$p.list
#    ssp.parent[base::is.na(ssp.parent)] <- 0
#    polygon.pos = base::matrix(list(c()), nrow = num.step, ncol = num.ssp)
#    total.pop = base::matrix(0, nrow = num.step, ncol = num.ssp)


#    base::cat("Processing time steps:\n")
    for (t in 1:num.step) {
        base::cat(time.points[t], " |")
        live.ssp <- 1:num.ssp
        orphan.list <- base::which(ssp.parent == 0)
        if (t > 1) {
            min.pos <- c()
            for (ssp in live.ssp) {
                if (!base::is.null(polygon.pos[t - 1, ssp][[1]])) {
                  min.pos <- c(min.pos, base::min(polygon.pos[t - 
                    1, ssp][[1]]))
                }
                else {
                  min.pos <- c(min.pos, 2)
                }
            }
            live.ssp <- live.ssp[base::order(min.pos, decreasing = FALSE)]
        }
        live.child <- base::rep(base::list(c()), num.ssp)
        for (ssp in live.ssp) {
            base::environment(FindLiveSuccessors) <- base::environment()
            live.child <- FindLiveSuccessors(ssp, live.child)
        }
        start.list <- c()
        for (ssp in live.ssp) {
            ssp.temp <- ssp
            while (!(ssp.parent[ssp.temp] %in% live.ssp) && ssp.parent[ssp.temp] > 
                0) {
                ssp.temp <- ssp.parent[ssp.temp]
            }
            if (ssp.temp %in% orphan.list) {
                start.list <- base::c(start.list, ssp)
            }
        }
        for (ssp in start.list) {
            base::environment(FindTotalpop) <- base::environment()
            total.pop <- FindTotalpop(ssp, total.pop)
        }
        list.poly <- list(0, polygon.pos)
        for (ssp in start.list) {
            list.poly[[2]][t, ssp][[1]] <- base::c(list.poly[[2]][t, 
                ssp][[1]], list.poly[[1]])
            base::environment(FindPolygonPos) <- base::environment()
            list.poly <- FindPolygonPos(ssp, list.poly)
        }
        polygon.pos <- list.poly[[2]]
    }
    
    
    
#    base::cat("\n")
#    polygon.pos.temp <- polygon.pos
#    polygon.pos = base::matrix(base::list(c()), nrow = num.step, 
#        ncol = num.ssp)
#    for (i in 1:num.ssp) {
#        t <- 1
#        if (!(base::all(!base::is.null(polygon.pos.temp[t, i][[1]]), 
#            (polygon.pos.temp[t, i][[1]][1] == polygon.pos.temp[t, 
#                i][[1]][2]), (polygon.pos.temp[t + 1, i][[1]][1] == 
#                polygon.pos.temp[t + 1, i][[1]][2])))) {
#            polygon.pos[t, i][[1]] <- polygon.pos.temp[t, i][[1]]
#        }
#        for (t in 2:(num.step - 1)) {
#            if (!(base::all(!base::is.null(polygon.pos.temp[t, 
#                i][[1]]), (polygon.pos.temp[t, i][[1]][1] == 
#                polygon.pos.temp[t, i][[1]][2]), (base::any(base::all((polygon.pos.temp[t + 
#                1, i][[1]][1] == polygon.pos.temp[t + 1, i][[1]][2]), 
#                (polygon.pos.temp[t - 1, i][[1]][1] == polygon.pos.temp[t - 
#                  1, i][[1]][2]))))))) {
#                polygon.pos[t, i][[1]] <- polygon.pos.temp[t, 
#                  i][[1]]
#            }
#        }
#        t <- num.step
#        if (!(base::all(!base::is.null(polygon.pos.temp[t, i][[1]]), 
#            (polygon.pos.temp[t, i][[1]][1] == polygon.pos.temp[t, 
#                i][[1]][2]), (polygon.pos.temp[t - 1, i][[1]][1] == 
#                polygon.pos.temp[t - 1, i][[1]][2])))) {
#            polygon.pos[t, i][[1]] <- polygon.pos.temp[t, i][[1]]
#        }
#    }
#    first.time <- base::rep(0, num.ssp)
#    last.time <- base::rep(0, num.ssp)
#    for (ssp in 1:num.ssp) {
#        first.time[ssp] <- base::which(time.points == base::min(pop.dyn$times[base::intersect(base::which(pop.dyn$sp.id == 
#            ssp), base::which(pop.dyn$abundances != 0))])) - 
#            1
#        last.time[ssp] <- base::which(time.points == base::max(pop.dyn$times[base::intersect(base::which(pop.dyn$sp.id == 
#            ssp), base::which(pop.dyn$abundances != 0))])) + 
#            1
#    }
#    first.time[first.time == 0] <- 1
#    last.time[last.time == num.step + 1] <- num.step
#    RelativeAbundance = base::c(0, 1)
#    Generation = base::c(0, num.step)
#    if (time.interval.method == "linear") {
#        x.lim.last <- time.points[num.step]
#        x.lim.first <- time.points[1]
#    }
#    else {
#        x.lim.last <- num.step
#        x.lim.first <- 1
#    }
#    graphics::plot(Generation, RelativeAbundance, xlim = c(x.lim.first, 
#        x.lim.last), ylim = c(0, 1), xaxs = "i", yaxs = "i", 
#        col = "white", ...)
#    orphan.list.main <- parent.id.list$s.list[base::which(base::is.na(parent.id.list$p.list))]
#    for (ssp in orphan.list.main) {
#        base::environment(PlotPopulationDynamic) <- base::environment()
#        PlotPopulationDynamic(ssp, color.list.sorted)
#    }
#    color.list.final <- base::as.data.frame(base::cbind(base::as.character(tempid.name$name[base::order(temp.to.real.id$sp.id)]), 
#        color.list.sorted))
#    base::colnames(color.list.final) <- c("name", "color")
#    base::return(color.list.final)
#}

