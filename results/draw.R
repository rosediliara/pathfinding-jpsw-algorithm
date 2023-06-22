draw_grid <- function(gridfile)
{
    tmp <- read.csv(gridfile, header=FALSE, sep="\n", stringsAsFactors=FALSE)
    tmp <- unlist(tmp)
    width <- as.numeric(gsub("width", "", tmp[grep("width", tmp)]))
    height <- as.numeric(gsub("height", "", tmp[grep("height", tmp)]))
    tmp <- tmp[grep("map|height|width|octile", tmp, invert=TRUE)]
    plot(NA, main=basename(gridfile), xlim=c(0, width), ylim=c(0, height), 
        ylab=height, xlab=width, xaxt="n", yaxt="n")

    for(i in seq(1, length(tmp)))
    {
        obs <- strsplit(tmp[i], "")
        xvals <- which(unlist(obs) == "@") 
        yvals <- rep(i, length(xvals))
        # R is 1-indexed, grids are 0-indexed. So we apply a -1 offset
        points(xvals-1, yvals-1, pch=".", col="black", bg="black", cex=0.5)
    }
}

draw_plan <- function(planfile, mapwidth, path_color="blue", path_lwd=2, text_labels=FALSE)
{
    plan <- unlist(read.csv(planfile, header=FALSE, sep="\n", stringsAsFactors=FALSE))
    pbegin <- grep("^p", plan)
    for(i in seq(1, length(pbegin)))
    {
        s_index <- pbegin[i]+1
        if(s_index > length(plan)) { break; }

        t_index <- length(plan)
        if(i < length(pbegin))
        {
            t_index <- pbegin[i+1] - 1;
        }
        if(t_index <= s_index) { next; }

        xvals <- rep(0, (t_index - s_index))
        yvals <- rep(0, (t_index - s_index))

        print(paste(paste("plan", i), paste(plan[s_index], plan[t_index])))
        for(j in seq(0, (t_index-s_index), by=1))
        {
            xy_id <- as.numeric(unlist(strsplit(plan[j+s_index], " "))[1])
            #print(paste(j, xy_id))
            x <- xy_id %% mapwidth
            y <- floor(xy_id / mapwidth)
            #print(paste(x, y))
            xvals[j+1] <- x
            yvals[j+1] <- y
        }
        #print(xvals)
        if(length(xvals > 1))
        {
            points(xvals, yvals, pch=".", col=path_color, lwd=path_lwd, cex=0.5)
            points(xvals[1], yvals[1], col=path_color, bg=path_color, pch=21)
            points(xvals[length(xvals)], yvals[length(yvals)], col=path_color, 
                bg=path_color, pch=22)
            if(text_labels)
            {
                text(xvals[1], yvals[1], labels=paste("S", i, sep=""), col=path_color, pos=1)
                text(xvals[length(xvals)], yvals[length(yvals)], 
                    labels=paste("T", i, sep=""), col=path_color, pos=1)
            }
            print(paste(xvals[1], yvals[1]))
        }
    }
}
