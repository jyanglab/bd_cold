
############
library("ggplot2")
library("plyr")
res <- read.csv("data/qpcr_rnaseq_cv.csv")
# Error bars represent standard error of the mean
part1 <- ddply(res, .(Assay, temp), summarize,
               exp = mean(exp.x),
               err = sd(exp.x))
part1$temp <- as.factor(part1$temp)
# Error bars represent standard error of the mean
part2 <- ddply(res, .(Assay, temp), summarize,
               exp = mean(exp.y),
               err = sd(exp.y))
part2$temp <- as.factor(part2$temp)

p1 <- ggplot(part1, aes(x=Assay, y=exp, fill=temp)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=exp-err, ymax=exp+err),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("RPKM Value") +
    ggtitle("RNA-seq") +
    labs(fill="Temp (°C)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

p2 <- ggplot(part2, aes(x=Assay, y=exp, fill=temp)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=exp-err, ymax=exp+err),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    xlab("") +
    ylab("Relative Expression") +
    ggtitle("qPCR") +
    labs(fill="Temp (°C)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

pdf("graphs/Fig_S1.pdf", width=10, height=5)
multiplot(p2, p1, cols=2)
dev.off()

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}