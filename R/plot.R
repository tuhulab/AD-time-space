plotVarPart_internal <- function( obj, col=c(ggColorHue(ncol(obj)-1), "grey85"), label.angle=30, main="", ylab='', convertToPercent=TRUE, ylim,...){

  # convert to data.frame
  obj = as.data.frame(obj, check.names=FALSE)

  if( length(col) < ncol(obj) ){
    stop("Not enough colors specified by col")
  }

  # get gene name of each row
  obj$gene <- rownames(obj)

  # variable levels
  obj_level <-
    obj[,-(ncol(obj)-1):-ncol(obj)] %>%
    as.matrix() %>%
    colMedians(useNames = TRUE) %>%
    sort(decreasing = TRUE) %>%
    names() %>%
    c(colnames(obj)[ncol(obj)-1])


  # convert to data.frame for ggplot
  data <- reshape::melt(obj, id="gene")

  if( min(data$value) < 0 ){
    warning("Some values are less than zero")
  }

  if( convertToPercent ){
    data$value <- data$value * 100

    if( missing(ylim) ){
      ylim = c(0, 100)
    }

  }else{
    if( missing(ylim)){
      ylim = c(0, max(data$value))
    }
  }

  # change level
  data$variable <- data$variable %>% forcats::fct_relevel(obj_level)

  # add to pass R CMD check
  variable <- 1
  value <- 1

  # violin plot
  fig = ggplot(data=data, aes(x=variable, y=value)) +
    geom_violin( scale="width", aes(fill = factor(variable))) +
    ylab(ylab) + xlab('') + ylim(ylim) + theme_bw() +
    geom_boxplot(width=0.07, fill="grey", outlier.colour='black') +
    scale_fill_manual(values=col) +
    theme(legend.position="none") +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x = element_text(size  = 13,
                                     angle = label.angle,
                                     hjust = 1,
                                     vjust = 1))

  fig = fig + theme(text 		= element_text(colour="black"),
                    axis.text 	= element_text(colour="black"),
                    legend.text = element_text(colour="black"))

  if( main != ""){
    fig = fig + ggtitle( main ) + theme(plot.title = element_text(lineheight=.8, face="bold"))
  }

  return( fig )
}
