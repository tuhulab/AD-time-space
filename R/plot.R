#' Plot Variance Partition Results
#'
#' Internal function to create violin plots showing the proportion of variance
#' explained by different factors in variance partition analysis.
#'
#' @param obj A data frame or matrix containing variance partition results where
#'   rows are genes and columns are variance components.
#' @param col Character vector of colors for each variance component. Default uses
#'   ggColorHue with grey85 for residuals.
#' @param label.angle Numeric value specifying the angle of x-axis labels in degrees.
#'   Default is 30.
#' @param main Character string for the plot title. Default is empty string.
#' @param ylab Character string for the y-axis label. Default is empty string.
#' @param convertToPercent Logical indicating whether to convert values to percentages.
#'   Default is TRUE.
#' @param ylim Numeric vector of length 2 specifying y-axis limits. If not provided,
#'   automatically set based on data range.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A ggplot2 object showing violin plots with boxplots overlay for each
#'   variance component.
#'
#' @details
#' This function creates publication-quality violin plots to visualize variance
#' partition results. Variance components are ordered by median value (descending)
#' and displayed with violin plots overlaid with boxplots. Values can be displayed
#' as proportions (0-1) or percentages (0-100).
#'
#' @examples
#' \dontrun{
#' # Assuming vp_results contains variance partition output
#' plotVarPart_internal(
#'   obj = vp_results,
#'   main = "Variance Partition",
#'   ylab = "Fraction of Variance Explained (%)"
#' )
#' }
#'
#' @import ggplot2
#' @importFrom reshape melt
#' @importFrom matrixStats colMedians
#' @export
plotVarPart_internal <- function(
    obj, 
    col = c(ggColorHue(ncol(obj) - 1), "grey85"), 
    label.angle = 30, 
    main = "", 
    ylab = "", 
    convertToPercent = TRUE, 
    ylim,
    ...
  ) {
  
  # Input validation
  if (!is.data.frame(obj) && !is.matrix(obj)) {
    stop("obj must be a data frame or matrix")
  }
  
  # Convert to data.frame
  obj <- as.data.frame(obj, check.names = FALSE)
  
  # Check color vector length
  if (length(col) < ncol(obj)) {
    stop("Not enough colors specified by col. Need at least ", ncol(obj), " colors.")
  }
  
  # Add gene identifiers
  obj$gene <- rownames(obj)
  
  # Order variance components by median value (descending)
  obj_level <-
    obj[, -(ncol(obj) - 1):-ncol(obj)] %>%
    as.matrix() %>%
    colMedians(useNames = TRUE) %>%
    sort(decreasing = TRUE) %>%
    names() %>%
    c(colnames(obj)[ncol(obj) - 1])  # Add residual component
  
  # Convert to long format for ggplot
  data <- reshape::melt(obj, id = "gene")
  
  # Validate data values
  if (min(data$value, na.rm = TRUE) < 0) {
    warning("Some values are less than zero. This may indicate issues with variance partition.")
  }
  
  # Convert to percentage if requested
  if (convertToPercent) {
    data$value <- data$value * 100
    
    if (missing(ylim)) {
      ylim <- c(0, 100)
    }
  } else {
    if (missing(ylim)) {
      ylim <- c(0, max(data$value, na.rm = TRUE))
    }
  }
  
  # Reorder factor levels based on median
  data$variable <- data$variable %>% forcats::fct_relevel(obj_level)
  
  # Suppress NOTE in R CMD check
  variable <- value <- NULL
  
  # Create violin plot with boxplot overlay
  fig <- ggplot(data = data, aes(x = variable, y = value)) +
    geom_violin(scale = "width", aes(fill = factor(variable))) +
    ylab(ylab) + 
    xlab("") + 
    ylim(ylim) + 
    theme_bw() +
    geom_boxplot(width = 0.07, fill = "grey", outlier.colour = "black") +
    scale_fill_manual(values = col) +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(
      axis.text.x = element_text(
        size = 13,
        angle = label.angle,
        hjust = 1,
        vjust = 1
      )
    )
  
  # Apply consistent text styling
  fig <- fig + theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    legend.text = element_text(colour = "black")
  )
  
  # Add title if provided
  if (main != "") {
    fig <- fig + 
      ggtitle(main) + 
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
  }
  
  return(fig)
}

