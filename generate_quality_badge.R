#!/usr/bin/env Rscript
# Generate Code Quality Badge for README
#
# This script calculates a code quality score and generates badge markdown
# that can be added to the README.md file.

# Function to count issues in R files
count_code_issues <- function() {
  r_files <- list.files(
    path = c("R", "."),
    pattern = "\\.R$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  # Simple heuristics for code quality
  issues <- 0
  total_lines <- 0
  
  for (file in r_files) {
    lines <- readLines(file, warn = FALSE)
    total_lines <- total_lines + length(lines)
    
    # Count potential issues
    # 1. Lines too long (>120 chars)
    long_lines <- sum(nchar(lines) > 120)
    
    # 2. TODO/FIXME comments (not necessarily bad, but worth tracking)
    todo_comments <- sum(grepl("TODO|FIXME", lines, ignore.case = TRUE))
    
    # 3. Missing documentation (lines starting with function without roxygen)
    func_lines <- grep("^[[:alnum:]_]+ <- function", lines)
    undocumented <- 0
    for (i in func_lines) {
      if (i > 1 && !grepl("^#'", lines[i-1])) {
        undocumented <- undocumented + 1
      }
    }
    
    issues <- issues + long_lines * 0.5 + todo_comments * 0.2 + undocumented * 2
  }
  
  list(
    total_files = length(r_files),
    total_lines = total_lines,
    estimated_issues = round(issues),
    score = max(0, min(100, 100 - issues))
  )
}

# Calculate quality metrics
cat("Calculating code quality metrics...\n")
metrics <- count_code_issues()

cat("\n=== Code Quality Metrics ===\n")
cat(sprintf("Total R files: %d\n", metrics$total_files))
cat(sprintf("Total lines of code: %d\n", metrics$total_lines))
cat(sprintf("Estimated issues: %d\n", metrics$estimated_issues))
cat(sprintf("Code Quality Score: %.0f/100\n", metrics$score))

# Generate badge
generate_badge <- function(score) {
  color <- if (score >= 90) {
    "brightgreen"
  } else if (score >= 75) {
    "green"
  } else if (score >= 60) {
    "yellow"
  } else if (score >= 40) {
    "orange"
  } else {
    "red"
  }
  
  badge_url <- sprintf(
    "https://img.shields.io/badge/code%%20quality-%.0f%%25-%s",
    score, color
  )
  
  return(badge_url)
}

badge_url <- generate_badge(metrics$score)

cat("\n=== Badge for README ===\n")
cat(sprintf("![Code Quality](%s)\n", badge_url))

cat("\nAdd this badge to your README.md file!\n")

# Save metrics
saveRDS(metrics, "code_quality_metrics.rds")
cat("\nMetrics saved to: code_quality_metrics.rds\n")
