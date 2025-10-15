# Code Quality Check Script
# This script runs lintr on all R files and generates a quality report

library(lintr)

# Function to calculate code quality score
calculate_quality_score <- function(lint_results) {
  if (length(lint_results) == 0) {
    return(100)
  }
  
  # Count lints by severity
  total_lints <- length(lint_results)
  
  # Simple scoring: 100 - (number of issues)
  # Cap at 0 minimum
  score <- max(0, 100 - total_lints)
  
  return(score)
}

# Function to generate quality badge
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
    "https://img.shields.io/badge/code%%20quality-%d%%25-%s",
    score, color
  )
  
  return(badge_url)
}

# Main execution
cat("Running code quality checks...\n\n")

# Lint all R files
r_files <- list.files(
  path = c("R", "."),
  pattern = "\\.R$|\\.Rmd$",
  full.names = TRUE,
  recursive = FALSE
)

all_lints <- list()
for (file in r_files) {
  cat(sprintf("Checking %s...\n", file))
  lints <- lint(file)
  if (length(lints) > 0) {
    all_lints[[file]] <- lints
    print(lints)
  }
}

# Calculate overall score
total_lints <- sum(sapply(all_lints, length))
score <- calculate_quality_score(unlist(all_lints))

# Generate report
cat("\n=== Code Quality Report ===\n")
cat(sprintf("Total files checked: %d\n", length(r_files)))
cat(sprintf("Files with issues: %d\n", length(all_lints)))
cat(sprintf("Total issues found: %d\n", total_lints))
cat(sprintf("Code Quality Score: %d/100\n", score))

# Generate badge
badge_url <- generate_badge(score)
cat(sprintf("\nBadge URL: %s\n", badge_url))
cat(sprintf("Badge Markdown: ![Code Quality](%s)\n", badge_url))

# Save report
report <- list(
  timestamp = Sys.time(),
  score = score,
  total_files = length(r_files),
  files_with_issues = length(all_lints),
  total_issues = total_lints,
  badge_url = badge_url,
  details = all_lints
)

saveRDS(report, "code_quality_report.rds")
write(sprintf("Code Quality Score: %d/100", score), "code_quality_score.txt")

cat("\nReport saved to: code_quality_report.rds\n")
cat("Score saved to: code_quality_score.txt\n")
