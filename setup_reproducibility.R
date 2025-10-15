# Setup script for reproducibility
# This script initializes renv and captures all package dependencies

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv (this will create renv.lock)
renv::init()

# Capture session info for documentation
session_info <- sessionInfo()
saveRDS(session_info, "session_info.rds")

cat("Session Info:\n")
print(session_info)

cat("\nR Environment initialized successfully!\n")
cat("renv.lock file created with all package dependencies.\n")
cat("Run 'renv::restore()' to install the exact package versions.\n")
