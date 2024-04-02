profile <- Sys.getenv("QUARTO_PROFILE")

if(profile == "public"){

# cp -r videos $(QUARTO_PROJECT_OUTPUT_DIR)/src/
#
# This is intended as a Quarto post-render script to move local
# video files into the correct website directory.
#
# This is a workaround to an apparent Quarto bug and should not be 
# necessary when Quarto > 1.4 is available to us.
#
# see https://github.com/quarto-dev/quarto-cli/issues/3892

here::i_am("src/move_videos.R")
library(here)

# this is supposed to work, but doesn't seem to do
# output_dir <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR")

# maybe needs a newer Quarto version??
# https://quarto.org/docs/projects/scripts.html#pre-and-post-render

# that environment variable is not set here:
# envvars <- Sys.getenv()
# print(envvars[grepl("QUARTO", names(envvars))])

# instead:
output_dir <- file.path("docs", profile)

source_dir <- here("videos")
target_dir <- file.path(output_dir, "src/")

cat("Copying from", here("videos"), "to", here(target_dir), "...")

stopifnot(dir.exists(output_dir))
if(!dir.exists(target_dir)) dir.create(target_dir)

file.copy(here("videos"), target_dir, recursive = TRUE)

}
