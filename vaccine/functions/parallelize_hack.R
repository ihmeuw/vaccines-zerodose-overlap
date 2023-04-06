#Note that if launching from RStudio, works best if RStudio session is launched on <IMAGE> image

parallelize_hack <- function(lv,
                             save_objs,
                             script,
                             script_dir,
                             prefix,
                             output_log_dir,
                             error_log_dir,
                             project,
                             queue,
                             memory,
                             cores,
                             time,
                             sing_image,
                             raking=F){


# Create filename using time stamp to save 'save objs' for reloading
tmp_dir <- file.path("<FILEPATH>")
# Append random string on the end to avoid overlapping filenames for
# rapidly-submitted jobs
frand <- paste0(user, "_", gsub("-|:| ", "_", Sys.time()), sample(1:100000, 1))
fname <- file.path(tmp_dir, paste0(frand, ".RData"))

###
lv <- data.table(expand.grid(lv, stringsAsFactors = F))
lv$jobname <- paste0(prefix, "_", apply(lv, 1, paste, collapse = "_"))

# Ensure all whitespace collapsed
lv[, jobname := gsub(" ", "", jobname, fixed = TRUE)]

values_to_save <- c("lv")
if (!is.null(save_objs)) values_to_save <- c(values_to_save, save_objs)
save(
  file = fname,
  list = values_to_save
)


script_file <- paste0(script_dir, "/", script)



# Sbatch over lv rows
for (i in 1:nrow(lv)) {
  job_name <- lv[i, jobname]
  
  if(raking==T) memory <- load_profiling(vax='raking',draws=500,reg=lv[i,region])[['mem']]

sbatch <- paste0('sbatch -o ',output_log_dir,' -e ', error_log_dir,
               ' -A ',project,' -J ',job_name,' -p ',queue,
               ' -C archive --mem ', memory,'G -c ',cores,' -t ',time, ' --wrap \"<IMAGE> ',
               paste(sing_image,'-s',script_file,frand,i,'"'))

returned <- system(sbatch, intern = T)
message(returned)
#system returns e.g. "Submitted batch job 13373"
job_id <- as.numeric(gsub("Submitted batch job ", "", returned))

lv[i, jobid := job_id]
lv[i, the_qsub := sbatch]
}
return(list(lv, fname, sbatch))
}
