
username <- 'USERNAME'

### get Google cooper tracker
key         <- "FILEPATH" 
sheet_name  <- FILEPATH
gs_download <- function(key, sheet_name=NULL) {
  if (!is.null(sheet_name)) sheet <- paste0("&gid=", sheet_name) else sheet <- ""
  return(fread(paste0("URLPATH", sheet)))
}
df <- gs_download(key, sheet_name=sheet_name)


### compare to clean cooper refresh (download directly from GHDx Imports and Exports > Project export)
file_name <- "2020-06-25_cooper_export_for_project_373684"  # Change date to read in most recent GHDx Cooper download
refresh <- fread(paste0("FILEPATH")) 
refresh <- refresh[, names(refresh)[names(refresh) %in% names(df)], with=FALSE]

### append the new nids in GHDX Cooper not in Google tracker
refresh_new <- refresh[!nid %in% unique(df$nid)]
df <- rbind(df, refresh_new, fill=TRUE)

### what are changes that happened in the GHDx Cooper project not reflected in Google?
invisible(lapply(names(refresh), function(x) refresh[, (x) := as.character(get(x))]))
invisible(lapply(names(df), function(x) df[, (x) := as.character(get(x))]))
invisible(lapply(names(refresh), function(x) refresh[is.na(get(x)), (x) := ""]))
invisible(lapply(names(df), function(x) df[is.na(get(x)), (x) := ""]))
# identify the NIDs of sources where data type changed
cols <- c("microdata_status", "tabulations_status")
changes <- fsetdiff(refresh[, c("nid", cols), with=FALSE], df[, c("nid", cols), with=FALSE])[, nid]
# for each NID with a microdata/tabulation_status change, update the status
for (change in changes) {
  new_1 <- refresh[nid %in% change, cols[1], with=FALSE]
  new_2 <- refresh[nid %in% change, cols[2], with=FALSE]
  new <- ""
  # id if change is either in microdata_status or tabulation_status col
  if (new_1 != df[nid %in% change, cols[1], with=FALSE]) new <- paste(new, paste0("Change: ", new_1), sep=ifelse(nchar(new) >= 1, "; ", ""))
  if (new_2 != df[nid %in% change, cols[2], with=FALSE]) new <- paste(new, paste0("Change: ", new_2), sep=ifelse(nchar(new) >= 1, "; ", ""))
  # update the appropriate column for the specific NID
  if (new_1 != df[nid %in% change, cols[1], with=FALSE]) df[nid==change, (cols[1]) := new_1]
  if (new_2 != df[nid %in% change, cols[2], with=FALSE]) df[nid==change, (cols[2]) := new_2]
  # Flag that this NID changed data status in the sheet
  df[nid %in% change, check_updates := paste0(check_updates, paste0("NEW ", format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y-%m-%d"), " - ", new), sep=ifelse(nchar(check_updates) >= 1, "; ", ""))]  #"%Y-%m-%d %H:%M"
}

### check some easy things, like making sure all "Excluded" sources are not "" in exclusion_status
if (length(nrow(df[acceptance_status=="Excluded" & exclusion_status==""])) > 1) stop("Hang on! You have Excluded rows missing an exclusion_status. Fix in Google df and re-prep before uploading.")


### save files to upload back to Cooper and to replace in Google
## File for re-upload to Google:
# fwrite(df, paste0("FILEPATH"), row.names=FALSE)
fwrite(df, paste0("FILEPATH"), row.names=FALSE)

## File for re-upload to GHDx:
invisible(lapply(c("processing_status", "keywords", "secondary_data_type", "ownership"), function(x) df[, (x) := ""]))
df[acceptance_status %in% c("Accepted but want to replace", "Accepted", "Need more information to assess"), processing_status := "Assigned"]
df[acceptance_status %in% c("Accepted but want to replace", "Accepted", "Need more information to assess", "To be reviewed"), exclusion_status := ""]
df[acceptance_status %in% c("To be reviewed", "Excluded"), assigned := ""]
df[acceptance_status=="Excluded", add_comment := exclusion_status]
df[processing_status=="Assigned", assigned := username]

fwrite(df[, .(nid, title, assigned, additional_data_needs, add_comment, acceptance_status, exclusion_status, microdata_status, tabulations_status,
              processing_status, keywords, secondary_data_type, ownership)], 
       paste0("FILEPATH"), row.names=FALSE)

