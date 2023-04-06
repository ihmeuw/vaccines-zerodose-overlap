## get_sp_hierarchy ######################################################

#' @title Get sp hierarchy
#' @description Pull the hierarchy of admin codes
#'
#' \code{get_sp_hierarchy()} takes the admin2 shapefile and pulls a list of
#' admin level codes to be used in merges, etc.
#'
#' @param shapefile_version character string indicating version of shapefile to pull
#'
#' @return a list of 3 data tables (adm0, adm1, adm2). Each data table includes
#'  the location codes of its administrative level and columns for the higher-
#'  order administrative levels as well.
#' @examples
#'

get_sp_hierarchy <- function(shapefile_version = 'current') {

  admin_level <- 2

  # Use most up to date admin2 files - in rds format for faster loading
  admin_shp <- rgdal::readOGR(get_admin_shapefile(admin_level, version = shapefile_version),
                              stringsAsFactors = FALSE, GDAL1_integer64_policy=TRUE)

  # Pull admin codes
  admin_data <- copy(data.table(admin_shp@data))
  admin_codes <- lapply(0:2, function(aa) {
    col_list <- sapply(0:aa, function(adm) {
      return(c(paste0("ADM", adm, "_CODE"),
               paste0("ADM", adm, "_NAME")))
      })
    col_list <- as.vector(col_list)
    data_subset <- subset(admin_data, select = col_list)
    data_subset <- unique(data_subset)
    return(data_subset)
  })

  names(admin_codes) <- c("ADM0", "ADM1", "ADM2")

  return(admin_codes)

}

#' @title Merge shape hierarchy
#' @description Takes a data frame or data table with a column of administrative codes
#' and merges on the accompanying admin code names, as well as any higher
#' order administrative codes within which that admin code is nested
#'
#' @param df Data frame or data table with a column of admin codes
#' @param admin_level The administrative level of the codes in the
#'   \code{df} object.  Integer; 0, 1, or 2.
#' @param idx_col which column in \code{df} contains the administrative
#'   codes?  A character.
#' @param sp_h Spatial hierarchy object generated with
#'   \code{get_sp_hierarchy()}. If \code{NULL}, then will be generated
#'   within the function. If using this function several times, it may
#'   save time to generate the sp_h object once, then pass it to this
#'   function
#' @param shapefile_version String specifying shapefile version to pull of is.null(sp_h)
#' @return Data table with sp_hierarchy added on.  Note that the
#'   original admin_code column will be renamed to conform with
#'   the "ADMX_CODE" convention.
#' @examples
#' # add spatial hierarchy information to df, which contains
#' # admin 2 codes in the column "spatial_idx", using a pre-made
#' # spatial hierarchy object stored in sp_h
#' df2 <- merge_sp_hierarchy(df = df,
#'                           admin_level = 2,
#'                           idx_col = "spatial_idx",
#'                           sp_h = sp_h)

merge_sp_hierarchy <- function(df, admin_level, idx_col, sp_h = NULL, shapefile_version = 'current') {

  # Grab the spatial hierarchy
  if (is.null(sp_h)) sp_h <- get_sp_hierarchy(shapefile_version = shapefile_version)

  df <- as.data.table(df) %>%
        setnames(., idx_col, paste0("ADM", admin_level, "_CODE"))

  # grab df names but not the spatial idx
  df_names <- names(df)[names(df) != paste0("ADM", admin_level, "_CODE")]

  sp_idx_table <- sp_h[[paste0("ADM", admin_level)]]

  df <- merge(sp_idx_table, df,
               by = c(paste0("ADM", admin_level, "_CODE")),
               all.y = T, all.x = F)

  idx_names <- sort(names(sp_idx_table))

  setcolorder(df, c(idx_names, df_names))
  setorderv(df, c(paste0("ADM", 0:admin_level, "_CODE")))

  return(df)
}

