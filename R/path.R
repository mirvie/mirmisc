#' Get the path to the folder containing the feather files.
#'
#' This function requires you to have set the environment variable
#' `MIRVIE_FEATHER_PATH`, which you can do in the `~/.Rprofile` file. It should
#' have a line like `Sys.setenv(MIRVIE_FEATHER_PATH =
#' "path/to/mirvie/feathers/dir")`. If the file doesn't exist, create it and
#' make this the only line in the file. If this is done correctly, this function
#' then forms a path with `MIRVIE_FEATHER_PATH` as the root directory.
#'
#' There's a whole vignette explaining this function. To find it, run
#' `browseVignettes(package = "mirmisc")`.
#'
#' @param ... Character vectors. Elements of the path. Mostly, you'll leave this
#'   blank.
#' @param use_dotenv A flag. If the `MIRVIE_FEATHER_PATH` environment variable
#'   isn't found by the usual R means, check the `~/.env` file used by
#'   `python-dotenv`.
#' @param verify A flag. Check that `MIRVIE_FEATHER_PATH` exists and contains at
#'   least one *.feather file. Error if check fails. Default `TRUE`.
#'
#' @return An [fs::path].
#'
#' @examples
#' \dontrun{
#' get_feather_path()
#' }
#'
#' @export
get_feather_path <- function(..., use_dotenv = TRUE, verify = TRUE) {
  checkmate::assert_flag(use_dotenv)
  checkmate::assert_flag(verify)
  feather_dir <- Sys.getenv("MIRVIE_FEATHER_PATH")
  if (isTRUE(feather_dir == "") && use_dotenv) {
    if (file.exists("~/.env")) {
      readRenviron("~/.env")
      feather_dir <- Sys.getenv("MIRVIE_FEATHER_PATH")
    }
  }
  # It's very important for this function to have informative error messages.
  # Most of the code inside of it is dedicated to erroring.
  err_msg_no_envvar_found <- paste(
    "No environment variable {code('MIRVIE_FEATHER_PATH')} found."
  ) # this and following strings will get `stringr::str_glue()` later
  err_msg_edit_your_rprofile <- paste(
    "Use {code('usethis::edit_r_profile()')} to edit your",
    "{path('~/.Rprofile')} file."
  )
  err_msg_add_setenv_line <- paste(
    "Add a line",
    "{code('Sys.setenv(MIRVIE_FEATHER_PATH = \"/path/to/dir\")')}",
    "to locate your mirvie repo."
  )
  err_msg_restart_r <- paste(
    "Restart R for these changes to take effect."
  )
  err_msg_file_not_dir <- paste(
    "The path {path(feather_dir)} is to a file,",
    "not a directory."
  )
  err_msg_make_sure <- paste(
    "Make sure it's setting the environment variable",
    "{code('MIRVIE_FEATHER_PATH')}",
    "to the location of your feather files."
  )
  err_msg_dir_not_exist <- paste(
    "The directory {path(feather_dir)} does not exist."
  )
  err_msg_no_feathers <- paste(
    "No *.feather files found in the directory {path(feather_dir)} pointed to",
    "by the MIRVIE_FEATHER_PATH environment variable."
  )
  if (identical(feather_dir, "")) {
    custom_stop(
      err_msg_no_envvar_found,
      err_msg_edit_your_rprofile,
      err_msg_add_setenv_line,
      err_msg_restart_r
    )
  }
  if (verify) {
    if (!fs::dir_exists(feather_dir)) {
      if (fs::file_exists(feather_dir)) {
        custom_stop(
          err_msg_file_not_dir,
          err_msg_edit_your_rprofile,
          err_msg_make_sure,
          err_msg_restart_r
        )
      }
      custom_stop(
        err_msg_dir_not_exist,
        err_msg_edit_your_rprofile,
        err_msg_make_sure,
        err_msg_restart_r
      )
    }
    feathers_present <- fs::dir_ls(feather_dir) %>%
      stringr::str_ends(stringr::fixed(".feather")) %>%
      any()
    if (!feathers_present) {
      custom_stop(
        err_msg_no_feathers,
        err_msg_edit_your_rprofile,
        err_msg_make_sure
      )
    }
  }
  fs::path(feather_dir, ...)
}
