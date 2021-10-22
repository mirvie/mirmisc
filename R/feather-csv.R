#' Make a directory of CSVs from a directory of feathers.
#'
#' Take a directory containing feather files and create a sibling directory with
#' corresponding CSVs.
#'
#' @param feather_dir_path The path to the directory containing the feathers.
#' @param new_dir_name The name of the new directory. This should _not_ be an
#'   absolute path and rather just a name; I.e. it should not contain '/'. The
#'   created directory will be a sibling of the one at `feather_dir_path`.
#'
#' @return The path to the output directory, invisibly.
#'
#' @export
convert_feather_dir_to_csvs <- function(feather_dir_path,
                                        new_dir_name = "feather-csvs") {
  checkmate::assert_string(feather_dir_path, min.chars = 1)
  checkmate::assert_directory_exists(feather_dir_path)
  checkmate::assert_string(new_dir_name, min.chars = 1)
  new_dir_name <- strex::str_trim_anything(new_dir_name, stringr::coll("/"))
  checkmate::assert_false(stringr::str_detect(new_dir_name, stringr::coll("/")))
  feather_dir_dir <- fs::path_dir(feather_dir_path)
  feather_csv_path <- fs::path(feather_dir_dir, new_dir_name)
  if (!fs::dir_exists(feather_csv_path)) fs::dir_create(feather_csv_path)
  feather_file_paths <- fs::dir_ls(feather_dir_path, glob = "*.feather")
  for (ffp in feather_file_paths) {
    df <- arrow::read_feather(ffp)
    csv_file_path <- fs::path(
      feather_csv_path,
      strex::str_give_ext(fs::path_file(ffp), "csv", replace = TRUE)
    )
    readr::write_csv(df, csv_file_path)
  }
  invisible(feather_csv_path)
}
