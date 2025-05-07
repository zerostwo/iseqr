#' @importFrom stringr str_split str_subset
#' @importFrom purrr discard
#' @importFrom dplyr filter mutate select group_by summarise case_when distinct all_of across
#' @importFrom stringr str_trim str_detect str_remove
#' @importFrom magrittr %>%
#' @importFrom utils URLencode tail
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", ".data", "filename", "md5"))
}