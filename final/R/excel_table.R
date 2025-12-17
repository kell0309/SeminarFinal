#' Excel plot script
#'
#' @description
#' Writes a differential expression results table to an excel file using `openxlsx`.
#' The output contains one sheet named "DEG"
#'
#' @param deg_table A data.frame  containing DEG results to export.
#' @param out_xlsx  Output Excel filename or path
#'
#' @return  The path/filename of the saved Excel file.
#' @export

export_deg_excel <- function(deg_table, out_xlsx = "DEG_results.xlsx") {

  # create a new excel workbook
  wb <- createWorkbook()

  # Write DEG results table into a worksheet called "DEG"
  addWorksheet(wb, "DEG")
  writeData(wb, "DEG", deg_table)
  freezePane(wb, "DEG", firstRow = TRUE)

  # Saves the workbook to disk
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)

  # return the output file path so that the caller knows where it was saved
  out_xlsx
}
