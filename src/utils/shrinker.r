setCurrentDirectory <- function() {
  setwd(getSrcDirectory(function(){})[1])
}


shrinkData <- function(input, output, delete_percentage = NULL, amount = 500 ) {
  print("shrinkData start...")
  
  # Read the CSV file
  data <- read.csv(input)
  
  print(paste0("Rows: ", nrow(data), " | Cols: ", ncol(data)))
  
  
  # Calculate the number of rows to delete (excluding the header)
  total_rows <- nrow(data)
  
  if(!is.null(delete_percentage)){
    rows_to_delete <- floor((delete_percentage / 100) * (total_rows - 1))
  }
  else{
    rows_to_delete <- (total_rows - amount)
  }

  # Randomly select rows to delete (excluding the first row)
  rows_to_keep <- c(1, sample(2:total_rows, total_rows - 1 - rows_to_delete))
  
  # Subset the data to keep only the selected rows
  data_shrunken <- data[rows_to_keep, ]
  
  # Write the shrunk data to a new CSV file
  write.csv(data_shrunken, output, row.names = FALSE)
}
setCurrentDirectory()

#datasetName <- "CRCCombinedA1A3InvAd"
datasetName <- "5KPBMC"

datasetPath <- paste0("../../datasets/csv/",datasetName,".csv") 
outputPath <- paste0("../../datasets/csv/",datasetName,"_small.csv") 

shrinkData(datasetPath, outputPath, amount = 500)

print("Done.")
