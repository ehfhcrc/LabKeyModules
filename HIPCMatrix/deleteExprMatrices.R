# vim: sw=4:ts=4:nu:nospell:fdc=4
#
#  Copyright 2013 Fred Hutchinson Cancer Research Center
#
#  Licensed under the Apache License, Version 2.0 (the 'License');
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an 'AS IS' BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# dependencies
library(data.table)
labkey.url.base <- "https://test.immunespace.org"
labkey.url.path <- "/Studies/"
schemaName = "assay.ExpressionMatrix.matrix"
queryName = "Runs"


# get run names and print out
runs <- data.table( labkey.selectRows( baseUrl = labkey.url.base,
                                       folderPath = labkey.url.path,
                                       schemaName = schemaName,
                                       queryName = queryName,
                                       showHidden = T))
print( runs[ , c("Name","Study", "Row ID") ] )


# Ask user which runs to remove and verify
runs_rm <- readline( prompt = "Please enter the `Row Id` for runs to remove separated by commas: ")

runs_rm <- as.integer( unlist( strsplit( runs_rm, split = "," ) ) )
chk <- all( runs_rm %in% runs$`Row Id` )
if(chk == TRUE){
  chkd_runs <- runs[ `Row Id` %in% runs_rm, ]
  cat( "Runs to be removed: " )
  print( chkd_runs  )
}else{
  bad_ids <- runs_rm[ !( runs_rm %in% runs$`Row Id` ) ]
  stop(paste0("Following Runs not in available row ids: ", bad_ids, ". Please Retry"))
}

sure <- readline( prompt = "Are you sure? Y/n")
sure <- ifelse( sure == "" | sure == "Y", TRUE, FALSE)

# Remove runs from the DB via deleteRows, but need to change path for each study
if(sure == TRUE){
  rm_paths <- list()
  for(i in 1:dim(chkd_runs)[1]){
    row_rm <- data.frame(RowId = chkd_runs$`Row Id`[[i]], 
                         stringsAsFactors = FALSE
                         )
    deletedRows <- labkey.deleteRows(
      toDelete = row_rm,
      queryName = queryName,
      baseUrl       = labkey.url.base,
      folderPath    = paste0( labkey.url.path, "/", chkd_runs$Study[[i]] ),
      schemaName    = schemaName)
    
    rm_paths[[ chkd_runs$Name[[i]] ]] <- deletedRows$rows[[1]]$filepathroot
  }
}else{
  stop("Rows not deleted.  Please retry when ready.")
}

cat("file paths to be deleted")
print(rm_paths)

rm_files <- readline( prompt = "Remove files? Y/n")
rm_files <- ifelse( rm_files == "" | rm_files == "Y", TRUE, FALSE)
if(rm_files == TRUE){
  # remove flat files from server
  results <- lapply(rm_paths, FUN = function(path){
    suppressMessages(tryCatch(
      file.remove(path),
      error = function(e) return(e),
      warning = function(w) return(w)
    ))
  })
  
  print(results)
}else{
  stop("Operations ended. File removal aborted.")
}




