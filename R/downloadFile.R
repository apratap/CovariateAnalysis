# Utility function to download tsv or csv file from synapse and load it in to memory
downloadFile <- function(id, ...){
  tmp = data.table::fread(synapseClient::synGet(id)@filePath, data.table=F, header=T, ...)
}
