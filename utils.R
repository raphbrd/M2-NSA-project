load_gene_data <- function(filename, output_name) {
  df <- read.csv(filename, sep = ";", header = FALSE, skip = 3)
  df <- t(df)
  df <- df[3:dim(df)[1],]
  
  genes_express_idxes <- seq(1, dim(df)[1], by = 2)
  df <- as.data.frame(df[genes_express_idxes,])
  df <- df %>% mutate(across(everything(), as.integer))
  dim(df)
  
  cls <- read.delim(str_replace(filename, ".csv", ".cls"), comment.char = "#", header = FALSE)
  header <- stringr::str_trim(cls[1, ])
  cls <- stringr::str_trim(cls[2,])
  cls <- as.numeric(stringr::str_split(cls, " ")[[1]])
  
  # vérification que le nombre de gènes est cohérent avec le header
  n <- as.numeric(stringr::str_split(header, " ")[[1]][1])
  
  if (n != dim(df)[1]) {
    stop("Inconsistent number of observations between the cls and the csv files")
  }
  
  write.csv(df, file=output_name)
  
  return(list(data=df, classes=cls))
}
