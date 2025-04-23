library(openxlsx)
pitt_meta=read.csv("/ix/djishnu/Marisa/ABomics/FINAL_MODELS/pitt_meta.csv")
HLAgen= read.xlsx("./input_data/PittCohort_kidneytransplant_ABMR_HLAgenotype.xlsx", 
                  sheet= "HLA rec and donor")

HLAgen= HLAgen[HLAgen$CUSTOMER_ID %in% pitt_meta$CUSTOMER_ID,]

donor.df= HLAgen[,c("CUSTOMER_ID",grep("^DONOR",colnames(HLAgen), value=TRUE))]
colnames(donor.df)= gsub("DONOR_","", colnames(donor.df))
donor.df[-1]=sapply(donor.df[-1], function(x) as.numeric(x))
recip.df= HLAgen[,c("CUSTOMER_ID",grep("^RECIP",colnames(HLAgen), value=TRUE))]
colnames(recip.df)= gsub("RECIPIENT_","", colnames(recip.df))
recip.df[-1]=sapply(recip.df[-1], function(x) as.numeric(x))
# Function to add leading zero to single-digit numbers
add_leading_zero <- function(x) {
  if (is.na(x) || x=="-") {
    return("-")  # Return as is for NA or "-"
  } else if (x < 10) {
    return(sprintf("%02d", as.numeric(x)))  # Format with leading zero
  } else {
    return(as.character(x))  # Return as character for other numbers
  }
}

add_serotype = function(df,col_index){
  
  col_name= colnames(df)[col_index]

  df[[col_index]]= sapply(df[[col_index]], function(x){
    if (is.na(x) || x =="-"){
      return("-")
    } else {
      return(paste0(col_name,x))
    }
    })
  
  return(df)
}

recip.df[,-1]= lapply(recip.df[,-1],(function(col) sapply(col,add_leading_zero)))
donor.df[,-1]= lapply(donor.df[,-1],(function(col) sapply(col,add_leading_zero)))

# Find mismatches as a logical matrix
mismatches <- mapply(function(row1, row2) row1 != row2, recip.df[-1], donor.df[-1])
mismatch_count <- rowSums(mismatches)
mismatch_count


# Find mismatches as a logical matrix
mismatches2 <- mapply(function(col1, col2) col1 != col2, recip.df[-1], donor.df[-1])

# Sum mismatches by column
mismatch_count_by_column <- colSums(mismatches2)

# Identify columns with the highest number of mismatches
max_mismatches <- mismatch_count_by_column[mismatch_count_by_column == max(mismatch_count_by_column)]

# Display results
mismatch_count_by_column
max_mismatches