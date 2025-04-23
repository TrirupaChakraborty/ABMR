
sdrf_txt <- read.table("/ix/djishnu/Trirupa/ABomics.Prj/SigmiR/E-MTAB-12051.sdrf.txt", 
                       header = TRUE, 
                       sep = "\t", 
                       fill = TRUE)
#colnames(sdrf_txt)
#unique(sdrf_txt$Characteristics.disease.)

sdrf_filt= sdrf_txt[sdrf_txt$Characteristics.disease %in% c("Non rejection DSA+","Antibody-mediated rejection"),]

file_types=c("I1","R1","R2")
dir_path="/ix/djishnu/Trirupa/ABomics.Prj/SigmiR/scRNAseq_EMBL/"
for (url in sdrf_filt$Comment.FASTQ_URI.){
  url_parts=strsplit(url,"_")[[1]]
  file_rootname= paste0(url_parts[1:4], collapse = "_")
  print(file_rootname)
  
  for (ftype in file_types){
    download_url= paste0(file_rootname,"_",ftype,"_",url_parts[6])
    #print(download_url)
    file_name=strsplit(download_url,"/")[[1]][8]
    destination= paste0(dir_path,file_name)
    download.file(url, destination, method = "wget")
  }
}
  

