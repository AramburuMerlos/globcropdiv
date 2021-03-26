library(data.table)
if(!getwd() %like%  "globcropdiv$"){
  thispc <- system('hostname', TRUE)
  if (thispc == "LAPTOP-IVSPBGCA") { 
    setwd("D:/gdrive/share/globcropdiv/") #I'm guessing
  } else if (thispc == "ESP-RH-9891" | thispc == "LAPTOP-ST129J47") { 
    setwd("G:/My Drive/globcropdiv/")
  } else {
    # setwd
  }
}
