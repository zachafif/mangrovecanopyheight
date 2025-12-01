parse_d12_line <- function(line) {
  # Remove trailing whitespace if any
  line <- trimws(line)
  
  # Numeric string is from start to position before "lg"
  numeric_str <- substr(line, 1, 72)
  numeric_str <- gsub(" ", "", numeric_str)  # remove spaces
  
  # Split into groups of 3 digits (24 hours)
  hourly_vals <- substring(numeric_str, seq(1, nchar(numeric_str), 3), seq(3, nchar(numeric_str), 3))
  hourly_vals <- as.numeric(hourly_vals)
  hourly_vals[hourly_vals == 999] <- NA  # missing values
  
  if(length(hourly_vals) != 24) {
    warning(paste("Hourly value count is not 24, found:", length(hourly_vals)))
  }
  
  df <- data.frame(
    Hour = sprintf("%02d:00", 0:(length(hourly_vals) - 1)),
    WaterLevel_cm = hourly_vals
  )
  return(df)
}

parse_d12_file <- function(filepath) {
  lines <- readLines(filepath)
  parsed_list <- lapply(lines, parse_d12_line)
  parsed_list <- parsed_list[!sapply(parsed_list, is.null)]
  do.call(rbind, parsed_list)
}


tide.18 <- parse_d12_file("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\LG2018.D12")
tide.18$X=1:nrow(tide.18)

year <- 2018
hour_of_year <- tide.18$X  # example hour of year (0 = Jan 1, 00:00)

start_date <- as.POSIXct(paste0(year, "-01-01 00:00:00"), tz = "UTC")

datetime <- start_date + (hour_of_year * 3600)  # 3600 seconds = 1 hour
tide.18$datetime=datetime

tide=rbind(tide.18,tide.19,tide.20)


write.table(tide,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\tidetrat.csv",sep=";",row.names = FALSE,quote = FALSE)


t1=min(tide$datetime)
t2=max(tide$datetime)

plot(ftide(tide$WaterLevel_cm,tide$datetime),t1,t2)
title("LAEMNGOP Tidal Station Observation 2018-2020 ")

dat=read.delim("H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\dataset_v10_fin.csv",sep=";")
dat$datetime=strptime(dat$time, format = "%Y-%m-%d %H:%M:%S")
dat$datehr=format(dat$datetime, "%Y-%m-%d %H")

tide$datehr=format(tide$datetime, "%Y-%m-%d %H")
tide$lotide[tide$WaterLevel_cm <=167]=1

join=left_join(dat,tide,join_by(datehr))

tide.filt=join[join$lotide==1,]

write.table(join,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\datjointide_167.csv",sep=";",row.names = FALSE,quote = FALSE)


