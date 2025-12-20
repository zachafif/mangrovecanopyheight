## Set the folder containing the CSV files
path <- "H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\Trat_BioHeight_30112025-20251209T123704Z-1-001\\Trat_BioHeight_30112025\\DataRevised_30112025"   # change this

## List all CSV files (full paths)
files <- list.files(path = path,
                    pattern = "\\.csv",
                    full.names = TRUE)

# read all as character first (optional but often safer)
df_list <- lapply(
  files,
  \(f) read_csv(f, col_types = cols(.default = col_character()))
)

# now bind all rows
combined_df <- bind_rows(df_list)

combined_df= type.convert(combined_df, as.is = TRUE)
combined_df$Lon=as.numeric(combined_df$Lon)

combined_df= combined_df %>%
  mutate(`Height(m)` = coalesce(`Height(m)`, `Height (m)`), Lon=coalesce(Lon,Long))

summ_df=combined_df %>%
  group_by(Plotnumber) %>%
  summarise(
    Lat = first(na.omit(Lat)),
    Lon = first(na.omit(Lon)),
    n   = n(),                            # count rows in each group
    mean_x = mean(`Height(m)`, na.rm = TRUE)        # example summary
  )%>% drop_na()

point_gcp=sf_multipoint(
  summ_df,
  x = "Lon",
  y = "Lat",
  multipoint_id = "Plotnumber",
  keep = TRUE
)

sf::st_crs(point_gcp)=4326

write_sf(point_gcp,"H:\\02 MSc\\MASTER 2025\\JGSEE\\Paper\\Trat_BioHeight_30112025-20251209T123704Z-1-001\\Trat_BioHeight_30112025\\DataRevised_30112025\\GCP.shp")