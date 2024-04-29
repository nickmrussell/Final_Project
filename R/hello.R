
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

rm(list=ls())

library(easypackages)

libraries(c("tidyverse",
            "tools",
            "ggrepel",
            "plyr",
            "rstudioapi",
            "pdftools"))

###
# choose folder where CSV files are located, make sure to only have the data converted CSV files (from ABF format) in this folder
#setwd(rstudioapi::selectDirectory())

# manually set which file you want to run on below
# uncomment the one you want or select one
showDialog(title = "Select Data file",
           message = "Select CSV Data File to run this Script on")
mycsvfile <- selectFile(caption = "Select CSV Data File to run this Script on",
                        filter = "CSV Files (*.csv)",
                        existing = TRUE)

# Select Experiment Identifiers file that should be a CSV file in directory above which contains four columns:
# ObsID, MouseID, CellID, Filename_IV
# for later merging with data output where Filename_IV contains names of your CSV files
# without the .csv extension added to the name
# This will obviously be specific for each experiment so edit this accordingly
# to match the length of the number of CSV files you are trying to process in this directory
showDialog(title = "Select Experiment Identifiers file",
                       message = "Select Experiment Identifiers file that corresponds to this CSV Data File")
identifiers <- read_csv(selectFile(caption = "Select Experiment Identifiers file that corresponds to this CSV Data File",
                                   filter = "CSV Files (*.csv)",
                                   existing = TRUE))

showDialog(title = "Select Output Folder",
           message = "Select Folder where you would like Output data to be saved.")
setwd(selectDirectory())


# get filename without csv extension
filename <- file_path_sans_ext(basename(mycsvfile))

dir.create(paste("Single Data Output ",
                 get("filename"),
                 sep = ""))

setwd(paste("Single Data Output ",
            get("filename"),
            sep = ""))

dir.create("Plots")

all_data <- NULL
finaldf <- NULL

all_APs_all_sweeps <- NULL
all_files_APs_all_sweeps <- NULL

# import data and get rid of top two rows since the way the data is output from abf to csv
# these are not necessary as we will transform and fix column names later
# and then this removes any empty columns that are only NAs
df <- read_csv(mycsvfile, col_names = FALSE, skip = 2) %>%
  select(where(~!all(is.na(.x))))

# count total number of sweeps and create vectors for new names of columns based on number of sweep
# assumes first column 1 is Time data and remaining columns 2 and onward are Sweeps data
sweeps_count <- df %>% select(2:last_col()) %>% ncol()
old_names <- df %>% select(2:last_col()) %>% names() %>% as.vector()
new_names <- as.character(seq(1:sweeps_count)) %>% as.vector()

# rename the columns based on Sweep number
# make this a column called Sweep instead and the value the number of the Sweep
# Group by Sweep so operations on run on each Sweep and not total data
# Transmute renames columns, converts Time from seconds to milliseconds,
# Fix column class type, and creates a new dV/dT column
# Create new column of the voltage values normalized for counting peaks later
# Arrange by Sweep so data is sorted Sweep 1 to Sweep 2 and so on
# Ungroup to add new total time elapsed column
# Regroup data back by Sweep for future operations
df1 <- df %>%
  rename_with(~ new_names[which(old_names == .x)], .cols = old_names) %>%
  pivot_longer(cols = 2:last_col(),
               values_to = "Voltage (mV)",
               names_to = "Sweep") %>%
  group_by(Sweep) %>%
  transmute(`Time (ms)` = as.numeric(X1)*1000,
            `Voltage (mV)` = as.numeric(`Voltage (mV)`),
            `dV/dT` = (`Voltage (mV)` - lag(`Voltage (mV)` )) / (`Time (ms)` - lag(`Time (ms)`)),
            Sweep = as.numeric(Sweep),
            `Normalised Voltage (mV)` = as.numeric(`Voltage (mV)` + abs(min(`Voltage (mV)`)))) %>%
  arrange(Sweep)

# remove variables we don't need later
rm(df, sweeps_count, old_names, new_names)

# create NA variables which are overwritten later if they exist
# If these NAs are not created then later parts will fail to run to
# completion and create the data output file
this_identifier <- NULL
this_identifier$ObsID <- NA
this_identifier$MouseID <- NA
this_identifier$CellID <- NA
this_identifier$Filename_IV <- NA

`Sweep 9 Resting Em (mV)` <- NA
`Sweep 9 -50 pA step Steady-State Em` <- NA
`Sweep 9 Resistance_in (MegaOhms)` <- NA
`Sweep # for 1st AP fired` <- NA
`Rheobase (pA)` <- NA
`Threshold Em (mV)` <- NA
`Peak Voltage (mV)` <- NA
`Halfwidth (ms)` <- NA
`AP Amplitude (mV)` <- NA
`Half Amplitude (mV)` <- NA
`Halfwidth (ms)` <- NA

`maxAHP_1stAP Voltage (mV)` <- NA
`maxAHP_1stAP Time (ms)` <- NA
`maxAHP_1stAP Time Relative to Threshold` <- NA
maxAHPdeltathreshold <- NA
`fAHP_1stAP Voltage (mV)` <- NA
`fAHP_1stAP Time (ms)` <- NA
`fAHP_1stAP Time Relative to Threshold` <- NA
fAHP_delta_threshold <- NA
mAHP_1stAP <- NA
mAHP_delta_threshold <- NA
sAHP_1stAP <- NA
sAHP_delta_threshold <- NA

`Sweep # for SFA` <- NA
`Sweep # for Peak to peak Frequency_Max` <- NA
`Peak to peak Frequency_Max` <- NA
`Interevent Interval_1` <- NA
`Interevent Interval_last` <- NA
`SFA (1st ISI/last ISI)` <- NA


`Resting Em` <- df1 %>%
  filter(between(`Time (ms)`, 0, 100) & Sweep == 9) %>%
  summarise(mean_0_100 = mean(`Voltage (mV)`)) %>%
  select(mean_0_100) %>%
  unlist() %>%
  as.numeric()

`-50 pA step Steady-State Em` <- df1 %>%
  filter(between(`Time (ms)`, 198.3, 215) & Sweep == 9) %>%
  summarise(mean_198_215 = mean(`Voltage (mV)`)) %>%
  select(mean_198_215) %>%
  unlist() %>%
  as.numeric()

# this is using Ohm's law of R = V/I
# where voltz is in millivolts and current (50) is in nanoamps
# therefore result is in megaOhms after multiplying by 1000
`Rin (MOhm)` <- ((`Resting Em` - `-50 pA step Steady-State Em`)/50)*1000

# get the AP data for all the sweeps
# this compares each peak to the points around it
# (specifically points 2.5 ms, 5 ms, and 10 ms before and after)
# and if any of these pass a certain height (15, 20, or 20, respectively),
# then they are considered an action potential.

sweeps1 <- df1 %>%
  arrange(Sweep) %>%
  filter(((`Normalised Voltage (mV)` - lag(`Normalised Voltage (mV)`, 25)) >= 15 &
            (`Normalised Voltage (mV)` - lead(`Normalised Voltage (mV)`, 25)) >= 15) |
           ((`Normalised Voltage (mV)` - lag(`Normalised Voltage (mV)`, 50)) >= 20 &
              (`Normalised Voltage (mV)` - lead(`Normalised Voltage (mV)`, 50)) >= 20) |
           ((`Normalised Voltage (mV)` - lag(`Normalised Voltage (mV)`, 100)) >= 20 &
              (`Normalised Voltage (mV)` - lead(`Normalised Voltage (mV)`, 100)) >= 20)) %>%
  mutate(Filename_IV = mycsvfile)

sweeps2 <- df1 %>%
  filter(lead(`dV/dT`) <= 0 &
           `dV/dT` > 0) %>%
  mutate(Filename_IV = mycsvfile)

all_APs_all_sweeps <- inner_join(sweeps1, sweeps2)

# filter only unique split peaks that occur at least 3 ms away from each other so we are not iterating through redundant peaks with our 3 ms window later
split_peaks <- all_APs_all_sweeps %>%
  group_by(Sweep) %>%
  filter(abs(`Time (ms)` - lag(`Time (ms)`)) <= 3 |
         abs(`Time (ms)` - lead(`Time (ms)`)) <= 3)

splitdf <- NULL
truepeaks <- NULL

# get rid of split peaks where there are two points very close together (within plus or minus 3 ms)
if (nrow(split_peaks) >= 1) {
  for (thisrow in 1:nrow(split_peaks)) {
    tmpdf = split_peaks[thisrow,]
    tmptime = tmpdf$`Time (ms)`
    tmpsplitdf = split_peaks %>%
      filter(Sweep == tmpdf$Sweep &
               between(`Time (ms)`, (tmptime - 3), (tmptime + 3)))
    maxsplit = tmpsplitdf %>% slice_max(n = 1, order_by = `Voltage (mV)`) %>% head(1)
    splitdf = rbind(splitdf, tmpsplitdf)
    truepeaks = rbind(truepeaks, maxsplit)
  }

  truepeaks <- truepeaks %>% distinct(.keep_all = TRUE)
  splitdf <- splitdf %>% distinct(.keep_all = TRUE)

  no_split_peaks <- anti_join(all_APs_all_sweeps, splitdf)

  all_APs_all_sweeps <- rbind(truepeaks, no_split_peaks) %>%
    group_by(Sweep) %>%
    arrange(Sweep, `Time (ms)`)
}

sweeps_AP_count <- all_APs_all_sweeps %>%
  select(Sweep) %>%
  table() %>%
  as.data.frame()

# match up filename to other identifiers in identifier data like ObsID, Mouse, and Cell
this_identifier <- identifiers %>% filter(Filename_IV == filename)

# remove variables we don't need later
rm(excessive_counts,
   firstfilter,
   secondfilter,
   all_APs_all_sweeps2,
   group1,
   group2)

########################################################################
########################################################################
############ IF LOOP BEGINS HERE FOR THOSE SWEEPS WITH ACTION POTENTIALS
########################################################################
########################################################################
########################################################################

if (nrow(sweeps_AP_count) > 1) {

  sweeps_AP_count <- sweeps_AP_count %>%
    transmute("Sweep" = Sweep,
              `AP Count` = Freq)

  # fix variable type for Sweep else error will occur later on
  sweeps_AP_count$Sweep <- as.character(sweeps_AP_count$Sweep)

  # if there are no Action Potentials in any Sweeps then all of these values below should be NA
  # get sweep with first AP
  `Sweep # for 1st AP fired` <- sweeps_AP_count %>%
    slice(1) %>%
    select(Sweep) %>%
    unlist()  %>%
    as.numeric()

  # create data frame to count number of APs for all sweeps 10 and onward including those with no APs
  # the subtracted 1 is so it doesn't include the first Sweep with an AP otherwise will be wrong
  sweeps_before_first_AP_count <- tibble(Sweep = c(10:(`Sweep # for 1st AP fired` - 1)),
                                         `AP Count` = 0)

  AP_count_all_sweeps <- rbind(sweeps_before_first_AP_count, sweeps_AP_count)

  firstAPsweep <- df1 %>%
    filter(Sweep == `Sweep # for 1st AP fired`)

  firstAPsweep_firstpeak <- all_APs_all_sweeps %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired`) %>%
    slice_head(n = 1)

  # get values of data from 15 ms before first peak and 20 ms after first peak
  Time1 <- (firstAPsweep_firstpeak %>% select(`Time (ms)`) %>% unlist() %>% as.numeric()) - 15
  Time2 <- (firstAPsweep_firstpeak %>% select(`Time (ms)`) %>% unlist() %>% as.numeric()) + 20

  firstAP_range <- firstAPsweep %>%
    filter(between(`Time (ms)`,Time1, Time2))


  `Threshold Em (mV)` <- firstAP_range %>%
    filter(lead(`dV/dT` > 10)) %>%
    slice(1) %>%
    ungroup() %>%
    select(`Voltage (mV)`) %>%
    unlist() %>%
    as.numeric()


  `Threshold Time (ms)` <- firstAP_range %>%
    filter(lead(`dV/dT` > 10)) %>%
    slice(1) %>%
    ungroup() %>%
    select(`Time (ms)`) %>%
    unlist() %>%
    as.numeric()


  `Peak Voltage (mV)` <- max(firstAP_range$`Voltage (mV)`)
  `AP Amplitude (mV)` <- `Peak Voltage (mV)` - `Threshold Em (mV)`
  `Half Amplitude (mV)` <- `AP Amplitude (mV)`/2

  # get lowest number of voltage (as in negative not small)

  `maxAHP_1stAP Voltage (mV)` <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             between(`Time (ms)`,
                     `Threshold Time (ms)`,
                     (`Threshold Time (ms)` + 15))) %>%
    slice_min(n = 1, order_by = `Voltage (mV)`) %>%
    select(`Voltage (mV)`) %>%
    unlist() %>%
    as.numeric() %>%
    head(1)

  # find the Time that corresponds to the `maxAHP_1stAP Voltage`
  `maxAHP_1stAP Time (ms)` <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             between(`Time (ms)`,
                     `Threshold Time (ms)`,
                     (`Threshold Time (ms)` + 15))) %>%
    slice_min(n = 1, order_by = `Voltage (mV)`) %>%
    select(`Time (ms)`) %>%
    unlist() %>%
    as.numeric() %>%
    head(1)

  `maxAHP_1stAP Time Relative to Threshold` <- abs(`Threshold Time (ms)` - `maxAHP_1stAP Time (ms)`)


  # calculate halfwidth
  `HalfWidth Volts` <- `Threshold Em (mV)` - `Half Amplitude (mV)`
  
  # dV = V2 - V1
  # so V2 = dV + V1
  V2 <- `Threshold Em (mV)` + `Half Amplitude (mV)`

  # filter times greater than Threshold Time and then find which two voltage points that V2 lies between
  # slice for just the first two in case there are more than one
  half_ap_range <- firstAP_range %>%
    filter(`Time (ms)` >= `Threshold Time (ms)` &
             V2 < lag(`Voltage (mV)`) &
             V2 > lead(`Voltage (mV)`)) %>%
    slice(1:2)

  # this finds formula of curve around point
  xs <- half_ap_range %>%
    ungroup() %>%
    pull(`Time (ms)`)

  ys <- half_ap_range %>%
    ungroup() %>%
    slice(1:2) %>%
    pull(`Voltage (mV)`)

  linearfit <- lm(ys ~ xs)

  slope <- linearfit$coefficients["xs"] %>% as.numeric()
  y_intercept <- linearfit$coefficients["(Intercept)"] %>% as.numeric()

  # Rearrange and use curve to solve for x which is T2
  T2 <- (V2 - y_intercept)/slope

  # find Halfwidth by subtracting Threshold Time from T2
  `Halfwidth (ms)` <- (T2 - `Threshold Time (ms)`)

  # remove variables we don't need later
  rm(half_ap_range,
     slope,
     y_intercept,
     xs,
     ys,
     V2,
     `HalfWidth Volts`,
     T2,
     Time1,
     Time2,
     linearfit)

  # create a dataframe of what rheobase value corresponds to what sweep
  # here we just said sweeps 10 through 20 with their corresponding values
  # from 50 to 550 with increasing 50 steps in between
  rheobase_df <- tibble(Sweep = c(10:20), Rheobase = c(seq(50, 550, 50)))

  `Rheobase (pA)` <- rheobase_df %>%
    filter(Sweep == `Sweep # for 1st AP fired`) %>%
    pull(Rheobase)


  `fAHP_1stAP Voltage (mV)` <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             between(`Time (ms)`,
                     `Threshold Time (ms)`,
                     (`Threshold Time (ms)` + 15))) %>%
    slice_min(n = 1, order_by = `Voltage (mV)`) %>%
    slice_head(n = 1) %>% 
    pull(`Voltage (mV)`)

  `fAHP_1stAP Time (ms)` <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             between(`Time (ms)`,
                     `Threshold Time (ms)`,
                     (`Threshold Time (ms)` + 5))) %>%
    slice_min(n = 1, order_by = `Voltage (mV)`) %>%
    slice_head(n = 1) %>% 
    pull(`Time (ms)`)

  `fAHP_1stAP Time Relative to Threshold` <- abs(`Threshold Time (ms)` - `fAHP_1stAP Time (ms)`)


  mAHP_1stAP <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             `Time (ms)` == as.character(`Threshold Time (ms)` + 10)) %>%
    pull(`Voltage (mV)`)

  sAHP_1stAP <- df1 %>%
    ungroup() %>%
    filter(Sweep == `Sweep # for 1st AP fired` &
             `Time (ms)` == as.character(`Threshold Time (ms)` + 15)) %>%
    pull(`Voltage (mV)`)

  `Sweep # for SFA` <- AP_count_all_sweeps %>%
    arrange(desc(`AP Count`)) %>%
    pull(Sweep) %>% 
    head(1)


  # this gets the Interevent Interval by subtracting current time from previous time
  # this gets Peak-to-Peak Frequency by taking inverse of the difference of the
  # current Peak time with the last Peak time which gives you mHZ which then needs to be
  # multiplied by 1000 to get Hz
  # so this is simplified by just dividing 1000 by time difference to get Hz

  if (max(AP_count_all_sweeps$`AP Count`) < 2) {
    `Peak to peak Frequency_Max` <- NA
    `Interevent Interval_1` <- NA
    `Interevent Interval_last` <- NA
    `SFA (1st ISI/last ISI)` <- NA

  } else if (max(AP_count_all_sweeps$`AP Count`) == 2) {

    `Peak to peak df` <- all_APs_all_sweeps %>%
      group_by(Sweep) %>%
      arrange(Sweep, `Time (ms)`) %>%
      mutate(`Peak-to-Peak Frequency (Hz)` = 1000/(`Time (ms)` - lag(`Time (ms)`))) %>%
      ungroup() %>%
      select(Sweep, `Peak-to-Peak Frequency (Hz)`) %>%
      drop_na() %>%
      slice_max(n = 1, order_by = `Peak-to-Peak Frequency (Hz)`)

    `Peak to peak Frequency_Max` <- `Peak to peak df` %>%
      pull(`Peak-to-Peak Frequency (Hz)`) %>%
      head(1)

    `Sweep # for Peak to peak Frequency_Max` <- `Peak to peak df` %>%
      pull(Sweep) %>%
      min()

    `Interevent Interval_1` <- NA
    `Interevent Interval_last` <- NA
    `SFA (1st ISI/last ISI)` <- NA

  } else if (max(AP_count_all_sweeps$`AP Count`) >= 3) {

    `Peak to peak df` <- all_APs_all_sweeps %>%
      group_by(Sweep) %>%
      arrange(Sweep, `Time (ms)`) %>%
      mutate(`Peak-to-Peak Frequency (Hz)` = 1000/(`Time (ms)` - lag(`Time (ms)`))) %>%
      ungroup() %>%
      select(Sweep, `Peak-to-Peak Frequency (Hz)`) %>%
      drop_na() %>%
      slice_max(n = 1, order_by = `Peak-to-Peak Frequency (Hz)`)

    `Peak to peak Frequency_Max` <- `Peak to peak df` %>%
      pull(`Peak-to-Peak Frequency (Hz)`) %>%
      head(1)

    `Sweep # for Peak to peak Frequency_Max` <- `Peak to peak df` %>%
      pull(Sweep) %>%
      min()

    `Interevent Interval_1` <- all_APs_all_sweeps %>%
      group_by(Sweep) %>%
      filter(Sweep == `Sweep # for SFA`) %>%
      ungroup() %>%
      mutate(`Interevent Interval (ms)` = (`Time (ms)` - lag(`Time (ms)`))) %>%
      drop_na(`Interevent Interval (ms)`) %>%
      pull(`Interevent Interval (ms)`) %>% 
      head(1)

    `Interevent Interval_last` <- all_APs_all_sweeps %>%
      group_by(Sweep) %>%
      filter(Sweep == `Sweep # for SFA`) %>%
      ungroup() %>%
      mutate(`Interevent Interval (ms)` = (`Time (ms)` - lag(`Time (ms)`))) %>%
      drop_na(`Interevent Interval (ms)`) %>%
      pull(`Interevent Interval (ms)`) %>% 
      tail(1)

    `SFA (1st ISI/last ISI)` <- (`Interevent Interval_1`/`Interevent Interval_last`)

  }

  maxAHPdeltathreshold <- `maxAHP_1stAP Voltage (mV)` - `Threshold Em (mV)`

  fAHP_delta_threshold <- `fAHP_1stAP Voltage (mV)` - `Threshold Em (mV)`

  mAHP_delta_threshold <- mAHP_1stAP- `Threshold Em (mV)`

  sAHP_delta_threshold <- sAHP_1stAP - `Threshold Em (mV)`

  `SFA (1st ISI/last ISI)` <- (`Interevent Interval_1`/`Interevent Interval_last`)

  wide_format_AP_count_all <- AP_count_all_sweeps %>%
    pivot_wider(names_from = Sweep,
                values_from = `AP Count`)

  colnames(wide_format_AP_count_all) <- paste("Sweep", colnames(wide_format_AP_count_all), "AP Count", sep = " ")

  all_files_APs_all_sweeps <- rbind(all_files_APs_all_sweeps, all_APs_all_sweeps)

  finaldf <- wide_format_AP_count_all %>%
    mutate(ObsID = this_identifier$ObsID,
           MouseID = this_identifier$MouseID,
           CellID = as.character(this_identifier$CellID),
           Filename_IV = this_identifier$Filename_IV,
           `Sweep 9 Resting Em (mV)` = `Resting Em`,
           `Sweep 9 -50 pA step Steady-State Em` = `-50 pA step Steady-State Em`,
           `Sweep 9 Resistance_in (MegaOhms)` = `Rin (MOhm)`,
           .before = `Sweep 10 AP Count`) %>%
    mutate("Sweep # for 1st AP fired" = `Sweep # for 1st AP fired`,
           "Rheobase (pA)" = `Rheobase (pA)`,
           "Threshold Em (mV)" = `Threshold Em (mV)`,
           "Peak Voltage (mV)" = `Peak Voltage (mV)`,
           "AP Amplitude (mV)" = `AP Amplitude (mV)`,
           "Half Amplitude (mV)" = `Half Amplitude (mV)`,
           "HalfWidth (ms)" = `Halfwidth (ms)`,
           "max AHP 1stAP (mV)" = `maxAHP_1stAP Voltage (mV)`,
           "max AHP delta Threshold"= maxAHPdeltathreshold,
           "max AHP Time Relative to Threshold" = `maxAHP_1stAP Time Relative to Threshold`,
           "fAHP_1stAP Voltage (mV)" = `fAHP_1stAP Voltage (mV)`,
           "fAHP delta Threshold" = fAHP_delta_threshold,
           "fAHP Time Relative to Threshold" = `fAHP_1stAP Time Relative to Threshold`,
           "mAHP_1stAP (mV)" = mAHP_1stAP,
           "mAHP delta Threshold" = mAHP_delta_threshold,
           "sAHP_1stAP" = sAHP_1stAP,
           "sAHP delta Threshold" = sAHP_delta_threshold,
           "Peak to peak Frequency_Max" = as.numeric(`Peak to peak Frequency_Max`),
           "Sweep # for Peak to peak Frequency_Max" = `Sweep # for Peak to peak Frequency_Max`,
           "Interevent Interval_1" = `Interevent Interval_1`,
           "Interevent Interval_last" = `Interevent Interval_last`,
           "SFA (1st ISI/last ISI)" = `SFA (1st ISI/last ISI)`,
           "Sweep # for SFA" = `Sweep # for SFA`)
} else {finaldf <- tibble(ObsID = this_identifier$ObsID,
                          MouseID = this_identifier$MouseID,
                          CellID = as.character(this_identifier$CellID),
                          Filename_IV = this_identifier$Filename_IV,
                          `Sweep 9 Resting Em (mV)` = `Resting Em`,
                          `Sweep 9 -50 pA step Steady-State Em` = `-50 pA step Steady-State Em`,
                          `Sweep 9 Resistance_in (MegaOhms)` = `Rin (MOhm)`,
                          "Sweep # for 1st AP fired" = `Sweep # for 1st AP fired`,
                          "Rheobase (pA)" = `Rheobase (pA)`,
                          "Threshold Em (mV)" = `Threshold Em (mV)`,
                          "Peak Voltage (mV)" = `Peak Voltage (mV)`,
                          "AP Amplitude (mV)" = `AP Amplitude (mV)`,
                          "Half Amplitude (mV)" = `Half Amplitude (mV)`,
                          "HalfWidth (ms)" = `Halfwidth (ms)`,
                          "max AHP 1stAP (mV)" = `maxAHP_1stAP Voltage (mV)`,
                          "max AHP delta Threshold"= maxAHPdeltathreshold,
                          "max AHP Time Relative to Threshold" = `maxAHP_1stAP Time Relative to Threshold`,
                          "fAHP_1stAP Voltage (mV)" = `fAHP_1stAP Voltage (mV)`,
                          "fAHP delta Threshold" = fAHP_delta_threshold,
                          "fAHP Time Relative to Threshold" = `fAHP_1stAP Time Relative to Threshold`,
                          "mAHP_1stAP (mV)" = mAHP_1stAP,
                          "mAHP delta Threshold" = mAHP_delta_threshold,
                          "sAHP_1stAP" = sAHP_1stAP,
                          "sAHP delta Threshold" = sAHP_delta_threshold,
                          "Peak to peak Frequency_Max" = `Peak to peak Frequency_Max`,
                          "Sweep # for Peak to peak Frequency_Max" = `Sweep # for Peak to peak Frequency_Max`,
                          "Interevent Interval_1" = `Interevent Interval_1`,
                          "Interevent Interval_last" = `Interevent Interval_last`,
                          "SFA (1st ISI/last ISI)" = `SFA (1st ISI/last ISI)`,
                          "Sweep # for SFA" = `Sweep # for SFA`)
}

all_data <- rbind.fill(all_data, finaldf) %>%
  arrange(ObsID) %>%
  relocate(ends_with("AP Count"),
           .before = `Sweep # for 1st AP fired`)

count_obs <- NULL
tmp_pl <- NULL
pdfname <- NULL
for (sweep in unique(all_APs_all_sweeps$Sweep)) {
  count_obs <- length(all_APs_all_sweeps$`Voltage (mV)`[all_APs_all_sweeps$Sweep == sweep])

  tmp_pl <- df1 %>%
    group_by(Sweep) %>%
    filter(any(Sweep == sweep)) %>%
    ggplot(aes(x = `Time (ms)` ,
               y = `Voltage (mV)`),
           lab) +
    geom_line(linewidth = 0.2,
              show.legend = FALSE) +
    geom_point(data = all_APs_all_sweeps %>% filter(Sweep == sweep),
               size = 2,
               show.legend = FALSE,
               aes(x = `Time (ms)`,
                   y = `Voltage (mV)`,
                   color = "#FB8072")) +
    geom_text_repel(data = all_APs_all_sweeps %>% filter(Sweep == sweep),
                             size = 2,
                             show.legend = FALSE,
                             aes(x = `Time (ms)`,
                                 y = `Voltage (mV)`,
                                 label = `Voltage (mV)`)) +
    theme(legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, face = "bold", hjust = c(0,1), colour = c("black", "#FB8072")),
          axis.title = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5)) +
    labs(subtitle = c(paste("Sweep Number:", sweep), paste("Total Number of Peaks:", count_obs))) +
    ggtitle(get("filename"))
  pdfname <- paste(get("filename"), get("sweep"), sep = "_Sweep_") %>% paste(".pdf", sep = "")
  ggsave(filename = pdfname, plot = tmp_pl, path = "Plots", width = 10, height = 6.5)
}

write_csv(finaldf,
          paste("Computer Data Analysis ",
                get("filename"),
                ".csv",
                sep = ""))

pdf_combine(input = list.files(path = "Plots", full.names=TRUE, pattern=".pdf"),
            output = paste("All Plots ",
                           get("filename"),
                           ".pdf",
                           sep = ""))

unlink("Plots", recursive = TRUE)
