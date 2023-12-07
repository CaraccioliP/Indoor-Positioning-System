#Clean data 

#create function
CleanFunct = function(df)
{
  variables = strsplit(df, "[;=,]")[[1]]
  #to fix an error -> one anomalous obsv. We can also remove it by hand.
  if (length(variables) == 10)
    return (NULL)
  #select composite variables for the device (pos, maciDs) -> the ones that have multiple data in one line
  comp_var = matrix(variables[ - (1:10)], ncol=4, byrow=TRUE)
  cbind(matrix(variables[c(2,4,6,7,8,10)], nrow(comp_var), 6, byrow=TRUE), comp_var)
}


# Function to round to nearest 45 degree angle
orient_round = function(orient) 
{
  binned_angles = seq(0, by = 45, length = 9)
  q = sapply(orient, function(theta) which.min(abs(theta - binned_angles)))
  c(binned_angles[1:8], 0) [q]
}


# This function cleans the df, and removes unwanted variables.
Clean_all = function(filename, macnames)
{
  #open txt file
  txt = readLines(filename)
  #eliminate the #
  lines = txt[ substr(txt, 1, 1) != "#" ]
  #clean the textfile
  tmpMat = lapply(lines, CleanFunct)
  df = as.data.frame(do.call("rbind", tmpMat), 
                          stringsAsFactors= FALSE) 
  #add names to df
  names(df) = c("time", "scanMac", 
                     "posX", "posY", "posZ", "orientation", 
                     "mac", "signal", "channel", "type")
  
  # keep only signals from access points
  df = df[ offline$type == "3", ]
  
  # drop scanMac, posZ, channel, and type
  dropVars = c("scanMac", "posZ", "channel", "type")
  df = df[ , !( names(df) %in% dropVars ) ]
  
  # convert numeric values
  numVars = c("time", "posX", "posY", "orientation", "signal")
  df[numVars] = lapply(df[numVars], as.numeric)
  
  # convert time to POSIX
  df$original_time = df$time
  df$time = df$time/1000
  class(df$time) = c("POSIXt", "POSIXct")
  
  # round orientations to nearest 45 degree reference angle
  df$rounded_orient = orient_round(df$orientation)
  
  # drop more unwanted access points
  df = df[ df$mac %in% macnames, ]
  
  return(df)
}


mac_locs <- readr::read_table("accessPointLocations.txt", col_names = c("Macs", "posX", "posY"))
mac_locs = as.data.frame(mac_locs)
#Extract the Macs name from mac_locs
macnames = mac_locs$Macs


#this results in a clean filtered df
offline = Clean_all("offline.final.trace.txt", macnames)
online = Clean_all("online.final.trace.txt", macnames)

# 
# 
# #organize the fd by xy combination
# offline$pos_XY = paste(offline$posX, offline$posY, sep = "-")
# 
# byLocAngleAP = with(offline, by(offline, list(pos_XY, angle, mac), function(x) x))
# 
# signalSummary = lapply(byLocAngleAP, function(oneLoc) {
#   ans = oneLoc[1, ]
#   ans$medSignal = median(oneLoc$signal)
#   ans$avgSignal = mean(oneLoc$signal)
#   ans$num = length(oneLoc$signal)
#   ans$sdSignal = sd(oneLoc$signal)
#   ans$iqrSignal = IQR(oneLoc$signal)
#   ans
# })
# 
# 
# offlineSummary = do.call("rbind", signalSummary)
# 
# 
# online$posXY = paste(online$posX, online$posY, sep = "-")
# keepVars = c("pos_XY", "posX","posY", "orientation", "angle")
# 
# byLoc = with(online, 
#              by(online, list(pos_XY), 
#                 function(x) {
#                   ans = x[1, keepVars]
#                   avgSS = tapply(x$signal, x$mac, mean)
#                   y = matrix(avgSS, nrow = 1, ncol = 6,
#                              dimnames = list(ans$posXY, names(avgSS)))
#                   cbind(ans, y)
#                 }))
# 
# byLocof = with(offline, 
#                by(offline, list(pos_XY), 
#                   function(x) {
#                     ans = x[1, keepVars]
#                     avgSS = tapply(x$signal, x$mac, mean)
#                     y = matrix(avgSS, nrow = 1, ncol = 6,
#                                dimnames = list(ans$pos_XY, names(avgSS)))
#                     cbind(ans, y)
#                   }))
# 
# 
# onlineSummary = do.call("rbind", byLoc)
# offlineSummary = do.call('rbind', byLocof)

# write.csv(offlineSummary, "C:\\Users\\Pascal\\Desktop\\PSU DS\\Fall 2023\\STAT 410\\Indoor-Positioning-System\\offline_summary.csv", row.names=TRUE)
# 
# write.csv(onlineSummary, "C:\\Users\\Pascal\\Desktop\\PSU DS\\Fall 2023\\STAT 410\\Indoor-Positioning-System\\online_summary.csv", row.names=TRUE)
