library(codetools)
library(lattice)
library(fields)

# Read data using readLines
offl <- file("offline.final.trace.txt", "r")
txt <- readLines(offl)
sum(substr(txt, 1, 1) == "#")
length(txt)
strsplit(txt[4], ";")[[1]]
unlist(lapply(strsplit(txt[4], ";")[[1]],
              function(x)
                sapply(strsplit(x, "=")[[1]], strsplit, ",")))
tokens <- strsplit(txt[4], "[;=,]")[[1]]
tokens[1:10]
# Extract values of variables
tokens[c(2, 4, 6:8, 10)]
tokens[ - (1:10)]

tmp <- matrix(tokens[ - (1:10 ) ], ncol = 4, byrow = TRUE)
mat <- cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp),
                    ncol = 6, byrow = TRUE), tmp)
# Check
dim(mat)
processLine =
  function(x)
  {
    tokens = strsplit(x, "[;=,]")[[1]]
    tmp = matrix(tokens[ - (1:10) ], ncol = 4, byrow = TRUE)
    cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp),
                 ncol = 6, byrow = TRUE), tmp)
  }
tmp = lapply(txt[4:20], processLine)
sapply(tmp, nrow)
offline = as.data.frame(do.call("rbind", tmp))
# Check it
dim(offline)
lines <- txt[ substr(txt, 1, 1) != "#" ]
tmp = lapply(lines, processLine)
options(error = recover, warn = 2)
tmp = lapply(lines, processLine)
processLine = function(x)
{
  tokens = strsplit(x, "[;=,]")[[1]]
  
  if (length(tokens) == 10) 
    return(NULL)
  
  tmp = matrix(tokens[ - (1:10) ], , 4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow(tmp), 6, 
               byrow = TRUE), tmp)
}

options(error = recover, warn = 1)
tmp <- lapply(lines, processLine)
offline <- as.data.frame(do.call("rbind", tmp), 
                         stringsAsFactors = FALSE)

dim(offline)

names(offline) = c("time", "scanMac", "posX", "posY", "posZ", 
                   "orientation", "mac", "signal", 
                   "channel", "type")
numVars = c("time", "posX", "posY", "posZ", 
            "orientation", "signal")
offline[ numVars ] =  lapply(offline[ numVars ], as.numeric)
offline = offline[ offline$type == "3", ]
offline = offline[ , "type" != names(offline) ]
dim(offline)
offline$rawTime = offline$time
offline$time = offline$time/1000
class(offline$time) = c("POSIXt", "POSIXct")
unlist(lapply(offline, class))
summary(offline[, numVars])
summary(sapply(offline[ , c("mac", "channel", "scanMac")], as.factor))
offline <- offline[ , !(names(offline) %in% c("scanMac", "posZ"))]
length(unique(offline$orientation))
roundOrientation = function(angles) {
  refs = seq(0, by = 45, length  = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}
# We now use our function to created the rounded angles
offline$angle <- roundOrientation(offline$orientation)

c(length(unique(offline$mac)), length(unique(offline$channel)))

subMacs <- names(sort(table(offline$mac), decreasing = TRUE))[1:7]
offline <- offline[ offline$mac %in% subMacs, ]
macChannel <- with(offline, table(mac, channel))
apply(macChannel, 1, function(x) sum(x > 0))
offline <- offline[ , "channel" != names(offline)]

locDF <- with(offline,
              by(offline, list(posX, posY), function(x) x))
# Check
length(locDF)

sum(sapply(locDF, is.null))
locDF <- locDF[ !sapply(locDF, is.null) ]
# Check
length(locDF)
locCounts <- sapply(locDF, nrow)

# If we want to keep the position information with the lcoation, we do this with
locCounts <- sapply(locDF,
                    function(df)
                      c(df[1, c("posX", "posY")], count = nrow(df)))
# Confirm we have a matrix
class(locCounts)
# Confirm 3 rows
dim(locCounts)
# Examine a few of the counts
locCounts[ , 1:8]


offline$posXY <- paste(offline$posX, offline$posY, sep = "-")
byLocAngleAP <- with(offline, 
                     by(offline, list(posXY, angle, mac), 
                        function(x) x))
signalSummary = 
  lapply(byLocAngleAP,            
         function(oneLoc) {
           ans = oneLoc[1, ]
           ans$medSignal = median(oneLoc$signal)
           ans$avgSignal = mean(oneLoc$signal)
           ans$num = length(oneLoc$signal)
           ans$sdSignal = sd(oneLoc$signal)
           ans$iqrSignal = IQR(oneLoc$signal)
           ans
         })

offlineSummary = do.call("rbind", signalSummary)


oneAPAngle <- subset(offline, mac == subMacs[5] & angle == 0)
oneAPAngle = subset(offlineSummary, 
                    mac == subMacs[5] & angle == 0)

smoothSS = Tps(oneAPAngle[, c("posX","posY")], 
               oneAPAngle$avgSignal)

vizSmooth = predictSurface(smoothSS)

plot.surface(vizSmooth, type = "C")

points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)

surfaceSS = function(data, mac, angle = 45) {
  require(fields)
  oneAPAngle = data[ data$mac == mac & data$angle == angle, ]
  smoothSS = Tps(oneAPAngle[, c("posX","posY")], 
                 oneAPAngle$avgSignal)
  vizSmooth = predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C", 
               xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5) 
}

parCur = par(mfrow = c(2,2), mar = rep(1, 4))

mapply(surfaceSS, mac = subMacs[ rep(c(5, 1), each = 2) ], 
       angle = rep(c(0, 135), 2),
       data = list(data = offlineSummary))

par(parCur)
offlineSummary <- subset(offlineSummary, mac != subMacs[2])
AP <- matrix( c( 7.5, 6.3, 2.5, -.8, 12.8, -2.8,  
                 1, 14, 33.5, 9.3,  33.5, 2.8),
              ncol = 2, byrow = TRUE,
              dimnames = list(subMacs[ -2 ], c("x", "y") ))

AP


diffs <- offlineSummary[ , c("posX", "posY")] -
  AP[ offlineSummary$mac, ]
offlineSummary$dist <- sqrt(diffs[ , 1]^2 + diffs[ , 2]^2)


macs <- unique(offlineSummary$mac)
# Read data using readLines
onli <- file("online.final.trace.txt", "r")
txt <- readLines(onli)

# Process the lines
lines <- txt[substr(txt, 1, 1) != "#"]

# Filter lines based on subMacs (replace with your actual logic)
subMacs <- c("00:0f:a3:39:e1:c0", "00:0f:a3:39:dd:cd", "00:14:bf:b1:97:8a",
             "00:14:bf:3b:c7:c6", "00:14:bf:b1:97:90", "00:14:bf:b1:97:8d",
             "00:14:bf:b1:97:81")  # Replace with specific MAC addresses
filtered_lines <- grep(paste(subMacs, collapse = "|"), lines, value = TRUE)

# Continue with the rest of your code
processLine <- function(x) {
  tokens <- strsplit(x, "[;=,]")[[1]]
  tmp <- matrix(tokens[-(1:10)], ncol = 4, byrow = TRUE)
  cbind(matrix(tokens[c(2, 4, 6:8, 10)], nrow = nrow(tmp), ncol = 6, byrow = TRUE), tmp)
}
tmp <- lapply(filtered_lines, processLine)
online <- as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)

# Set column names
names(online) <- c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac", "signal", "channel", "type")

# Create a unique location identifier
online$posXY <- paste(online$posX, online$posY, sep = "-")

# Check
length(unique(online$posXY))

roundOrientation = function(angles) {
  refs = seq(0, by = 45, length  = 9)
  q = sapply(angles, function(o) which.min(abs(o - refs)))
  c(refs[1:8], 0)[q]
}

# Convert orientation to numeric
online$orientation <- as.numeric(online$orientation)

# Check for any non-numeric values
non_numeric_values <- online$orientation[is.na(online$orientation)]
if (length(non_numeric_values) > 0) {
  warning("Non-numeric values found in orientation column. Setting them to NA.")
  online$orientation[is.na(online$orientation)] <- NA
}

# Now, use the roundOrientation function
online$angle <- roundOrientation(online$orientation)

# Create the table
tabonlineXYA <- table(online$posXY, online$angle)
tabonlineXYA[1:6, ]


keepVars <- c("posXY", "posX", "posY", "orientation", "angle")

byLoc <- with(online, 
              by(online, list(posXY), 
                 function(x) {
                   ans <- x[1, keepVars]
                   avgSS <- tapply(x$signal, x$mac, mean)
                   y <- matrix(avgSS, nrow = 1, ncol = length(avgSS),
                               dimnames = list(ans$posXY, names(avgSS)))
                   cbind(ans, y)
                 }))

# Ensure that all data frames in byLoc have the same columns
colnames_list <- lapply(byLoc, colnames)
common_cols <- Reduce(intersect, colnames_list)

# Filter only common columns in each data frame
byLoc <- lapply(byLoc, function(df) df[, common_cols, drop = FALSE])

onlineSummary <- do.call("rbind", byLoc)

dim(onlineSummary)



m <- 3; angleNewObs <- 230
refs <- seq(0, by = 45, length = 8)
nearestAngle <- roundOrientation(angleNewObs)

if (m %% 2 == 1) {
  angles = seq(-45 * (m-1) /2, 45 * (m-1) /2, length = m)
} else {
  m = m+1
  angles = seq(-45 * (m-1) /2, 45 * (m-1) /2, length = m)
  if (sign(angleNewObs - nearestAngle) > -1)
    angles = angles[ -1 ]
  else
    angles = angles[ -m ]
}

## m odd and even are handled seperatly. Map angles to refs
## Adjustments
angles <- angles + nearestAngle
angles[angles < 0] = angles[ angles < 0 ] + 360
angles[angles > 360] = angles[ angles > 360 ] - 360

offlineSubset <-
  offlineSummary[ offlineSummary$angle %in% angles, ]

reshapeSS <- function(data, varSignal = "signal",
                      keepVars = c("posXY", "posX","posY")) {
  byLocation <-
    with(data, by(data, list(posXY),
                  function(x) {
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6,
                               dimnames = list(ans$posXY,
                                               names(avgSS)))
                    cbind(ans, y)
                  }))
  
  newDataSS <- do.call("rbind", byLocation)
  return(newDataSS)
}


trainSS = reshapeSS(offlineSubset, varSignal = "avgSignal")

selectTrain = function(angleNewObs, signals = NULL, m = 1){
  # m is the number of angles to keep between 1 and 5
  refs = seq(0, by = 45, length  = 8)
  nearestAngle = roundOrientation(angleNewObs)
  
  if (m %% 2 == 1) 
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  else {
    m = m + 1
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    if (sign(angleNewObs - nearestAngle) > -1) 
      angles = angles[ -1 ]
    else 
      angles = angles[ -m ]
  }
  angles = angles + nearestAngle
  angles[angles < 0] = angles[ angles < 0 ] + 360
  angles[angles > 360] = angles[ angles > 360 ] - 360
  angles = sort(angles) 
  
  offlineSubset = signals[ signals$angle %in% angles, ]
  reshapeSS(offlineSubset, varSignal = "avgSignal")
}

train130 <- selectTrain(130, offlineSummary, m = 3)
head(train130)
# Should be 166
length(train130[[1]])


findNN <- function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1, 
                function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)) )
  closest = order(dists)
  return(trainSubset[closest, 1:3 ])
}

predXY <- function(newSignals, newAngles, trainData, 
                   numAngles = 1, k = 3){
  
  closeXY = list(length = nrow(newSignals))
  
  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]),
                          trainSS)
  }
  
  estXY = lapply(closeXY, function(x)
    sapply(x[ , 2:3], 
           function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

