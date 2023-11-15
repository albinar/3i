# Written by Albina Rahim - Last updated on September 2019
# Updated by Marjan Najafabadipour on November 2019

#flowGroup helper functions

#################################################################################
# Prints out the time since start_time. Used for optimizing code.
TimeOutput <- function(start_time) {
    start_time <- as.POSIXct(start_time)
    dt <- difftime(Sys.time(), start_time, units="secs")
    # Since you only want the H:M:S, we can ignore the date...
    # but you have to be careful about time-zone issues
    format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}


load_flowFrame <- function(Fcs_file_name){
  flowFrame <- read.FCS(Fcs_file_name)
    return(flowFrame)
}


find_number <- function(y, order_of_plotting){
    location <- unlist(lapply(order_of_plotting, function(x) {
        temp <- which(x == y)
        if (length(temp) == 0){
            temp <- "nope"
        }
        return(temp)
    }))
    names(location) <- 1:length(location)
    location <- location[which(location != "nope")]
    return(location)
}