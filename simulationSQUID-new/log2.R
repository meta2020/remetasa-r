logfile <- "output.txt"

sink(logfile)
msg_con <- file(logfile, open = "a")
sink(msg_con, type = "message")

# on.exit({
#   sink(type = "message")
#   close(msg_con)
#   sink()
# }, add = TRUE)

message("this is a message")


sink(type = "message")
sink()
