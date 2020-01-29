table <- read.table('wing_area_calculations.txt', h=T)

table$R_difference <- abs(table$areaM_E-table$areaM_M)

head(table)

tableN <- c(table$A_difference, table$R_difference)

tableN <- cbind(tableN, c(rep('Absolute difference', nrow(table)), rep('Size difference', nrow(table))))

tableN <- as.data.frame(tableN)

colnames(tableN) <- c('value', 'type')
tableN$value <- as.numeric(as.character(tableN$value))

library(ggplot2)

box_plot <- ggplot(tableN, aes(x = type, y = value))

box_plot +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               dotsize = 1,
               stackdir = 'center') +
  theme_classic()

# box_plot +
#   geom_boxplot() +
#   geom_jitter(shape = 19,
#               color = "steelblue",
#               cex=3,
#               position = position_jitter(width = 0.21)) +
#   theme_classic()
