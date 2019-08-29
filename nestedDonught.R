library('plotly')
library('dplyr')

data <- data.frame(c('cats', 'monkeys', 'dogs'), c(30, 10, 20), c(20, 10, 10))
colnames(data) <- c('animal', 'street', 'home')

p <- plot_ly(data) %>%
  add_pie(labels = ~animal, values = ~street, type = 'pie', hole = 0.75, sort = F) %>%
  add_pie(data, labels = ~animal, values = ~home, 
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F, hole = 0.3)
  
p