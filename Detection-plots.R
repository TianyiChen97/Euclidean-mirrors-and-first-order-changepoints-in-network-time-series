analyze_power_curve <- function(ts_data, p_string, plot_hist = TRUE, plot_power = TRUE) {
  
  # 1. Construct the column name for the critical value (cv)
  cv_col_name <- paste0("q_", p_string)
  
  # Check if this column exists in the data
  if (!cv_col_name %in% names(ts_data)) {
    stop(paste("Error: Column", cv_col_name, "not found in the provided data."))
  }
  
  # 2. Calculate the critical value (95th percentile)
  message(paste("Calculating CV from column:", cv_col_name))
  cv <- quantile(ts_data[[cv_col_name]], 0.95, na.rm = TRUE)
  
  # 3. [Optional] Plot histogram
  if (plot_hist) {
    hist(ts_data[[cv_col_name]],
         breaks = 50,
         main = paste("Histogram for", cv_col_name, "(CV =", round(cv, 3), ")"),
         xlab = "Test Statistic Value")
    abline(v = cv, col = "red", lwd = 2)
  }
  
  # 4. Calculate power values AND error bars
  # Sapply will return a matrix where each column is a q_value
  # and the rows are 'power', 'lower_bound', 'upper_bound'
  power_stats_matrix <- sapply(ts_data, function(stats_vector) {
    
    # Create the binary vector (1 = reject, 0 = fail to reject)
    binary_vector <- stats_vector > cv
    
    # Get nmc (handle NAs)
    nmc <- sum(!is.na(binary_vector))
    
    if (nmc == 0) {
      return(c(power = NA, lower_bound = NA, upper_bound = NA))
    }
    
    # Calculate power (the mean)
    p <- mean(binary_vector, na.rm = TRUE)
    
    # Calculate standard deviation of the binary vector
    sd_binary <- sd(binary_vector, na.rm = TRUE)
    
    # If p=0 or p=1, sd_binary will be NA. In this case, the error is 0.
    if (is.na(sd_binary)) {
      sd_binary <- 0
    }
    
    # Calculate margin of error as per your formula
    margin <- 1.96 * sd_binary / sqrt(nmc)
    
    # Calculate bounds
    lower <- p - margin
    upper <- p + margin
    
    # Return the three statistics
    c(power = p, lower_bound = lower, upper_bound = upper)
  })
  
  # Convert the matrix to a more usable data frame
  # Transpose (t()) so rows are q_values and columns are stats
  power_stats_df <- as.data.frame(t(power_stats_matrix))
  
  
  # 5. Extract and process x-axis values from the names
  q_names <- rownames(power_stats_df) # Get names from row names
  
  if (any(grepl("p", q_names))) {
    message("Detected 'p' in names. Using 'p' as decimal.")
    x_q_values <- as.numeric(sub("p", ".", sub("q_", "", q_names)))
  } else {
    message("No 'p' detected. Using standard numeric conversion.")
    x_q_values <- as.numeric(sub("q_", "", q_names))
  }
  
  # 6. Order the values for plotting
  plot_data <- cbind(
    x_q_values = x_q_values,
    power_stats_df # This now has 'power', 'lower_bound', 'upper_bound'
  )
  plot_data_ordered <- plot_data[order(plot_data$x_q_values), ]
  
  # 7. [Optional] Plot the power curve
  if (plot_power) {
    # Set ylim to ensure error bars fit
    ylim_min <- min(0, plot_data_ordered$lower_bound, na.rm = TRUE)
    ylim_max <- max(1, plot_data_ordered$upper_bound, na.rm = TRUE)
    
    plot(plot_data_ordered$x_q_values, plot_data_ordered$power,
         type = "b", # "b" for both points and lines
         pch = 19,
         xlab = "q",
         ylab = "Power",
         main = paste("Power Curve (CV from", p_string, ")"),
         ylim = c(ylim_min, ylim_max) # Apply ylim
    )
    
    # --- Add Error Bars ---
    arrows(
      x0 = plot_data_ordered$x_q_values,
      y0 = plot_data_ordered$lower_bound,
      x1 = plot_data_ordered$x_q_values,
      y1 = plot_data_ordered$upper_bound,
      length = 0.03, # Length of the cross-bar
      angle = 90,    # 90 degrees makes it a flat bar
      code = 3       # Code 3 means draw bars at both ends
    )
    
    ppp <- as.numeric(sub("p", ".", p_string))
    
    # Add vertical lines
    abline(v = ppp, lwd = 1.5)
    abline(v = ppp - 0.1, col = "gray60", lty = 2)
    abline(v = ppp + 0.1, col = "gray60", lty = 2)
  }
  
  # 8. Return the ordered data (now with error bounds)
  return(plot_data_ordered)
}


load(".../all_n_varying_results_m16_nmc300_tstar8_p0p4_finished_20251023_233709.Rdata")
all_n_results_graph <- all_n_results
rm(all_n_results)


library(dplyr)
library(ggplot2)
library(tidyverse)


dodge_width <- 0

all_n500_graph <- all_n_results_graph$n_500
df_hist <- data.frame(val = as.numeric(unlist(all_n500_graph$q_0p4)))

ggplot(df_hist, aes(x = val)) +
  geom_histogram(bins = 50, fill = 'black', color = "white") +
  geom_vline(xintercept = cv, color = "red", size = 1) +
  labs(x = "test statistics", y = "Frequency") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )





cv = quantile(all_n500_graph$q_0p4,0.95)



all_n1000_graph <-   all_n_results_graph$n_1000
all_n1500_graph <-   all_n_results_graph$n_1500


data1 <- analyze_power_curve(all_n500_graph, '0p4')
data2 <- analyze_power_curve(all_n1000_graph, '0p4')
data3 <- analyze_power_curve(all_n1500_graph, '0p4')

# 1. Combine the three data frames into one.
# We create a new column 'n_group' to label the data.
all_data_graph <- bind_rows(
  data1 %>% mutate(n_group = "n = 500"),
  data2 %>% mutate(n_group = "n = 1000"),
  data3 %>% mutate(n_group = "n = 1500")
) %>%
  # 2. Convert 'n_group' to a factor.
  # This ensures the legend is ordered correctly (500, 1000, 1500)
  # instead of alphabetically (1000, 1500, 500).
  mutate(n_group = factor(n_group, levels = c("n = 500", "n = 1000", "n = 1500"))) %>%
  mutate(cat = factor("graph"))




filter_data = all_data_graph[all_data_graph$n_group == "n = 1500" |all_data_graph$n_group == "n = 500" ,]
filter_data = all_data_graph[all_data_graph$n_group == "n = 500" ,]

ggplot(filter_data, 
       aes(x = x_q_values, 
           y = power)) +
  
  geom_errorbar(
    aes(ymin = lower_bound, ymax = upper_bound),
    width = 0.002, # Sets the width of the horizontal "caps"
    position = position_dodge(dodge_width) # Dodges the bars
  ) +
  
  geom_line(
    linewidth = 1.1, 
    position = position_dodge(dodge_width)
  ) +
  geom_point(
    size = 2.5, 
    position = position_dodge(dodge_width)
  ) +
  
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "black", alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red", alpha = 0.7) +
  
  labs(
    x = "q",
    y = "empirical power",
  ) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  
  # --- Theme and Legend (as before) ---
  theme_minimal() +
  theme(
    # --- Increase Axis Font Sizes ---
    axis.title = element_text(size = 16, face = "bold"), # Large labels ("q", "empirical power")
    axis.text = element_text(size = 14),               # Large numbers on the ticks
    
    # --- Legend and Layout ---
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank()                 # Cleans up the look for larger text
  )
