# generates a confusion matrix for 8 variable
generate_conf_matrix <- function(true_labels, predicted_labels){
  conf_matrix <- confusionMatrix(predicted_labels, true_labels)

  cm_df <- as.data.frame(conf_matrix$table)

  # Plot the confusion matrix using ggplot2
  ggplot(cm_df, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq), color = "black", size = 5) +
    scale_fill_gradient(low = "lightyellow", high = "purple") +
    theme_minimal() +
    labs(title = "Confusion Matrix", x = "Predicted Class", y = "True Class") +
    scale_x_discrete(limits = as.character(1:8)) +
    scale_y_discrete(limits = as.character(1:8))
}
