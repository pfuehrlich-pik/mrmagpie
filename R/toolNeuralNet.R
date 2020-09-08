#' Neural Network Reconstruction
#'
#'Reconstructs and evaluate a neural network from the weights and biases provided as arguments
#'
#' @param inputs_ml Neural Network input features properly scaled with
#' the scale and center attributes of the scaled training set.
#' @param weights The learned weights and biases in a list format as outputed by the function \code{keras::get_weights()}. This function works with magpie objects
#' @return The evaluated result of the neural network for the \code{input_ml} parameter.
#' @author Marcos Alves
#' @export toolNeuralNet

toolNeuralNet <- function(inputs_ml, weights) {
  x <- paste0(".f <- function(input) {")
  x <- append(x, paste0("y <- {input %*% weights[[1]] + t(weights[[2]])} %>%"))
  x <- append(x, paste0("{log(1+exp(.))} %>%"))
  for (i in seq(3, length(weights), 2)) {
    if (i < length(weights) - 2) {
      x <- append(x, paste0("{. %*% weights[[", i, "]] + t(weights[[", i + 1, "]])} %>%"))
      x <- append(x, paste0("{log(1+exp(.))} %>%"))
    } else {
      x <- append(x, paste0("{. %*% weights[[", i, "]] + t(weights[[", i + 1, "]])}"))
    }
  }
  x <- append(x, "return(y)}")
  func <- eval(parse(text = x))
  out <- apply(inputs_ml, 1, func)
  return(out)
}
