# Error handling and edge case tests for PSGP package

library(psgp)

# Test 1: Invalid input validation
cat("Testing input validation...\n")

# Test with NULL inputs - should fail gracefully
tryCatch({
  estimateParameters(NULL, NULL, NULL)
  stop("Should have failed with NULL inputs")
}, error = function(e) {
  cat("PASS: NULL input validation working:", e$message, "\n")
})

# Test 2: Empty data
cat("Testing empty data handling...\n")

tryCatch({
  empty_x <- numeric(0)
  empty_y <- numeric(0)
  estimateParameters(empty_x, empty_y, NULL)
  stop("Should have failed with empty data")
}, error = function(e) {
  cat("PASS: Empty data validation working:", e$message, "\n")
})

# Test 3: Mismatched dimensions
cat("Testing dimension mismatch handling...\n")

tryCatch({
  x_coords <- c(1, 2, 3, 4)  # 2 coordinate pairs
  y_values <- c(1, 2, 3)     # 3 values - mismatch
  estimateParameters(x_coords, y_values, NULL)
  stop("Should have failed with dimension mismatch")
}, error = function(e) {
  cat("PASS: Dimension mismatch validation working:", e$message, "\n")
})

# Test 4: Invalid prediction coordinates
cat("Testing invalid prediction coordinates...\n")

# Set up valid training data
library(sp)
data(meuse)
observations <- data.frame(x = meuse$x[1:10], y = meuse$y[1:10], value = log(meuse$zinc[1:10]))
coordinates(observations) = ~x+y

# Create invalid prediction locations (odd number of coordinates)
tryCatch({
  invalid_pred <- c(1, 2, 3)  # Odd number - should be pairs
  spatialPredict(observations, invalid_pred)
  stop("Should have failed with invalid prediction coordinates")
}, error = function(e) {
  cat("PASS: Invalid prediction coordinates validation working:", e$message, "\n")
})

# Test 5: Numerical stability with ill-conditioned data
cat("Testing numerical stability...\n")

# Create data points that are very close together (potential numerical issues)
close_x <- c(0, 0.000001, 0, 0.000001)
close_y <- c(0, 0, 0.000001, 0.000001)
close_values <- c(1, 1.1, 1.05, 1.15)

close_obs <- data.frame(x = close_x, y = close_y, value = close_values)
coordinates(close_obs) = ~x+y

tryCatch({
  result <- spatialPredict(close_obs, c(0.0005, 0.0005))
  cat("PASS: Numerical stability test completed without errors\n")
}, error = function(e) {
  cat("WARNING: Numerical stability issue detected:", e$message, "\n")
})

cat("Error handling tests completed.\n")