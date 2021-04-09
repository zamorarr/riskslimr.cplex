test_that("loss and loss_grad work with simple example", {
  x <- matrix(c(1, 5), nrow = 1)
  y <- c(-1)
  c0 <- 0
  weights <- c(1) # n x 1
  lambda <- c(2, 1)

  # test
  loss <- test_loss_cpp(x, y, weights, c0, lambda)
  loss_grad <- as.vector(test_loss_grad_cpp(x, y, weights, c0, lambda))
  expect_equal(loss, log(1 + exp(7)))
  expect_equal(loss_grad, c(1,5)/(1 + exp(-7)))
})

test_that("loss and loss_grad work with weights", {
  x <- matrix(c(1, 10, 100, 5, 50, 500), nrow = 3) # n x d
  y <- c(-1, 1, -1) # n x 1
  c0 <- 0
  weights <- c(1, 0, 0) # n x 1
  lambda <- c(2, 1) # d x 1

  # test
  loss <- test_loss_cpp(x, y, weights, c0, lambda)
  loss_grad <- as.vector(test_loss_grad_cpp(x, y, weights, c0, lambda))
  expect_equal(loss, log(1 + exp(7)))
  expect_equal(loss_grad, c(1,5)/(1 + exp(-7)))
})
