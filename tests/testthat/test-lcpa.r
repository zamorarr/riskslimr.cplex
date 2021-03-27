titanic <- readRDS("titanic.rds")
x <- titanic$x
y <- titanic$y

test_that("titanic dataset works", {
  # actual
  r <- lcpa_cpp(x, y, logfile = tempfile(), R_max = 3L, time_limit = 10L)

  # expected
  alpha <- as.integer(c(1,1,1,0,0,0))
  lambda <- as.integer(c(-3, 1, 2, 0, 0, 0))
  optimality_gap <- 0

  # test
  expect_identical(r$alpha, alpha)
  expect_identical(r$lambda, lambda)
  expect_identical(r$optimality_gap, optimality_gap)
})
