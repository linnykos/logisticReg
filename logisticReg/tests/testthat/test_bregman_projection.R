context("Test bregman projections")

## .projection_bregman is correct

test_that(".projection_bregman works", {
  z <- rep(.5, 2)
  f <- .conjugate_gaussian; grad_f <- .conjugate_grad_gaussian
  plane <- .plane(basis = matrix(c(.7,.3), nrow = 2, ncol = 1))

  res <- .projection_bregman(z, f, grad_f, plane)
})
