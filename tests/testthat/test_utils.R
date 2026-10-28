test_that(".with_seed restores the stream, or leaves it absent", {
  set.seed(5)
  before <- .Random.seed
  a <- .with_seed(1, runif(3))
  expect_identical(.Random.seed, before)
  expect_identical(a, .with_seed(1, runif(3)))

  rm(".Random.seed", envir = globalenv())
  .with_seed(2, runif(1))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(5)
})

test_that(".fh_formula replaces the response and adds the area term", {
  f <- .fh_formula(log(y) ~ s(x1) + x2)
  expect_identical(deparse(f), ".flexsae_y ~ s(x1) + x2 + s(.flexsae_area, bs = \"re\")")
})

test_that(".get_variable rejects a vector of the wrong length", {
  expect_error(.get_variable(data.frame(a = 1:3), 1:2, "vardir"), "length 3")
})
