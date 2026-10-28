test_that("the package loads and documents itself", {
  expect_true("flexsae" %in% loadedNamespaces())
  expect_s3_class(utils::packageDescription("flexsae"), "packageDescription")
})
