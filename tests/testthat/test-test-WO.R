test_that("WO for KHCE, SAS", {
  expect_equal(calcWO(as_hce(KHCE))$WO, 1.330881591215)
})

test_that("WP for KHCE, SAS", {
  expect_equal(calcWO(as_hce(KHCE))$WP, 0.570977777778)
})

test_that("SE for WP for KHCE, SAS", {
  expect_equal(calcWO(as_hce(KHCE))$SE_WP, 0.014734461746)
})

test_that("p-value for WP for KHCE, SAS", {
  expect_equal(calcWO(as_hce(KHCE))$Pvalue, 0.000001456397859)
})

test_that("gamma for WP for KHCE, SAS", {
  expect_equal(calcWINS(as_hce(KHCE))$gamma$gamma, 0.142053035505)
})

test_that("SE for gamma for WP for KHCE, SAS", {
  expect_equal(calcWINS(as_hce(KHCE))$SE$gamma_SE, 0.029488695888)
})


