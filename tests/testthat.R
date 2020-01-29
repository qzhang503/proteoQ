library(testthat)
library(ProteoQ)

test_check("ProteoQ")

context("The index of TMT sets")

test_that("set index is numberic", {
	expect_is(set_index, "is.numberic")
	expect_output(channelInfo(label_scheme, set_idx), "List of 3")
})



