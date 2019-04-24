library(testthat)
library(ProteoQ)

test_check("ProteoQ")

context("The index of TMT sets")
# library(stringr)

test_that("set index is numberic", {
	expect_is(set_index, "is.numberic")
	expect_output(channelInfo(lab_schm, set_index), "List of 3")
  # expect_equal(str_length("ab"), 2)
  # expect_equal(str_length("abc"), 3)
})

test_that("single or multiple PSM files", {
	expect_is(set_index, "is.numberic")
	expect_output(channelInfo(lab_schm, set_index), "List of 3")
})

test_that("filenames in Windows OS versus in Xcalibur", {
	expect_is(set_index, "is.numberic")
	expect_output(channelInfo(lab_schm, set_index), "List of 3")
})

