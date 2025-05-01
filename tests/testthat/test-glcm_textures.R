m<- matrix(c(0,0,0,2,0,0,2,2,1,1,2,3,1,1,2,3), nrow=4) # test matrix from Hall-Beyer tutorial
h_glcm<- matrix(c(4,2,1,0,2,4,0,0,1,0,6,1,0,0,1,2), nrow = 4) # horizontal GLCM from Hall-Beyer tutorial
h_glcm_norm<- matrix(c(0.166, 0.083, 0.042,0,0.083,0.166,0,0,0.042,0,0.25,0.042,0,0,0.042,0.083), nrow=4) # normalized horizontal GLCM from Hall-Beyer tutorial
v_glcm<- matrix(c(6,0,2,0,0,4,2,0,2,2,2,2,0,0,2,0), nrow=4) # vertical GLCM from Hall-Beyer tutorial
v_glcm_norm<- matrix(c(0.249, 0, 0.083, 0, 0, 0.166,0.083,0, 0.083,0.083,0.083,0.083, 0,0,0.083,0), nrow=4) # normalized vertical GLCM from Hall-Beyer tutorial

# Ansers to exercises in Hall-Beyer Guide
h_metrics<- c(glcm_contrast= 0.586, glcm_dissimilarity= 0.418, glcm_homogeneity = 0.807,
              glcm_ASM = 0.145, glcm_entropy = 2.095, glcm_mean= 1.292, glcm_variance= 1.039,
              glcm_correlation = 0.718)
# In Hall-Beyer correlation is 0.691 ibut they appear to have an error in denominator where forgot to take square root even though it is shown.
# [(1.039067)(1.039067)]^0.5 is 1.039067 (the variance) not 1.07966 which is 1.039067^2 but they appear to have an error in denominator where forgot to take square root even though it is shown.
#[(1.039067)(1.039067)]^0.5 is 1.039067 (the variance) not 1.07966 which is 1.039067^2
v_metrics<- c(glcm_contrast= 0.996, glcm_dissimilarity = 0.664, glcm_mean = 1.162, glcm_variance = 0.968) # Hall-Beyer says var=0.970 but they accidentally put 0.50 instead of 0.249 in v_glcm

test_that("Test that horizontal GLCM counts works", {
  expect_equal(h_glcm, make_glcm(m, n_levels = 4, shift = c(1,0), normalize = FALSE))
})

test_that("Test that horizontal GLCM norm works", {
  expect_equal(sum(h_glcm_norm - round(make_glcm(m, n_levels = 4, shift = c(1,0), normalize = TRUE),3)), -0.002) #rounding issue so will be off in two spots by 0.001
})

test_that("Test that vertical GLCM counts works", {
  expect_equal(v_glcm, make_glcm(m, n_levels = 4, shift = c(0,1), normalize = FALSE))
})

test_that("Test that vertical GLCM norm works", {
  expect_equal(sum(v_glcm_norm - round(make_glcm(m, n_levels = 4, shift = c(0,1), normalize = TRUE),3)), -0.002) #rounding issue so will be off in two spots by 0.001
})

test_that("Test glcm_metrics horizontal", {
  expect_equal(round(glcm_metrics(h_glcm_norm, metrics = names(h_metrics)),3), h_metrics) #rounding issue so will be off in two spots by 0.001
})

test_that("Test glcm_metrics vertical", {
  expect_equal(round(glcm_metrics(v_glcm_norm, metrics = names(v_metrics)),3), v_metrics) #rounding issue so will be off in two spots by 0.001
})

test_that("glcm_textures ep works", {
  txt_ep_extpected <- readRDS(system.file("testdata", "txt_ep.RDS", package = "GLCMTextures"))
  r<- rast(volcano, extent= ext(2667400, 2667400 + ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), crs = "EPSG:27200") #Use preloaded volcano dataset as a raster
  txt_ep<- glcm_textures(r, w = c(3,5), n_levels = 16, quant_method = "prob", shift=list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)),
                         metrics = c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"))
  txt_ep<- values(txt_ep, mat=TRUE)
  expect_equal(txt_ep, txt_ep_extpected)
})

test_that("glcm_textures er works", {
  txt_er_extpected <- readRDS(system.file("testdata", "txt_er.RDS", package = "GLCMTextures"))
  r<- rast(volcano, extent= ext(2667400, 2667400 + ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), crs = "EPSG:27200") #Use preloaded volcano dataset as a raster
  txt_er<- glcm_textures(r, w = c(3,5), n_levels = 16, quant_method = "range", shift=list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)),
                         metrics = c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"))
  txt_er<- values(txt_er, mat=TRUE)
  expect_equal(txt_er, txt_er_extpected)
})

test_that("glcm_textures NA handling and ordering of metrics works", {
  txt_ep32_NA_expected <- readRDS(system.file("testdata", "txt_ep32_NA.RDS", package = "GLCMTextures"))
  r<- rast(volcano, extent= ext(2667400, 2667400 + ncol(volcano)*10, 6478700, 6478700 + nrow(volcano)*10), crs = "EPSG:27200") #Use preloaded volcano dataset as a raster
  set.seed(5)
  r[sample(x = 1:ncell(r), size = 100)]<- NA
  txt_ep32_NA<- glcm_textures(r, w = c(3,5), n_levels = 32, quant_method = "prob", shift=list(c(1, 0), c(1, 1), c(0, 1), c(-1, 1)),
                              metrics = rev(c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation")),
                              na.rm=TRUE, impute_corr = TRUE)
  txt_ep32_NA<- values(txt_ep32_NA, mat=TRUE)
  expect_equal(txt_ep32_NA, txt_ep32_NA_expected)
})
