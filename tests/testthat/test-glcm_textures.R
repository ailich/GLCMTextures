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
