####################################################################
# pairwise sample distributions functions
####################################################################

pairwise_distances <- function(df,
                               sample_ids,
                               result_column_name = "estbinres",
                               id_column_name = "id") {
  # browser()

  sample_lat_long_id <- df[
    which(df[, id_column_name, drop = TRUE] %in% sample_ids),
    c("id", "longitude", "latitude")
  ]
  samp_ids <- sample_lat_long_id[, c("id"), drop = TRUE]
  sample_lat_long <- sample_lat_long_id[, c("longitude", "latitude")]

  all_lat_long <- df[, c("longitude", "latitude")]


  mat <- distm(all_lat_long, sample_lat_long, fun = distGeo) / 1000

  list(
    sample_ids = samp_ids,
    pairwise_distance_matrix = mat
  )
}

####################################################################
# nearest neighbour function
####################################################################


nn_func <- function(results_data,
                    ordered_pairwise_distance_matrix,
                    pairwise_distance_samples,
                    num_nearest_neighbours,
                    sample_position) {
  # browser()
  # subset the ordered pd matrix, +1 because distance to self will be removed
  nn_mat <- ordered_pairwise_distance_matrix[1:(num_nearest_neighbours + 1), ]

  # find sample id corresponding to the matrix column
  sample_id <- pairwise_distance_samples[sample_position]

  # find nn samples
  nn_samp_id <- nn_mat[, sample_position]
  nn_samp_id <- nn_samp_id[!nn_samp_id %in% sample_id]

  dat <- results_data[which(results_data$id %in% nn_samp_id), ]

  sum(dat$estbinres) / length(dat$estbinres)
}



nn_func2 <- function(results_data,
                     ordered_pairwise_distance_matrix,
                     pairwise_distance_samples,
                     num_nearest_neighbours,
                     sample_position) {
  # browser()
  # subset the ordered pd matrix, +1 because distance to self will be removed
  nn_mat <- ordered_pairwise_distance_matrix[1:(num_nearest_neighbours + 1), ]

  # find sample id corresponding to the matrix column
  sample_id <- pairwise_distance_samples[sample_position]

  # find nn samples
  nn_samp_id <- nn_mat[, sample_position]
  nn_samp_id <- nn_samp_id[!nn_samp_id %in% sample_id]

  dat <- results_data[which(results_data$id %in% nn_samp_id), ]

  prev <- overall_prev(dat)

  prev$Lower
}

####################################################################
# nearest neighbour distance function
####################################################################

nd_func <- function(nd_ordered_pairwise_distance_matrix,
                    pairwise_distance_samples,
                    sample_position,
                    results_data) {
  # browser()

  # find sample id corresponding to the matrix column
  sample_id <- pairwise_distance_samples[sample_position]

  # find nn samples
  nn_samp_id <- nd_ordered_pairwise_distance_matrix[, sample_position]
  nn_samp_id <- nn_samp_id[!nn_samp_id %in% sample_id]

  dat <- results_data[which(results_data$id %in% nn_samp_id), ]

  sum(dat$estbinres) / length(dat$estbinres)
}



# ####################################################################
# # spatial analysis functions
# ####################################################################
#
# percent_positive <- function(region, df, region_name_column){
#
#   #subset data by region
#   samples <- df[df[ , region_name_column] == region,]
#   n_samples <- nrow(samples)
#
#   #calculate the percentage of positive results
#   pos <- samples[which(samples$estbinres == 1),]
#   n_pos <- nrow(pos)
#   percent_pos <- n_pos/n_samples * 100
#
#   #calculate 95% confidence intervals
#   standard_error <- sqrt(percent_pos * (100 - percent_pos)/n_samples)
#   upper_ci <- min(percent_pos + 1.96 * standard_error, 100)
#   lower_ci <- max(percent_pos - 1.96 * standard_error, 0)
#
#   # return the percent pos, confidence intervals, and region
#   data.frame(percent_pos = percent_pos,
#              lower_ci = lower_ci,
#              upper_ci = upper_ci,
#              regions = region)
# }
#
#
#
# # plots
# plot_samples_per_regions <- function(df, out_path,
#                                      plot_theme = plot_theme){
#   #distribution of samples per region
#   df$region <- as.factor(df$region)
#   p <- ggplot(df, aes(x = region, fill = region)) +
#     geom_bar() +
#     plot_theme +
#     theme(legend.position = "none") +
#     ggtitle("Samples per Region") +
#     xlab("Region") +
#     ylab("#Samples")
#
#
#     return(p)
#
# }
#
#
# plot_positives_per_regions <- function(df, out_path,
#                                        plot_theme = plot_theme){
#   #distribution of positive samples per region
#   p <- ggplot(df2, aes(x = regions, y = percent_pos, fill = regions)) +
#     geom_bar(stat = "identity") +
#     geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci)) +
#     plot_theme +
#     theme(legend.position = "none") +
#     xlab("Region") +
#     ylab("Percent Positive")
#
#     return(p)
#
# }
#
#
# # map plots
#
# map_all_results <- function(d_sh, england_sf, out_path, plot_theme){
#   warning("Have you jittered latitude and longitude?")
#
#   #map of all results
#   p <- ggplot() +
#     geom_sf(data = england_sf, fill = "grey", colour = "grey") +
#     geom_sf(data = d_sh[d_sh$estbinres == FALSE,], colour = "black",
#             size = 0.2) +
#     geom_sf(data = d_sh[d_sh$estbinres == TRUE,], colour = "red",
#             size = 0.8) +
#     ggtitle(sprintf("n = %s",nrow(d_sh))) +
#     plot_theme
#
#     return(p)
#
# }
#
#
# # daily_maps <- function(d_sh, england_sf, out_path, plot_theme){
# #   warning("Have you jittered latitude and longitude?")
# #
# #   d_sh <- d_sh[!is.na(d_sh$d_comb),]
# #   dates <- as.character(sort(unique(d_sh$d_comb)))
# #   d_sh$d_comb <- as.character(d_sh$d_comb)
# #   for (date in dates){
# #     data <- d_sh[d_sh$d_comb == date,]
# #
# #     n <- nrow(data)
# #
# #     p <- ggplot() +
# #       geom_sf(data = england_sf, fill = "grey", colour = "grey") +
# #       geom_sf(data = data[data$estbinres == FALSE,], colour = "black",
# #               size = 0.2) +
# #       ggtitle(sprintf("%s, n = %s",date, n)) +
# #       plot_theme
# #
# #     #plot points in two layers so detected can be clearly seen
# #   #   if (table(data$estbinres)[[TRUE]] > 0){
# #   #     p <- p +
# #   #       geom_sf(data = data[data$estbinres == TRUE,], colour = "red",
# #   #               size = 0.8)
# #   #   }
# #   #   ggsave(sprintf("%spng_files/%s_results.png", out_path, date))
# #   #   saveRDS(p, sprintf("%srds_files/%s_results_plot.rds", out_path, date))
# #   # }
# #   # message("Plots saved to %srds_files and %spng_files", out_path, out_path)
# # }
#
#
# map_results_to_regions <- function(df, england_regions, out_path, plot_theme){
#   #browser()
#   to_plot <- merge(england_regions, df, by.x = "rgn17nm", by.y = "regions")
#
#   p <- ggplot() +
#     geom_sf(data = to_plot, aes(fill = percent_pos)) +
#     plot_theme +
#     ggtitle("Percent Prevalence per Region") +
#     theme(legend.title = element_text(),
#           legend.position = "right") +
#     labs(fill = "Prevalence (%)") +
#     scale_fill_gradient(low = "#FEF0F0" , high = "#450202")
#   # scale_fill_continuous(trans = "reverse") +
#   # guides(colour = guide_legend(reverse = FALSE))
#
#   return(p)
#
# }
#
#
# ####################################################################
# # pairwise sample distributions functions
# ####################################################################
#
# pairwise_distances <- function(df,
#                                sample_ids,
#                                result_column_name = "res",
#                                id_column_name = "id"){
#   #browser()
#   res <- unique(df[ , result_column_name])
#
#   if(res == 0){
#     covid_result <- "negative"
#   } else if(res ==1){
#     covid_result <- "positive"
#   } else {
#     stop("non-unique test result in sample")
#   }
#
#   lat_long <- df[which(df[ , id_column_name, drop = TRUE] %in% sample_ids), c("longitude", "latitude")]
#   mat <- distm(lat_long, fun = distGeo)/1000
#   distances <- mat[which(lower.tri(mat, diag = FALSE))]
#
#   list(covid_result = covid_result,
#        sample_ids = sample_ids,
#        pairwise_distances = distances,
#        pairwise_distance_matrix = mat)
# }
#
#
#
# list_to_tibble <- function(out_list){
#   #browser()
#
#   tibble(pairwise_distance = out_list$pairwise_distances,
#          covid_result = rep(out_list$covid_result, length = length(out_list$pairwise_distances)))
# }
#
#
#
# sample_id_frequency <- function(sample_results,
#                                  max_distance,
#                                  min_distance){
#
#   #browser()
#
#   mat <- sample_results$pairwise_distance_matrix
#
#   # set main diagonal and upper tri to NA
#   # set values above max dist to NA
#   mat[which(upper.tri(mat, diag = TRUE))] <- NA
#   mat[which(mat > max_distance)] <- NA
#   mat[which(mat < min_distance)] <- NA
#
#   # find the matrix indices for where we have non-NA values
#   values_indices <- which(!(is.na(mat)), arr.ind = TRUE)
#   # make these values a vector
#   values_indices <- as.vector(values_indices)
#
#   sample_ids <- sample_results$sample_ids[values_indices]
#   sample_ids <- tibble(id = sample_ids)
#
#   samp_frequency <- sample_ids %>%
#     group_by(id) %>%
#     dplyr::count() %>%
#     rename(frequency = n)
#
#   # find the samples with 0 frequency
#   zero_samps <- sample_results$sample_ids[which(!(sample_results$sample_ids %in% samp_frequency$id))]
#
#   zero_freq <- data.frame(id = zero_samps,
#                           frequency = rep(0L, times = length(zero_samps)))
#
#   ret <- bind_rows(samp_frequency, zero_freq)
#
#   ret
# }
#
#
# ##################################################################
# # investigating clusters
# ###################################################################
#
# calc_freq_freq <- function(samp_freq_tibble){
#
#   #browser()
#   freq_freq <- samp_freq_tibble %>%
#     rename(id_frequency = frequency) %>%
#     group_by(id_frequency) %>%
#     dplyr::count() %>%
#     rename(frequency_of_id_frequency = n)
#
#   # need to ensure that we have an count value for all integers - even if the
#   #  id_frequency does not appear in our data.  So we need to introduce some 0
#   #  counts for id_frequency values not in our data
#
#   # create a new tibble for this
#   temp <- tibble(id_frequency = seq(0, max(freq_freq$id_frequency)))
#
#   freq_freq <- left_join(temp,
#                          freq_freq,
#                          by = c("id_frequency" = "id_frequency"))
#
#   freq_freq$frequency_of_id_frequency[is.na(freq_freq$frequency_of_id_frequency)] <- 0
#
#   freq_freq$cumu_sum <- cumsum(freq_freq$frequency_of_id_frequency)
#
#   freq_freq
# }
#
#
# find_neg_samps <- function(number_of_points,
#                            high_negs_df){
#
#   neg_rows <- sample(x = 1:nrow(high_negs_df), size = number_of_points, replace = FALSE)
#   high_negs_df[neg_rows, ]
# }
#
# calc_convex_hull <- function(df){
#   #browser()
#   ch_rows <- chull(x = df$latitude,
#                    y = df$longitude)
#
#   ch_points <- df[ch_rows, c("latitude", "longitude")]
#
#   ch_area <- polyarea(x = ch_points$longitude, y = ch_points$latitude)
#
#   list(convex_hull_points = ch_points, convex_hull_area = ch_area)
# }
#
#
#
# ##################################################################
# # old functions
# ##################################################################
#
#
# # dist <- function(ids, d_sh){
# #   #browser()
# #   lat_lons <- d_sh[d_sh$id %in% ids, c("longitude", "latitude")] %>%
# #     as.matrix()
# #   lat_lons_list <- lapply(seq_len(nrow(lat_lons)), function(x) lat_lons[x,])
# #   dist_km <- geosphere::distGeo(lat_lons)[1]/1000
# # }
# #
# #
# #
# # sampler <- function(sample_number,
# #                     data,
# #                     number_test_results,
# #                     result_type,
# #                     result_column_name,
# #                     sample_draw_replacement = FALSE){
# #
# #   #browser()
# #
# #   results_values <- unique(data[, result_column_name, drop = TRUE])
# #
# #   if(any(!(results_values %in% c(0,1)))){
# #
# #     stop("invalid results values. values should be in c(0,1)")
# #   }
# #
# #   if(result_type == "positive"){
# #     res_value <- 1
# #   } else if(result_type == "negative") {
# #     res_value <- 0
# #   } else {
# #     message("incorrect result_type entered: should be 'positive' or 'negative'")
# #   }
# #
# #   # draw a sample from the results
# #   ids <- base::sample(x = data[data[, result_column_name] == res_value, ]$id,
# #                       size = number_test_results,
# #                       replace = sample_draw_replacement)
# #
# #   combos <- t(combn(ids, 2))
# #
# #   distances <- apply(combos, 1, dist, d_sh = data)
# #
# #   data.frame("pairwise_distance" = distances,
# #              "test_result" = rep(result_type,
# #                                  length(distances)),
# #              "sample_number" = paste0("sample ", sample_number))
# #
# #
# # }
#
# #
# # sampler <- function(i, data, n = 5){
# #   #browser()
# #   neg_ids <- base::sample(data[data$estbinres == 0, ]$id, n) #random sample of equal length
# #
# #   neg_combos <- t(combn(neg_ids, 2))
# #
# #   neg_distances <- apply(neg_combos, 1, dist, d_sh = data)
# #
# #   new_dists <- cbind(neg_distances, rep(sprintf("All Tests Sample %s", i),
# #                                         length(neg_distances))) %>%
# #     as.data.frame()
# #   names(new_dists) <- c("Pairwise_Distance", "pos_neg")
# #   return(new_dists)
# # }
# #
# #
# #
# # pos_sampler <- function(data){
# #   browser()
# #   pos_ids <- data[data$estbinres == 1, ]$id
# #
# #   pos_combos <- t(combn(pos_ids, 2))
# #
# #   pos_distances <- apply(pos_combos, 1, dist, d_sh = data)
# #
# #   new_dists <- cbind(pos_distances, rep("Positives",
# #                                         length(pos_distances))) %>%
# #     as.data.frame()
# #   names(new_dists) <- c("Pairwise_Distance", "pos_neg")
# #   return(new_dists)
# # }
