library(magick)

# 1. Read your saved PNG files from disk
img1 <- image_read("./figures/genus_level_rel_ab.png")
img2 <- image_read("./figures/genus_level_rel_ab_AGS.png")

# 2. Combine the images into an array vector
combined_vector <- c(img1, img2)

# 3. Append them together 
# Use stack = FALSE for horizontal juxtaposition, stack = TRUE for vertical layout
horizontal_merge <- image_append(combined_vector, stack = FALSE)
vertical_merge   <- image_append(combined_vector, stack = TRUE)

# 4. Save the final combined image back to your file path
image_write(horizontal_merge, path = "./figures/abundance_wide.png", format = "png")
