
library(tidyverse)
library(stars)
library(furrr)

options(future.fork.enable = T)
plan(multicore)


# Function to get cell position index from coords
fn_get_cell_pos <- 
  
  function(star_obj, bbox) {
    
    index_x <-
      bbox[c(1,3)] %>% 
      map(~ .x - st_get_dimension_values(star_obj, 1)) %>% 
      map(abs) %>% 
      map_int(which.min)
    
    index_y <- 
      bbox[c(2,4)] %>% 
      map(~ .x - st_get_dimension_values(star_obj, 2)) %>% 
      map(abs) %>% 
      map_int(which.min)
    
    c(index_x[1], index_y[1],
      index_x[2], index_y[2])
  }





# Load cmip

cmip <- 
  "/mnt/pers_disk/daily_cmip/" %>% 
  fs::dir_ls() %>%
  .[c(2,5)] %>%                                                                 # *************
  future_map(read_ncdf, ncsub = cbind(start = c(133, 43, 1),
                                      count = c(11, 14, NA)))


"/mnt/pers_disk/daily_cmip/" %>% 
  fs::dir_ls() %>%
  .[c(2,5)] %>% 
  map(read_ncdf, proxy = T) %>% 
  map(st_get_dimension_values, 3) -> time_vector

time_vector %>% do.call(c, .) %>% unname() -> time_vector


cmip <- 
  do.call(c, cmip)



# What cells should be processed?
# not all cells overlap with w5e5 cells with data

foo <- 
  "/mnt/pers_disk/nz_south_chelsa-w5/" %>% 
  fs::dir_ls(regexp = ".nc") %>% 
  .[1] %>% 
  read_ncdf(ncsub = cbind(start = c(1, 1, 1),
                          count = c(NA, NA, 1)))

bar <- 
  foo %>%
  adrop() %>%
  st_warp(cmip %>% slice(time, 1), 
          method = "max", 
          use_gdal = T, 
          no_data_value = 1)

cmip[is.na(bar)] <- NA

cmip <- 
  cmip %>% 
  mutate(pr = pr %>% units::set_units(kg/m^2/d)) %>% 
  units::drop_units()


# Table of cells to process

cells_tb <- 
  cmip %>% 
  slice(time, 1) %>% 
  as_tibble(center = F) %>% 
  filter(!is.na(pr)) %>% 
  select(lon, lat)


delta_lon <- st_dimensions(cmip)$lon$delta
delta_lat <- st_dimensions(cmip)$lat$delta

length_time_cmip <- unname(dim(cmip)[3])


# TILE (CELL) LOOP ************************************************************

for (i in 22:nrow(cells_tb)) {
  # for (i in 2:4) {
  
  print(str_glue("PROCESSING CELL {i}"))
  
  plan(multicore)
  
  lon <- cells_tb$lon[i]
  lat <- cells_tb$lat[i]
  
  # Define extent to process
  # 3 x 3
  
  ext_to_process <- 
    c(st_point(c(lon-delta_lon, lat-delta_lat)),
      st_point(c(lon+2*delta_lon, lat+2*delta_lat))) %>% 
    st_sfc(crs = 4326) %>% 
    st_bbox()
  
  ext_to_process_index <- 
    fn_get_cell_pos(foo, ext_to_process)
  
  
  # Load w5e5
  # only ext to process
  
  print(str_glue("   Importing w5e5"))
  
  w5e5 <- 
    "/mnt/pers_disk/nz_south_chelsa-w5/" %>% 
    fs::dir_ls(regexp = ".nc") %>% 
    future_map(read_ncdf, proxy = F, ncsub = cbind(start = c(ext_to_process_index[1],
                                                             ext_to_process_index[2],
                                                             1),
                                                   count = c(diff(ext_to_process_index[c(1,3)]),
                                                             diff(ext_to_process_index[c(2,4)]),
                                                             NA))) %>% 
    suppressMessages() 
  
  w5e5 <-   
    do.call(c, w5e5)
  
  length_time_w5e5 <- unname(dim(w5e5)[3])
  
  # align to cmip
  w5e5 <- 
    w5e5 %>% 
    st_warp(st_as_stars(ext_to_process, dx = 0.025))
  
  # convert to mm
  w5e5 <- 
    w5e5 %>% 
    mutate(pr = pr %>% units::set_units(kg/m^2/d)) %>% 
    units::drop_units()
  
  # aggregate to coarse res
  
  print(str_glue("   Aggregating w5e5 spatially"))
  
  w5e5_coarse <-
    w5e5 %>%
    st_warp(st_as_stars(ext_to_process, dx = 1.25, dy = 1), 
            method = "average",
            use_gdal = T,
            no_data_value = -1) %>% 
    setNames("pr")
  
  
  # prepare predictors
  
  print(str_glue("   Prepare predictors (w5e5)"))
  
  d <- vector("list", 10)
  
  d[[10]] <- w5e5_coarse %>% 
    setNames("center")
  
  d[[2]] <- d[[10]]
  st_dimensions(d[[2]])$y$offset <- st_dimensions(d[[10]])$y$offset - 1
  d[[2]] <- setNames(d[[2]], "north")
  
  d[[3]] <- d[[10]]
  st_dimensions(d[[3]])$y$offset <- st_dimensions(d[[10]])$y$offset + 1
  d[[3]] <- setNames(d[[3]], "south")
  
  d[[4]] <- d[[10]]
  st_dimensions(d[[4]])$x$offset <- st_dimensions(d[[10]])$x$offset - 1
  d[[4]] <- setNames(d[[4]], "east")
  
  d[[5]] <- d[[10]]
  st_dimensions(d[[5]])$x$offset <- st_dimensions(d[[10]])$x$offset + 1
  d[[5]] <- setNames(d[[5]], "west")
  
  d[[6]] <- d[[10]]
  st_dimensions(d[[6]])$y$offset <- st_dimensions(d[[10]])$y$offset - 1
  st_dimensions(d[[6]])$x$offset <- st_dimensions(d[[10]])$x$offset - 1
  d[[6]] <- setNames(d[[6]], "northeast")
  
  d[[7]] <- d[[10]]
  st_dimensions(d[[7]])$y$offset <- st_dimensions(d[[10]])$y$offset - 1
  st_dimensions(d[[7]])$x$offset <- st_dimensions(d[[10]])$x$offset + 1
  d[[7]] <- setNames(d[[7]], "northwest")
  
  d[[8]] <- d[[10]]
  st_dimensions(d[[8]])$y$offset <- st_dimensions(d[[10]])$y$offset + 1
  st_dimensions(d[[8]])$x$offset <- st_dimensions(d[[10]])$x$offset - 1
  d[[8]] <- setNames(d[[8]], "southeast")
  
  d[[9]] <- d[[10]]
  st_dimensions(d[[9]])$y$offset <- st_dimensions(d[[10]])$y$offset + 1
  st_dimensions(d[[9]])$x$offset <- st_dimensions(d[[10]])$x$offset + 1
  d[[9]] <- setNames(d[[9]], "southwest")
  
  d[2:9] <-
    map(d[2:9], st_warp, d[[10]]) 
  
  
  # crop to tile (cell)
  
  ext_to_train <- 
    c(st_point(c(lon, lat)),
      st_point(c(lon+delta_lon, lat+delta_lat))) %>% 
    st_sfc(crs = 4326) %>% 
    st_bbox()
  
  d[2:10] <- 
    map(d[2:10], function(s) {
      
      s <- s[ext_to_train]
      return(s)
      
    })
  
  # add w5e5 (response)
  d[[1]] <- w5e5[ext_to_train]
  
  # regrid to fine res
  
  print(str_glue("   Regridding predictors (w5e5)"))
  
  d[2:10] <- 
    map(d[2:10], st_warp, d[[1]])
  
  # homogeneize dimensions
  d[2:10] <-
    map(d[2:10], function(s) {
      
      st_dimensions(s) <- st_dimensions(d[[1]])
      return(s)
      
    })
  
  d <-
    do.call(c, c(d, along = "v"))
  
  gc()
  
  
  # now the same with cmip *****************************************************
  
  cmip_sub <- 
    cmip[ext_to_process]
  
  cmip_sub <- 
    cmip_sub %>% 
    st_set_dimensions(1, names = "x") %>% 
    st_set_dimensions(2, names = "y")
  
  # prepare predictors
  
  print(str_glue("   Prepare predictors (cmip)"))
  
  e <- vector("list", 10)
  
  e[[10]] <- cmip_sub %>% 
    setNames("center")
  
  e[[2]] <- e[[10]]
  st_dimensions(e[[2]])$y$offset <- st_dimensions(e[[10]])$y$offset - 1
  e[[2]] <- setNames(e[[2]], "north")
  
  e[[3]] <- e[[10]]
  st_dimensions(e[[3]])$y$offset <- st_dimensions(e[[10]])$y$offset + 1
  e[[3]] <- setNames(e[[3]], "south")
  
  e[[4]] <- e[[10]]
  st_dimensions(e[[4]])$x$offset <- st_dimensions(e[[10]])$x$offset - 1
  e[[4]] <- setNames(e[[4]], "east")
  
  e[[5]] <- e[[10]]
  st_dimensions(e[[5]])$x$offset <- st_dimensions(e[[10]])$x$offset + 1
  e[[5]] <- setNames(e[[5]], "west")
  
  e[[6]] <- e[[10]]
  st_dimensions(e[[6]])$y$offset <- st_dimensions(e[[10]])$y$offset - 1
  st_dimensions(e[[6]])$x$offset <- st_dimensions(e[[10]])$x$offset - 1
  e[[6]] <- setNames(e[[6]], "northeast")
  
  e[[7]] <- e[[10]]
  st_dimensions(e[[7]])$y$offset <- st_dimensions(e[[10]])$y$offset - 1
  st_dimensions(e[[7]])$x$offset <- st_dimensions(e[[10]])$x$offset + 1
  e[[7]] <- setNames(e[[7]], "northwest")
  
  e[[8]] <- e[[10]]
  st_dimensions(e[[8]])$y$offset <- st_dimensions(e[[10]])$y$offset + 1
  st_dimensions(e[[8]])$x$offset <- st_dimensions(e[[10]])$x$offset - 1
  e[[8]] <- setNames(e[[8]], "southeast")
  
  e[[9]] <- e[[10]]
  st_dimensions(e[[9]])$y$offset <- st_dimensions(e[[10]])$y$offset + 1
  st_dimensions(e[[9]])$x$offset <- st_dimensions(e[[10]])$x$offset + 1
  e[[9]] <- setNames(e[[9]], "southwest")
  
  e[2:9] <-
    map(e[2:9], st_warp, e[[10]])
  
  # crop to tile (cell)
  
  e[2:10] <- 
    map(e[2:10], function(s) {
      
      s <- s[ext_to_train]
      return(s)
      
    })
  
  # regrid to fine res
  
  print(str_glue("   Regridding predictors (cmip)"))
  
  e[2:10] <- 
    map(e[2:10], st_warp, st_as_stars(ext_to_train, dx = 0.025))
  
  # add first v (necessary for merging later)
  e[[1]] <- e[[2]]
  
  # # homogeneize dimensions
  # d[2:10] <-
  #   map(d[2:10], function(s) {
  #     
  #     st_dimensions(s) <- st_dimensions(d[[1]])
  #     return(s)
  #     
  #   })
  
  e <-
    do.call(c, c(e, along = "v"))
  
  gc()
  
  
  # merge both 
  
  print(str_glue("   Merging w5e5 and cmip"))
  
  e <- c(d, e, along = "time")
  
  gc()
  
  
  # Downscale
  
  print(str_glue("   Downscaling"))
  
  p <- 
    e %>% 
    st_apply(c(1,2), function(x) {
      
      narm <- which(is.na(x[1,]))
      
      if (1 %in% narm) {
        
        p <- rep(NA_real_, length_time_cmip)
        
      } else if (length(narm) > (10-2)) {
        
        p <- rep(NA_real_, length_time_cmip)
        
      } else {
        
        if (length(narm) > 0) {
          x <- as.data.frame(x[,-narm])
        } else {
          x <- as.data.frame(x)
        }
        
        # m <- lm(V1 ~ ., data = x[1:length_time_w5e5,])
        # p <-
        #   predict(m, x[(length_time_w5e5+1):nrow(x), -1]) %>%
        #   unname() %>%
        #   {if_else(. < 0, 0, .)}
        
        m <- randomForest::randomForest(V1 ~ ., data = x[1:length_time_w5e5,], ntree = 30)
        p <-
          predict(m, x[(length_time_w5e5+1):nrow(x), -1]) %>%
          unname()
        
        
        # cut(x$V1, breaks = 10, labels = F)
        # 
        # tictoc::tic()
        # m <- Cubist::cubist(x = x[,-1], y = x[,1])
        # tictoc::toc()
        # 
        # p <- 
        #   predict(m, x[(length_time_w5e5+1):nrow(x), -1]) %>%
        #   unname()
        
      }
      
      gc()
      
      return(p)
      
    },
    .fname = "time",
    FUTURE = T)
  
  
  # save tile
  
  print(str_glue("   Saving"))
  
  p %>% 
    aperm(c(2,3,1)) %>% 
    write_stars(str_glue("/mnt/pers_disk/temp/cell_{str_pad(i, 2, 'left', '0')}_rf.tif"))
  
  
  rm(p, e, d, w5e5, w5e5_coarse)
  gc()
  plan(sequential)

  
  
  tempdir() %>% 
    fs::dir_ls() %>% 
    walk(file.remove)
  
  
}


"/mnt/pers_disk/temp/" %>% 
  fs::dir_ls(regexp = "rf", invert = T) %>% 
  map(read_stars, RasterIO = list(bands = 6999)) -> a



"/mnt/pers_disk/temp/" %>% 
  fs::dir_ls(regexp = "rf", invert = T) %>% 
  map(read_stars, RasterIO = list(bands = 1:4015)) -> a

a %>% 
  map(st_apply, c(1,2), mean) -> a

do.call(st_mosaic, a) -> aa

aa %>% mapview::mapview()



"/mnt/pers_disk/temp/" %>% 
  fs::dir_ls(regexp = "rf", invert = T) %>% 
  map(read_stars, RasterIO = list(bands = 9126:12775)) -> a

a %>% 
  map(st_apply, c(1,2), mean) -> a

do.call(st_mosaic, a) -> aa

aa %>% mapview::mapview()





cmip %>% 
  filter(year(time) >= 1990,
         year(time) < 2000 ) %>% 
  st_apply(c(1,2), mean) %>% 
  mapview::mapview()
  


tempdir() %>% 
  fs::dir_ls() %>% 
  walk(unlink, recursive = T)

















for (lon in 50:60) {
  for (lat in 20:30) {
    
    print(str_glue("{lon} - {lat}"))
    
    x <- 
      d[,lon,lat,,] %>% 
      pull() %>% 
      adrop(c(1,2))
    
    narm <- which(is.na(x[1,]))
    
    if (1 %in% narm) {
      
      p <- rep(NA_real_, dim(x)[1])
      
    } else {
      
      x <- as.data.frame(x[,-narm])
      m <- lm(V1 ~ ., data = x)
      p <- predict(m) %>% unname() %>% {if_else(. < 0, 0, .)}
      
    }
    
    print(str_glue("{lon} - {lat} - {narm}"))
    
  }
  
}








e[,25,25,,] %>% pull() %>% adrop(c(1,2)) -> x

e[,2,20,,] %>% pull() %>% adrop(c(1,2)) -> x

narm <- which(is.na(x[1,]))

as.data.frame(x[,-narm]) -> x

m <- lm(V1 ~ ., data = x)

predict(m) %>% unname() %>% min()

broom::augment(m) %>% 
  ggplot(aes(V1, .fitted)) + 
  geom_point()

m <- lm(pr ~ ., data = tb %>% select(-x,-y,-time))






ggplot(tb, aes(pr, center)) + geom_point()
















