################################################################################
## Lake Sampling Station Map (v5)
##
## Purpose: Plot sampling stations at three Westchester County lakes with
##          OSM-derived shoreline polygons, a New York State map (high-res
##          local shapefile) with GPS coordinate grid, and a world map
##          highlighting the US/NY location with dashed connector lines.
##
##
## Requirements: osmdata, sf, ggplot2, cowplot, dplyr, grid,
##               rnaturalearth, rnaturalearthdata
##
## Local data:
##   NY state shoreline: data/NYS_Civil_Boundaries.shp/State_Shoreline.shp
##
## Data sources:
##   Lake shorelines: OpenStreetMap via Overpass API (osmdata).
##     Data (c) OpenStreetMap contributors, ODbL 1.0.
##   World map: Natural Earth via rnaturalearth. Public domain.
##
## NOTE: This script requires internet access for the OSM queries.
##       Run from the project root where the data/ directory is accessible.
################################################################################

# --- Load Libraries ---
library(osmdata)
library(sf)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)

# rnaturalearth: world country boundaries for the context map
# Install if needed: install.packages(c("rnaturalearth", "rnaturalearthdata"))
library(rnaturalearth)
library(rnaturalearthdata)


# ==============================================================================
# 1. Define Sampling Stations
# ==============================================================================

# Station coordinates from lakesites.csv
sites <- data.frame(
  Lake = c("Lake Mohegan", "Teatown Lake", "Lake Rippowam"),
  Latitude  = c(41.3175, 41.2087, 41.2995),
  Longitude = c(-73.8529, -73.8343, -73.56),
  stringsAsFactors = FALSE
)

# Assign lake-specific colors
lake_colors <- c(
  "Lake Mohegan"  = "#3498DB",
  "Teatown Lake"  = "#999999",
  "Lake Rippowam" = "#E69F00"
)

# Convert to sf object (WGS84)
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)


# ==============================================================================
# 2. Query Lake Shoreline Polygons from OpenStreetMap
# ==============================================================================

# Helper function: query OSM for a lake polygon by name within a bounding box
# centered on the sampling point. Returns an sf polygon object or NULL.
get_lake_polygon <- function(lake_name, lon, lat, buffer = 0.02) {
  
  # Define bounding box around the station
  bbox <- c(lon - buffer, lat - buffer, lon + buffer, lat + buffer)
  
  # Query OSM for natural=water features
  lake_data <- opq(bbox = bbox) %>%
    add_osm_feature(key = "natural", value = "water") %>%
    osmdata_sf()
  
  # Lake polygons may appear in osm_polygons or osm_multipolygons.
  # Check both and try to match by name.
  poly <- NULL
  
  # Check osm_multipolygons first (larger features often stored here)
  if (!is.null(lake_data$osm_multipolygons) && nrow(lake_data$osm_multipolygons) > 0) {
    mp <- lake_data$osm_multipolygons
    if ("name" %in% names(mp)) {
      match_idx <- grep(gsub("Lake ", "", lake_name), mp$name, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        poly <- mp[match_idx[1], ]
      }
    }
  }
  
  # If not found in multipolygons, check osm_polygons
  if (is.null(poly) && !is.null(lake_data$osm_polygons) && nrow(lake_data$osm_polygons) > 0) {
    sp_data <- lake_data$osm_polygons
    if ("name" %in% names(sp_data)) {
      match_idx <- grep(gsub("Lake ", "", lake_name), sp_data$name, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        poly <- sp_data[match_idx[1], ]
      }
    }
  }
  
  # Fallback: use all water polygons in the bbox if no name match
  if (is.null(poly)) {
    message(paste0("No exact name match for '", lake_name,
                   "'. Using all water polygons within bbox as fallback."))
    all_polys <- NULL
    if (!is.null(lake_data$osm_polygons) && nrow(lake_data$osm_polygons) > 0) {
      all_polys <- lake_data$osm_polygons
    }
    if (!is.null(lake_data$osm_multipolygons) && nrow(lake_data$osm_multipolygons) > 0) {
      if (is.null(all_polys)) {
        all_polys <- lake_data$osm_multipolygons
      } else {
        all_polys <- st_as_sf(
          rbind(st_geometry(all_polys), st_geometry(lake_data$osm_multipolygons))
        )
      }
    }
    poly <- all_polys
  }
  
  return(poly)
}

# Query each lake
message("Querying OSM for Lake Mohegan...")
poly_mohegan <- get_lake_polygon("Lake Mohegan", -73.8529, 41.3175, buffer = 0.015)

message("Querying OSM for Teatown Lake...")
poly_teatown <- get_lake_polygon("Teatown Lake", -73.8343, 41.2087, buffer = 0.015)

message("Querying OSM for Lake Rippowam...")
poly_rippowam <- get_lake_polygon("Lake Rippowam", -73.56, 41.2995, buffer = 0.02)


# ==============================================================================
# 3. Create Individual Lake Panels
# ==============================================================================

# All three lake panels use the same axis text angle (45 degrees) and identical
# theme/margin settings. This is critical for plot_grid alignment — panels with
# different margin structures will not align properly even with align = "vh".
make_lake_panel <- function(lake_polygon, lake_name, lon, lat, color) {
  
  # Create station point sf
  pt <- st_as_sf(
    data.frame(x = lon, y = lat),
    coords = c("x", "y"), crs = 4326
  )
  
  # Build the plot
  p <- ggplot()
  
  # Add shoreline polygon if available
  if (!is.null(lake_polygon)) {
    p <- p +
      geom_sf(data = lake_polygon, fill = "white", color = "grey30",
              linewidth = 0.4)
  }
  
  # Add station point and formatting
  p <- p +
    geom_sf(data = pt, color = color, size = 3, shape = 16) +
    labs(title = lake_name) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      # All panels use 45-degree x-axis rotation for consistent height
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title = element_blank(),
      # Explicit uniform margins: top, right, bottom, left
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    )
  
  return(p)
}

# Create the three panels (all with identical theme structure)
panel_mohegan <- make_lake_panel(
  poly_mohegan, "Lake Mohegan", -73.8529, 41.3175,
  lake_colors["Lake Mohegan"]
)

panel_teatown <- make_lake_panel(
  poly_teatown, "Teatown Lake", -73.8343, 41.2087,
  lake_colors["Teatown Lake"]
)

panel_rippowam <- make_lake_panel(
  poly_rippowam, "Lake Rippowam", -73.56, 41.2995,
  lake_colors["Lake Rippowam"]
)


# ==============================================================================
# 4. Create New York State Map (High-Res Local Shapefile, with Label and Grid)
# ==============================================================================

# Read the high-resolution NYS shoreline shapefile
message("Reading NY state shoreline shapefile...")
# NOT YET RESOLVED: this shapefile (and its accompanying .dbf/.shx/.prj
# siblings) is not present anywhere in data/ in this repository. Add the
# NYS Civil Boundaries shapefile set to data/NYS_Civil_Boundaries.shp/
# before this script can run.
if (!file.exists("data/NYS_Civil_Boundaries.shp/State_Shoreline.shp")) {
  stop("Figure1.R: missing 'data/NYS_Civil_Boundaries.shp/State_Shoreline.shp'. ",
       "Add the shapefile (and its .dbf/.shx/.prj siblings) to data/ before ",
       "running this script.")
}
ny_state <- st_read("data/NYS_Civil_Boundaries.shp/State_Shoreline.shp")

# Ensure CRS is WGS84 (EPSG:4326) for consistent plotting with station coords
ny_state <- st_transform(ny_state, crs = 4326)

# Build the NY state map with GPS coordinate grid and station dots
ny_map <- ggplot() +
  # NY state outline (high-resolution local shapefile)
  geom_sf(data = ny_state, fill = "white", color = "black", linewidth = 0.5) +
  
  # Sampling station points
  geom_point(
    data = sites,
    aes(x = Longitude, y = Latitude, color = Lake),
    size = 3, shape = 16
  ) +
  scale_color_manual(values = lake_colors) +
  
  # GPS coordinate grid via coord_sf
  coord_sf(datum = st_crs(4326)) +
  scale_x_continuous(breaks = seq(-80, -72, by = 2)) +
  scale_y_continuous(breaks = seq(40, 46, by = 1)) +
  
  # Theme with visible coordinate labels and gridlines
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text  = element_text(size = 7, color = "grey30"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3,
                                    linetype = "dotted"),
    panel.grid.minor = element_blank(),
    # Explicit uniform margins matching the bottom-row partner
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )


# ==============================================================================
# 5. Create World Map with US/NY Highlighted
# ==============================================================================

# Load medium-resolution world country boundaries from Natural Earth
message("Loading world map from Natural Earth...")
world <- ne_countries(scale = "medium", returnclass = "sf")

# Separate US from other countries for differential styling
us <- world %>% filter(admin == "United States of America")
other_countries <- world %>% filter(admin != "United States of America")

# Approximate centroid of NY for the marker on the world map
ny_lon <- -75.5
ny_lat <- 43.0

# Build the world map
world_map <- ggplot() +
  # All other countries in light grey
  geom_sf(data = other_countries, fill = "white", color = "grey70",
          linewidth = 0.1) +
  
  # Highlight the US
  geom_sf(data = us, fill = "white", color = "grey50", linewidth = 0.2) +
  
  # Mark NY location with a red point
  annotate("point", x = ny_lon, y = ny_lat,
           color = "red", size = 2, shape = 16) +
  
  # Standard projection, show all countries
  coord_sf(datum = st_crs(4326)) +
  
  # Theme matching the NY state map for visual consistency
  theme_bw(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 7, color = "grey30"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.2,
                                    linetype = "dotted"),
    panel.grid.minor = element_blank(),
    # Explicit uniform margins matching the bottom-row partner
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )


# ==============================================================================
# 6. Build Aligned Rows with plot_grid, Then Compose with ggdraw
# ==============================================================================

# --- Step 1: Align the top row (three lake panels) ---
# plot_grid with align = "vh" and axis = "tblr" forces all panels to share
# the same plot-area dimensions by aligning top, bottom, left, and right
# axes. This ensures identical rendered heights regardless of label content.
top_row <- plot_grid(
  panel_mohegan, panel_rippowam, panel_teatown,
  nrow = 1,
  align = "vh",       # Align both vertically and horizontally
  axis  = "tblr",     # Align on all four sides (top, bottom, left, right)
  rel_widths = c(1, 1, 1)  # Equal width for each panel
)

# --- Step 2: Align the bottom row (NY state map + world map) ---
# Same alignment approach ensures the two bottom panels share plot-area height.
bottom_row <- plot_grid(
  ny_map, world_map,
  nrow = 1,
  align = "vh",
  axis  = "tblr",
  rel_widths = c(1, 1.1)  # World map slightly wider (it spans more longitude)
)

# --- Step 3: Stack the two rows with equal height ---
# plot_grid again to stack the aligned rows vertically.
# rel_heights = c(1, 1) ensures equal row height.
stacked_rows <- plot_grid(
  top_row, bottom_row,
  ncol = 1,
  rel_heights = c(1, 1)  # Equal height for both rows
)

# --- Step 4: Overlay dashed connector lines using ggdraw ---
# ggdraw places the stacked plot as the base, then draws lines on top.
# Coordinates are normalized (0-1) relative to the full canvas.
#
# With equal-height rows:
#   Top row occupies approximately y = [0.50, 1.0]
#   Bottom row occupies approximately y = [0.0, 0.50]
#
# Within the bottom row:
#   NY state map is on the left (~x = [0.0, 0.47])
#   World map is on the right (~x = [0.47, 1.0])
#
# Within the top row:
#   Lake Mohegan is top-left (~x = [0.0, 0.33])
#   Lake Rippowam is top-center (~x = [0.33, 0.67])
#   Teatown Lake is top-right (~x = [0.67, 1.0])

final_plot <- ggdraw(stacked_rows) +
  
  # --- Dashed lines: World map -> NY state map ---
  # From NY marker on world map (approx. x=0.62, y=0.28) to
  # the right edge of the NY state map panel (x=0.47).
  # Two lines bracket the NY state panel vertically.
  draw_line(
    x = c(0.60, 0.47), y = c(0.28, 0.38),
    color = "black", linetype = "dashed", linewidth = 0.5
  ) +
  draw_line(
    x = c(0.60, 0.47), y = c(0.28, 0.18),
    color = "black", linetype = "dashed", linewidth = 0.5
  ) +
  
  # --- Dashed lines: NY state map station dots -> lake panels above ---
  # The three lakes cluster in the SE corner of NY on the state map.
  # The exact x/y of each dot depends on the rendered state geometry.
  # These approximate positions connect to the horizontal center of each
  # lake panel at the bottom edge of the top row (y ~ 0.50).
  #
  # Lake Mohegan: SE corner of NY map -> center of top-left panel
  draw_line(
    x = c(0.38, 0.17), y = c(0.27, 0.50),
    color = lake_colors["Lake Mohegan"],
    linetype = "dashed", linewidth = 0.6
  ) +
  
  # Lake Rippowam: SE corner of NY map -> center of top-center panel
  draw_line(
    x = c(0.39, 0.50), y = c(0.26, 0.50),
    color = lake_colors["Lake Rippowam"],
    linetype = "dashed", linewidth = 0.6
  ) +
  
  # Teatown Lake: SE corner of NY map -> center of top-right panel
  draw_line(
    x = c(0.37, 0.83), y = c(0.25, 0.50),
    color = lake_colors["Teatown Lake"],
    linetype = "dashed", linewidth = 0.6
  ) +
  
  # Attribution label at the very bottom
  draw_label(
    "NY boundary: NYS Civil Boundaries | Lake shorelines: OSM (ODbL 1.0) | World map: Natural Earth",
    x = 0.50, y = 0.005, hjust = 0.5, vjust = 0,
    size = 7, color = "grey50"
  )

# Display
print(final_plot)


# ==============================================================================
# 7. Save the Final Figure
# ==============================================================================

dir.create("figs", showWarnings = FALSE, recursive = TRUE)
ggsave("figs/Figure1.png", final_plot,
       width = 14, height = 10, dpi = 300, bg = "white")

message("Map saved to: figs/Figure1.png")

# Optional: save as PDF for publication
# ggsave("figs/lake_sampling_stations_map.pdf", final_plot,
#        width = 14, height = 10, bg = "white")


################################################################################
## NOTES ON FINE-TUNING
##
## 1. Plot was saved as an SVG file and manually modified so that plots aligned
## 2. New connectors were drawn
## 3. A state label was added
## 4. A black trim was added to the points.
################################################################################