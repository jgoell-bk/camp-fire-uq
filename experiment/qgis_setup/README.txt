Simulation-setup layers for QGIS (EPSG:26910, NAD83 / UTM zone 10N)
Written by tools/qgis_setup_export.py on 2026-10-20.
Preview: ../figures/setup_preview.png

setup.gpkg layers
  domain                ELMFIRE domain, 1457 x 1476 cells of 30 m (about 44 x 44 km)
  ignition              Fixed ignition cell, 39.810283 N 121.437222 W, 06:30 PST 8 Nov 2018
  stations              Jarbo Gap (041214), Openshaw (041215), Cohasset (041211) RAWS.
                        Fields: elev_ft, role, mean_ws_mph, mean_wind_from_deg (vector
                        mean of the observed wind over the 24 h window, 06:00 PST 8 Nov
                        to 06:00 PST 9 Nov). All three are inside the domain.
  paradise_town_2018    Paradise town limits, US Census TIGER/Line 2018 (GEOID 0655520, 47.5 km2)
  magalia_cdp_2018      Magalia census-designated place, same source (GEOID 0645120), context
  paradise_town_cells   Outline of the model cells whose centres are inside the town limits
  paradise_border_cells Wildland border band: burnable, non-urban model cells whose centres are
                        within 1 km outside the town limits (33,758 cells, 65 polygons; gaps
                        are building or nonburnable cells). Preview: ../figures/paradise_border_preview.png
  paradise_2km_standin  The 2 km circle used as a stand-in in the earlier pilot summaries

Rasters (same 30 m grid as the model)
  paradise_town_mask.tif  1 = cell centre inside Paradise town limits (52,734 cells, 11,392
                          of them with buildings). Usable as the Paradise mask:
                            python -m sobolfire assemble JG --mask qgis_setup/paradise_town_mask.tif
  paradise_border_mask.tif  1 = cell in the wildland border band above (same cells)
  building_cells.tif      1 = cell with buildings in the case's inputs (urban fuel model 91)
  hillshade.tif           Shaded relief from dem.tif (sun from the north-west, 45 degrees up)
  dem.tif                 The case's elevation raster (metres)

Sources
  Station coordinates: NIFC "Public View - Interagency Remote Automatic Weather Stations
    (RAWS)", https://data-nifc.opendata.arcgis.com/datasets/public-view-interagency-remote-automatic-weather-stations-raws
    retrieved 2026-10-20. These are current coordinates; DRI/WRCC and Iowa Mesonet list
    the same locations to within about 15 m.
      Jarbo Gap 39.73591, -121.48898, 2535 ft ("on a ridge above the Feather River canyon")
      Openshaw  39.58983, -121.63516,  268 ft
      Cohasset  39.87184, -121.76898, 1733 ft
  Town boundaries: US Census Bureau TIGER/Line 2018, California places,
    https://www2.census.gov/geo/tiger/TIGER2018/PLACE/tl_2018_06_place.zip
    (sha256 e95c233451b461fb246966edb935553ecb505105036a9bc24170a7f5d6b8d1b4),
    kept in ../external/tl_2018_06_place/
  Domain, ignition, buildings, elevation: the ELMFIRE case inputs (../config.yaml).
