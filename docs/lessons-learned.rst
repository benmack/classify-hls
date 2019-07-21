Lessons Learned
^^^^^^^^^^^^^^^

... or: **Things I would probably make different...**

* Do not keep the tile identifier at the beginning or the filename of a single feature, e.g.

  *data/processed/L3/raster/32UNU/scoll01/32unu__scoll01__vts4w__2018-01-07__NIR.vrt*

  mainly because this is taken over for the extracted data

  *data/processed/L3/extracted/32UNU/clc2018_lte50ha_32UNU/32unu__scoll01__vts4w__2018-01-07__NIR.npy*.

  When concatenating the latter later over mutliple tiles in a dataframe we need to remove the tile identifier anyway.
  Else the column names differ and we cannot concatenate.

