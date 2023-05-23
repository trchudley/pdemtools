import tempfile
import urllib.request
import geopandas as gpd

print("Downloading index file (this may take some time)")

url = "https://data.pgc.umn.edu/elev/dem/setsm/REMA/indexes/REMA_Strip_Index_s2s041_gpkg.gpkg"

with tempfile.NamedTemporaryFile(suffix=".gpkg") as temp:
    urllib.request.urlretrieve(url, temp.name)

    print("Exporting to parquet file in current directory...")
    gdf = gpd.read_file(temp.name)
    gdf.to_parquet("REMA_Strip_Index_s2s041.parquet")

print("Finished.")
