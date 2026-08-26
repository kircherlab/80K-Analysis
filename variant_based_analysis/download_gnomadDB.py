from gnomad_db.database import gnomAD_DB
download_link = "https://zenodo.org/record/6818606/files/gnomad_db_v3.1.2.sqlite3.gz?download=1"
output_dir = "data/gnomadDB/v3.1.2" # database_location
gnomAD_DB.download_and_unzip(download_link, output_dir)
