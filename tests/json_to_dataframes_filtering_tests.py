from json_to_dataframes import json_to_dataframes
from filtering_files import filtering_files

dfs = json_to_dataframes("data/data.json")
dfs_filtered = filtering_files(dfs, 4, 30, 2, 8, 100.00, 500.00, 0.0, 10000.0, [False, True, False, False, True, False])

dfs[0].to_excel("unfiltered_1.xlsx")
dfs[1].to_excel("unfiltered_2.xlsx")

dfs_filtered[0].to_excel("filtered_1.xlsx")
dfs_filtered[1].to_excel("filtered_2.xlsx")
