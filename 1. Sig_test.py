"""
This component of the Empirical Neighbourhood Calibration method conducts the
significance test on the specified data.
"""
import gdal
import numpy as np
import math
from considered_distances import considered_distances
from neighbourhood_evaluator import neighbourhood_evaluator
from MWU_test import mwu_test


# Specify the base path to the directory containing the empirical neighbourhood
# calibration tool-pack.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\P3_sampling\\"
# Select an example case study application. Specify the name below:
case_study_path = "Randstad_simpler\\Randstad\\"
# Set the paths to the directories and relevant data
data_path = base_path + case_study_path + "Data\\"
output_path = base_path + "Output\\"
map1_path = data_path + "lu1989.asc"
map2_path = data_path + "lu2000.asc"
mask_path = data_path + "region.asc"
# Set the land-use class names.
luc_names = ["Agriculture", "Greenhouses", "Residential", "Industry",
             "Services", "Socio-cultural uses", "Nature",
             "Recreation areas", "Airport", "Water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 2
act = luc - (pas + fea)
# Specify the maximum neighbourhood size distance considered
max_distance = 8

# Read in the map for the data at time slice 1.
src_map = gdal.Open(map1_path)
omap = np.array(src_map.GetRasterBand(1).ReadAsArray())
# Read in the map for the data at time slice 2.
src_map = gdal.Open(map2_path)
amap = np.array(src_map.GetRasterBand(1).ReadAsArray())
# Read in the masking map.
src_map = gdal.Open(mask_path)
mask = np.array(src_map.GetRasterBand(1).ReadAsArray())

# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]

# Determine the distances that will be analysed,
# use module: considered_distances.

temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])

# Conduct the significance testing using the module 'neighbourhood_evaluator.'
dummy = neighbourhood_evaluator(
    luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols
)
# Store the requisite output into specific dictionaries
all_cells_baseline = dummy[0]
no_new_cells_ci_baseline = dummy[1]
no_cells_ci_baseline = dummy[2]
transition_dictionary = dummy[3]
ef = dummy[4]
# Log scale the enrichment factor values.
log_ef = np.zeros(shape=(max_distance, luc, luc))
for p in range(0, luc):
    for q in range(0, luc):
        for c in range(0, max_distance):
            # If the enrichment factor is not evaluated a value of 0 is given.
            # Hence, a did not evaluate value of -9999 is used.
            if ef[c, p, q] == 0:
                log_ef[c, p, q] = -9999
            else:
                log_ef[c, p, q] = math.log(ef[c, p, q], 10)

# Specify what the significance limit is (recommended value is 1.96 i.e. 95%).
z_limit = 1.96
# Conduct the MWU test using the module 'mwu_test.'
z_scores = mwu_test(
    max_distance, luc, transition_dictionary, all_cells_baseline, N
)

# Determine which rules have meaningful over-representation at nearby distances.
sig_distances = np.zeros(shape=(act, luc))
c = 1
for i in range(0, act):
    for j in range(0, luc):
        if abs(z_scores[c, i+pas, j]) > z_limit and log_ef[c, i+pas, j] > 0:
            sig_distances[i, j] = sig_distances[i, j] + 1
c = 2
for i in range(0, act):
    for j in range(0, luc):
        if abs(z_scores[c, i + pas, j]) > z_limit and log_ef[c, i+pas, j] > 0:
            sig_distances[i, j] = sig_distances[i, j] + 1
att_rules = np.zeros(shape=(act, luc))
for i in range(0, act):
    for j in range(0, luc):
        if i + pas == j:
            att_rules[i, j] = 1
        elif sig_distances[i, j] > 1:
            att_rules[i, j] = 1
att_rules = np.transpose(att_rules)

# Save the attraction rules to a file.
att_rules_file_path = output_path + "att_rules.txt"
np.savetxt(att_rules_file_path, att_rules, fmt="%d")
