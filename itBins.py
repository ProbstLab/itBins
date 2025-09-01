###############################################################################
###############################################################################
# HEADER
###############################################################################


###############################################################################
###############################################################################
# IMPORTS
###############################################################################

# from core
import argparse
import time
import copy
import sys
import math
import json
# import warnings
import traceback

# from packages
import pandas as pd
import numpy as np


###############################################################################
###############################################################################
# settings
###############################################################################

# warnings.filterwarnings('ignore')
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", 2000)
pd.set_option("display.width", 160)
pd.set_option("display.max_colwidth", 300)
# pd.options.mode.chained_assignment = "raise"
pd.options.mode.chained_assignment = None
np.set_printoptions(threshold=10000)

###############################################################################
###############################################################################
# GLOBALS
###############################################################################

BACTERIAL_REFERENCE_GENES = ["B_Histidyl-tRNA_synthetase",
                             "B_Phenylalanyl-tRNA_synthetase_alpha",
                             "B_Preprotein_translocase_subunit_SecY",
                             "B_Valyl-tRNA_synthetase",
                             "B_alanyl_tRNA_synthetase",
                             "B_arginyl_tRNA_synthetase",
                             "B_aspartyl_tRNA_synthetase",
                             "B_gyrA",
                             "B_leucyl-tRNA_synthetase",
                             "B_recA",
                             "B_ribosomal_protein_L1",
                             "B_ribosomal_protein_L10",
                             "B_ribosomal_protein_L11",
                             "B_ribosomal_protein_L13",
                             "B_ribosomal_protein_L14",
                             "B_ribosomal_protein_L15",
                             "B_ribosomal_protein_L16-L10E",
                             "B_ribosomal_protein_L17",
                             "B_ribosomal_protein_L18",
                             "B_ribosomal_protein_L19",
                             "B_ribosomal_protein_L2",
                             "B_ribosomal_protein_L20",
                             "B_ribosomal_protein_L21",
                             "B_ribosomal_protein_L22",
                             "B_ribosomal_protein_L23",
                             "B_ribosomal_protein_L24",
                             "B_ribosomal_protein_L27",
                             "B_ribosomal_protein_L29",
                             "B_ribosomal_protein_L3",
                             "B_ribosomal_protein_L30",
                             "B_ribosomal_protein_L4",
                             "B_ribosomal_protein_L5",
                             "B_ribosomal_protein_L6P-L9E",
                             "B_ribosomal_protein_S10",
                             "B_ribosomal_protein_S11",
                             "B_ribosomal_protein_S12",
                             "B_ribosomal_protein_S13",
                             "B_ribosomal_protein_S15",
                             "B_ribosomal_protein_S16",
                             "B_ribosomal_protein_S17",
                             "B_ribosomal_protein_S18",
                             "B_ribosomal_protein_S19",
                             "B_ribosomal_protein_S2",
                             "B_ribosomal_protein_S20",
                             "B_ribosomal_protein_S3",
                             "B_ribosomal_protein_S4",
                             "B_ribosomal_protein_S5",
                             "B_ribosomal_protein_S6",
                             "B_ribosomal_protein_S7",
                             "B_ribosomal_protein_S8",
                             "B_ribosomal_protein_S9"]

BACTERIAL_REFERENCE_GENE_NUMBER = 51

ARCHAEAL_REFERENCE_GENES = ["A_CCA-adding_enzyme",
                            "A_DNA-directed_RNA_polymerase",
                            "A_DNA-directed_RNA_polymerase_subunit_N",
                            "A_Dimethyladenosine_transferase",
                            "A_Diphthamide_biosynthesis_protein",
                            "A_Fibrillarin-like_rRNA/tRNA_2'-O-methyltransferase",
                            "A_Glycyl-tRNA_synthetase",
                            "A_KH_type_1_domain_protein",
                            "A_Methionyl-tRNA_synthetase",
                            "A_Non-canonical_purine_NTP_pyrophosphatase",
                            "A_PUA_domain_containing_protein",
                            "A_Phenylalanyl-tRNA_synthetase_alpha_subunit",
                            "A_Phenylalanyl-tRNA_synthetase_beta_subunit",
                            "A_Pre-mRNA_processing_ribonucleoprotein",
                            "A_Prolyl-tRNA_synthetase",
                            "A_Protein_pelota_homolog",
                            "A_Ribosomal_protein_L10e",
                            "A_Ribosomal_protein_L13",
                            "A_Ribosomal_protein_L18e",
                            "A_Ribosomal_protein_L21e",
                            "A_Ribosomal_protein_L3",
                            "A_Ribosomal_protein_L7Ae/L8e",
                            "A_Ribosomal_protein_S13",
                            "A_Ribosomal_protein_S15",
                            "A_Ribosomal_protein_S19e",
                            "A_Ribosomal_protein_S2",
                            "A_Ribosomal_protein_S28e",
                            "A_Ribosomal_protein_S3Ae",
                            "A_Ribosomal_protein_S6e",
                            "A_Ribosomal_protein_S7",
                            "A_Ribosomal_protein_S9",
                            "A_Ribosome_maturation_protein_SDO1_homolog",
                            "A_Signal_recognition_particle_54_kDa_protein",
                            "A_Transcription_elongation_factor_Spt5",
                            "A_Translation_initiation_factor_5A",
                            "A_Translation_initiation_factor_IF-2_subunit_gamma",
                            "A_Valyl-tRNA_synthetase",
                            "A_tRNA_N6-adenosine_threonylcarbamoyltransferase"]

ARCHAEAL_REFERENCE_GENE_NUMBER = 38


#dt_timer = 0
#dt_count = 0
#dt_np_timer = 0
#dt_pd_timer = 0
#dt_np_count = 0

#GC_timer = 0
#FAST_GC_timer = 0
#Cov_timer = 0
#FAST_Cov_timer = 0

#task_m1_timer = 0
#task_0_timer = 0
#task_1_timer = 0
#task_2_timer = 0
#task_3_timer = 0
#task_4_timer = 0
#task_5_timer = 0
#task_6_timer = 0
#task_7_timer = 0
#task_71_timer = 0
#task_8_timer = 0
#task_9_timer = 0


###############################################################################
###############################################################################
# FUNCTIONS
###############################################################################

###############################################################################
# FAST DASTOOL BIN SCORE CALCULATION
###############################################################################

def FAST_dt_score(array_of_SCGs):
    """
    calculates the DAStool bin score for the supplied np.array of singel copy
    genes, agnostic of taxonomy and without applying filters
    returns float bin_score
    """
    array_of_SCG_sums = np.sum(array_of_SCGs, axis=1)
    reference_SCG_number = np.shape(array_of_SCG_sums)[0]
    total_SCG_number = np.sum(array_of_SCG_sums)
    unique_SCG_number = np.shape(array_of_SCG_sums[array_of_SCG_sums > 0])[0]
    duplicated_SCG_number = np.shape(array_of_SCG_sums[array_of_SCG_sums > 1])[0]
    try:
        bin_score = ((unique_SCG_number/reference_SCG_number) -
                     0.5*(duplicated_SCG_number/unique_SCG_number) -
                     0.5*((total_SCG_number-unique_SCG_number)/reference_SCG_number))
    except ZeroDivisionError:
        bin_score = 0
    return bin_score


###############################################################################
# FAST COMPLETENESS CALCULATION
###############################################################################

def FAST_completeness(array_of_SCGs):
    """
    calculates SCG completeness for the supplied np.array of singel copy
    genes, agnostic of taxonomy and without applying filters
    returns float bin_completeness
    """
    array_of_SCG_sums = np.sum(array_of_SCGs, axis=1)
    reference_SCG_number = np.shape(array_of_SCG_sums)[0]
    unique_SCG_number = np.shape(array_of_SCG_sums[array_of_SCG_sums > 0])[0]
    bin_completeness = unique_SCG_number / reference_SCG_number
    return bin_completeness


###############################################################################
# FAST COMPLETENESS CALCULATION
###############################################################################

def FAST_contamination(array_of_SCGs):
    """
    calculates SCG contamination for the supplied np.array of singel copy
    genes, agnostic of taxonomy and without applying filters
    returns float bin_contamination
    """
    array_of_SCG_sums = np.sum(array_of_SCGs, axis=1)
    reference_SCG_number = np.shape(array_of_SCG_sums)[0]
    unique_SCG_number = np.shape(array_of_SCG_sums[array_of_SCG_sums > 0])[0]
    excess_SCG_number = np.sum(array_of_SCG_sums) - unique_SCG_number
    bin_contamination = excess_SCG_number / reference_SCG_number
    return bin_contamination


###############################################################################
# FAST FRAME FILTERING
###############################################################################

def FAST_filter(input_array,
                array_of_one_or_three_filters,
                filter_frame_by_GC=True,
                filter_frame_by_Coverage=True,
                filter_frame_by_Taxonomy=True):
    """
    filters an array by one or three filters, returning a copy that is sliced
    to wherever all the filters != 0
    returns np.array
    """
    if array_of_one_or_three_filters.ndim > 1:
        array_of_filter = np.ones((np.shape(input_array)[-1]), dtype="int8")
        if (filter_frame_by_GC | filter_frame_by_Coverage | filter_frame_by_Taxonomy):
            if filter_frame_by_GC:
                array_of_filter = np.logical_and(array_of_filter, array_of_one_or_three_filters[0, :])
            if filter_frame_by_Coverage:
                array_of_filter = np.logical_and(array_of_filter, array_of_one_or_three_filters[1, :])
            if filter_frame_by_Taxonomy:
                array_of_filter = np.logical_and(array_of_filter, array_of_one_or_three_filters[2, :])
        if np.ndim(input_array) > 1:
            output_array = input_array[:, array_of_filter == 1]
        else:
            output_array = input_array[array_of_filter == 1]
    else:
        if np.ndim(input_array) > 1:
            output_array = input_array[:, array_of_one_or_three_filters == 1]
        else:
            output_array = input_array[array_of_one_or_three_filters == 1]
    return output_array


###############################################################################
# FAST ENTRY (PROTOTYPE)
###############################################################################

def PROT_entry(frame):
    """
    creates an entrypoint for FAST functions
    returns dict of np.arrays
    """
    array_filters = np.transpose(frame.loc[:, ["GC_Filter",
                                               "Cov_Filter",
                                               "Tax_Filter"]].to_numpy(dtype="bool"))
    array_GC = frame.loc[:, "GC"].to_numpy()
    array_Cov = frame.loc[:, "coverage"].to_numpy()
    array_length = frame.loc[:, "length"].to_numpy()
    array_bac = np.transpose(frame.loc[:, BACTERIAL_REFERENCE_GENES].to_numpy())
    array_arc = np.transpose(frame.loc[:, ARCHAEAL_REFERENCE_GENES].to_numpy())
    array_tax = np.transpose(frame.loc[:, ["taxL1",
                                           "taxL2",
                                           "taxL3",
                                           "taxL4",
                                           "taxL5",
                                           "taxL6"]].to_numpy())
    array_tax[array_tax == ""] = "unclassified"
    dict_of_arrays = {"filters": array_filters,
                      "GC": array_GC,
                      "Cov": array_Cov,
                      "Tax": array_tax,
                      "length": array_length,
                      "bacterial_SCGs": array_bac,
                      "archaeal_SCGs": array_arc}
    return dict_of_arrays


###############################################################################
# FAST COMP FILTER (PROTOTYPE)
###############################################################################

def PROT_completeness_cutoff(dict_of_arrays,
                             cutoff):
    """
    checks for higher completeness and compares it to the cutoff
    returns dict of np.arrays, weather it did something and the completenesses
    """
    bacterial_completeness = FAST_completeness(dict_of_arrays["bacterial_SCGs"])
    archaeal_completeness = FAST_completeness(dict_of_arrays["archaeal_SCGs"])
    higher_completeness = max(bacterial_completeness, archaeal_completeness)
    if higher_completeness < cutoff:
        dict_of_arrays["filters"][0:2, :] = False
        return dict_of_arrays, True, bacterial_completeness, archaeal_completeness
    return dict_of_arrays, False, bacterial_completeness, archaeal_completeness


###############################################################################
# FAST CONT FILTER (PROTOTYPE)
###############################################################################

def PROT_contamination_cutoff(dict_of_arrays,
                              cutoff):
    """
    checks for higher completeness and compares the higher completeness contamination to the cutoff
    returns dict of np.arrays, weather it did something and the completenesses
    """
    bacterial_completeness = FAST_completeness(FAST_filter(dict_of_arrays["bacterial_SCGs"], dict_of_arrays["filters"]))
    archaeal_completeness = FAST_completeness(FAST_filter(dict_of_arrays["archaeal_SCGs"], dict_of_arrays["filters"]))
    relevant_contamination = FAST_contamination(FAST_filter(dict_of_arrays[("bacterial_SCGs" if (bacterial_completeness > archaeal_completeness) else "archaeal_SCGs")], dict_of_arrays["filters"]))
    if relevant_contamination > cutoff:
        dict_of_arrays["filters"][:, :] = False
        return dict_of_arrays, True, relevant_contamination, ("bacterial_SCGs" if (bacterial_completeness > archaeal_completeness) else "archaeal_SCGs")
    return dict_of_arrays, False, relevant_contamination, ("bacterial_SCGs" if (bacterial_completeness > archaeal_completeness) else "archaeal_SCGs")


###############################################################################
# FAST REMOVE TAX (PROTOTYPE)
###############################################################################

def PROT_remove_taxonomy(dict_of_arrays,
                         remove_viruses=True,
                         remove_eukaryota=True,
                         remove_bacteria=False,
                         remove_archaea=False,
                         remove_unclassified=False,
                         remove_other=False):
    """
    removes specified taxonomies from dict_of_arrays
    returns dict of np.arrays
    """
    if remove_viruses:
        dict_of_arrays["filters"][2, np.where(dict_of_arrays["Tax"][0, :] == "Viruses")[0]] = False
    if remove_eukaryota:
        dict_of_arrays["filters"][2, np.where(dict_of_arrays["Tax"][0, :] == "Eukaryota")[0]] = False
    if remove_bacteria:
        dict_of_arrays["filters"][2, np.where(dict_of_arrays["Tax"][0, :] == "Bacteria")[0]] = False
    if remove_archaea:
        dict_of_arrays["filters"][2, np.where(dict_of_arrays["Tax"][0, :] == "Archaea")[0]] = False
    if remove_unclassified:
        dict_of_arrays["filters"][2, np.where(dict_of_arrays["Tax"][0, :] == "unclassified")[0]] = False
    if remove_other:
        dict_of_arrays["filters"][2, (np.where((dict_of_arrays["Tax"][0, :] != "Viruses") &
                                               (dict_of_arrays["Tax"][0, :] != "Eukaryota") &
                                               (dict_of_arrays["Tax"][0, :] != "Bacteria") &
                                               (dict_of_arrays["Tax"][0, :] != "Archaea") &
                                               (dict_of_arrays["Tax"][0, :] != "unclassified")))[0]] = False
    return dict_of_arrays


###############################################################################
# FAST PEAKLIST (PROTOTYPE)
###############################################################################

def PROT_create_peaklist(dict_of_arrays, type_of_peaklist):
    """
    create a peaklist
    returns dictionarys of np.arrays
    """
    if type_of_peaklist == "GC":
        a_values = dict_of_arrays["GC"]
    elif type_of_peaklist == "Cov":
        a_values = dict_of_arrays["Cov"]
    else:
        exception_message = "FAST_create_peaklist: Illegal argument \"" + str(type_of_peaklist) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\"."
        raise Exception(exception_message)
        traceback.print_exc()
    a_filters = dict_of_arrays["filters"]
    a_length = dict_of_arrays["length"]
    a_bac = dict_of_arrays["bacterial_SCGs"]
    a_arc = dict_of_arrays["archaeal_SCGs"]

    af_values = FAST_filter(a_values, a_filters)
    af_length = FAST_filter(a_length, a_filters)
    af_bac = FAST_filter(a_bac, a_filters)
    af_arc = FAST_filter(a_arc, a_filters)

    maximum_value = int(np.max(af_values))+2
    minimum_value = int(np.min(af_values))-1
    value_range = maximum_value-minimum_value

    h_list_of_values = np.arange(minimum_value, maximum_value, dtype="int")
    h_list_of_lengths = np.empty(value_range, dtype="int")

    for i, item in np.ndenumerate(h_list_of_values):
        h_list_of_lengths[i] = np.sum(af_length[(af_values >= item) & (af_values < (item+1))])
    h_lookahead_list = np.roll(h_list_of_lengths, -1)
    h_lookahead_list[-1] = 0

    h_list_of_deltas = h_lookahead_list - h_list_of_lengths

    h_list_of_parts = np.empty(value_range, dtype="int")
    h_list_of_parts[h_list_of_deltas < 0] = 3
    h_list_of_parts[h_list_of_deltas > 0] = 2
    h_list_of_parts[h_list_of_deltas == 0] = 1
    h_list_of_parts[(h_list_of_deltas == 0) & (h_list_of_lengths == 0)] = 0
    if h_list_of_lengths[1] > 0:
        h_list_of_parts[0] = 2
    else:
        h_list_of_parts[0] = 0

    h_list_of_parts_lookbehind = np.roll(h_list_of_parts, 1)
    h_list_of_peaks = np.zeros(value_range, dtype="int")
    h_list_of_peaks[h_list_of_parts_lookbehind == 0] = -2
    h_list_of_peaks[(h_list_of_parts_lookbehind == 3) & (h_list_of_parts == 2)] = -1
    h_list_of_peaks[(h_list_of_parts_lookbehind == 2) & (h_list_of_parts == 3)] = 1
    h_list_of_peaks[h_list_of_lengths == 0] = -2
    h_list_of_peaks[0] = -2

    h_list_of_inclusions = np.zeros(value_range, dtype="int")

    histogram = {type_of_peaklist: h_list_of_values,
                 "length": h_list_of_lengths,
                 "delta": h_list_of_deltas,
                 "parts": h_list_of_parts,
                 "peaks": h_list_of_peaks,
                 "included": h_list_of_inclusions}
    p_list_of_positions = h_list_of_values[h_list_of_peaks == 1]

    bacterial_score = FAST_dt_score(af_bac)
    archaeal_score = FAST_dt_score(af_arc)

    num_of_peaks = np.shape(p_list_of_positions)[0]

    p_list_of_peak_starts = np.empty(num_of_peaks, dtype="int")
    p_list_of_peak_ends = np.empty(num_of_peaks, dtype="int")
    p_list_of_left_separators = np.empty(num_of_peaks, dtype="int")
    p_list_of_right_separators = np.empty(num_of_peaks, dtype="int")
    p_list_of_peak_areas = np.empty(num_of_peaks, dtype="int")
    p_list_of_bacterial_score_improvements_on_drop = np.empty(num_of_peaks)
    p_list_of_archaeal_score_improvements_on_drop = np.empty(num_of_peaks)
    p_list_of_inclusions = np.zeros(num_of_peaks, dtype="int")

    for i, item in np.ndenumerate(p_list_of_positions):
        p_list_of_peak_starts[i] = np.max(h_list_of_values[(h_list_of_values < item) & (h_list_of_peaks < 0)])
        p_list_of_peak_ends[i] = np.min(h_list_of_values[(h_list_of_values > item) & (h_list_of_peaks < 0)])

    for i, item in np.ndenumerate(p_list_of_peak_starts):
        p_list_of_left_separators[i] = h_list_of_peaks[item-minimum_value]

    for i, item in np.ndenumerate(p_list_of_peak_ends):
        p_list_of_right_separators[i] = h_list_of_peaks[item-minimum_value]

    for i in range(num_of_peaks):
        p_list_of_peak_areas[i] = np.sum(h_list_of_lengths[(p_list_of_peak_starts[i] - minimum_value + 1):(p_list_of_peak_ends[i] - minimum_value)])  # + 1 makes both ends
        p_list_of_bacterial_score_improvements_on_drop[i] = FAST_dt_score(af_bac[:, (af_values <= p_list_of_peak_starts[i]) | (af_values >= p_list_of_peak_ends[i])]) - bacterial_score

        p_list_of_archaeal_score_improvements_on_drop[i] = FAST_dt_score(af_arc[:, (af_values <= p_list_of_peak_starts[i]) | (af_values >= p_list_of_peak_ends[i])]) - archaeal_score

    p_list_of_inclusions = np.array([0] * num_of_peaks)

    array_peaklist = {"peak_position": p_list_of_positions,
                      "start": p_list_of_peak_starts,
                      "end": p_list_of_peak_ends,
                      "area": p_list_of_peak_areas,
                      "sep_left": p_list_of_left_separators,
                      "sep_right": p_list_of_right_separators,
                      "bacdropimp": p_list_of_bacterial_score_improvements_on_drop,
                      "arcdropimp": p_list_of_archaeal_score_improvements_on_drop,
                      "included": p_list_of_inclusions}

    return array_peaklist, histogram


###############################################################################
# FAST PEAK PICKER (PROTOTYPE)
###############################################################################

def PROT_pick_peak(array_peaklist, histogram, type_of_peaklist):
    """
    select highest area peak down to sep
    returns array peaklist and histogram
    """
    index_largest_peak = np.where(array_peaklist["area"] == np.max(array_peaklist["area"]))[0]
    if np.shape(index_largest_peak)[0] == 1:
        # left_edge = array_peaklist["start"][index_largest_peak]
        # right_edge = array_peaklist["end"][index_largest_peak]
        array_peaklist["included"][index_largest_peak] = 1
        if type_of_peaklist == "GC":
            histogram["included"][(histogram["GC"] >= array_peaklist["start"][index_largest_peak]) & (histogram["GC"] <= array_peaklist["end"][index_largest_peak])] = 1
        elif type_of_peaklist == "GC":
            histogram["included"][(histogram["Cov"] >= array_peaklist["start"][index_largest_peak]) & (histogram["Cov"] <= array_peaklist["end"][index_largest_peak])] = 1
        return array_peaklist, histogram
    else:
        print("Well shit, i cannot deal with multiple peaks with the same area yet. (PROT_GC_pick)")


###############################################################################
# FAST PEAK EXPANDER (PROTOTYPE)
###############################################################################

def PROT_expand_to_baseline(array_peaklist, histogram, type_of_peaklist):
    """
    returns expanded peaklist and histogram
    """
    leftest_inclusion = np.min(np.where(histogram["included"] == 1)[0])
    closest_left_baseline = np.max(np.where((histogram["peaks"] == -2) & (histogram[type_of_peaklist] <= [histogram[type_of_peaklist][leftest_inclusion]]))[0])
    if leftest_inclusion != closest_left_baseline:
        histogram["included"][closest_left_baseline: (leftest_inclusion)] = 1
    rightest_inclusion = np.max(np.where(histogram["included"] == 1)[0])
    closest_right_baseline = np.min(np.where((histogram["peaks"] == -2) & (histogram[type_of_peaklist] >= [histogram[type_of_peaklist][rightest_inclusion]]))[0])
    if leftest_inclusion != closest_right_baseline:
        histogram["included"][(rightest_inclusion + 1):closest_right_baseline+1] = 1
    return array_peaklist, histogram


###############################################################################
# FAST PEAK LIMITER (PROTOTYPE)
###############################################################################

def PROT_limit_peak(dict_of_arrays,
                    array_peaklist,
                    histogram,
                    fixed_radius,
                    linear_peak_position_radius,
                    linear_bin_size_radius,
                    filter_type):
    selected_peak = array_peaklist["peak_position"][array_peaklist["included"] == 1]  # gc value of peak
    bin_size = np.shape(dict_of_arrays["length"])[0]
    right_edge = (selected_peak +
                  fixed_radius +
                  linear_peak_position_radius * selected_peak +
                  linear_bin_size_radius * bin_size)
    left_edge = (selected_peak -
                 fixed_radius -
                 linear_peak_position_radius * selected_peak -
                 linear_bin_size_radius * bin_size)
    if filter_type == "GC" :
        dict_of_arrays["filters"][0, (dict_of_arrays["GC"] < left_edge) | (dict_of_arrays["GC"] > right_edge)] = False
    elif filter_type == "Cov" :
        dict_of_arrays["filters"][1, (dict_of_arrays["Cov"] < left_edge) | (dict_of_arrays["Cov"] > right_edge)] = False
    else :
        exception_message = "FAST_limit_peak_GC: Illegal argument \"" + str(filter_type) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\"."
        raise Exception(exception_message)
    return dict_of_arrays


###############################################################################
# FAST FILTER UPDATER (PROTOTYPE)
###############################################################################

def PROT_filter_updater(array_peaklist, dict_of_arrays, filter_to_update):
    if filter_to_update == "GC":
        array_peaklist
    elif filter_to_update == "Cov":
        pass
    elif filter_to_update == "Tax":
        pass
    frame.loc[(frame["GC"] < ledge) | (frame["GC"] > redge), "GC_Filter"] = 0
    frame.loc[(frame["GC"] >= ledge) & (frame["GC"] <= redge), "GC_Filter"] = 1
    return dict_of_arrays


###############################################################################
# FAST TAX_DROP
###############################################################################

def PROT_tax_drop(dict_of_arrays,
                  drop_level=1,
                  max_length=1.0,
                  always_drop_if_score_improves=True):
    """
    drops a tax
    """
    # check if making a taxtable is even neccessary
    if np.shape(FAST_filter(dict_of_arrays["length"], dict_of_arrays["filters"]))[0] is False:
        return dict_of_arrays, False
    tax_table, bacterial_score, archaeal_score, empty = PROT_tax_table(dict_of_arrays, tax_level=drop_level)
    # check if tax_table is empty
    if empty:
        iprint("No tax dropped")
        return dict_of_arrays, False
    bacterial_completeness = FAST_completeness(FAST_filter(dict_of_arrays["bacterial_SCGs"], dict_of_arrays["filters"]))
    archaeal_completeness = FAST_completeness(FAST_filter(dict_of_arrays["archaeal_SCGs"], dict_of_arrays["filters"]))
    drop_length = max_length * np.sum(FAST_filter(dict_of_arrays["length"], dict_of_arrays["filters"]))  # check boolification
    chosen_domain = "bacterial_score_improvement" if (bacterial_completeness > archaeal_completeness) else "archaeal_score_improvement"
    highest_score_improvement = np.max(tax_table[chosen_domain])
    choose_by_score_improvement = np.where(tax_table[chosen_domain] == highest_score_improvement)
    if highest_score_improvement == 0:
        for item in np.nditer(choose_by_score_improvement):
            if (tax_table["length"][item] <= drop_length) | always_drop_if_score_improves:
                tax_to_drop = tax_table["tax"][item]
                iprint("Drop: ", tax_to_drop)
                dict_of_arrays["filters"][2, (dict_of_arrays["Tax"][drop_level] == tax_to_drop)] = False
            else:
                iprint("No tax dropped")
        iprint("Done dropping")
        return dict_of_arrays, False

    choose_by_length = np.where(tax_table["length"] == np.max(tax_table["length"][choose_by_score_improvement]))[0]
    if tax_table[chosen_domain][choose_by_length] < 0:
        iprint("No tax dropped")
        return dict_of_arrays, False
    if (tax_table["length"][choose_by_length] <= drop_length) | always_drop_if_score_improves:
        tax_to_drop = tax_table["tax"][choose_by_length]
        iprint("Drop: ", tax_to_drop)
        dict_of_arrays["filters"][2, (dict_of_arrays["Tax"][drop_level] == tax_to_drop)] = False
        return dict_of_arrays, True
    else:
        iprint("No tax dropped")
        return dict_of_arrays, False


###############################################################################
# FAST TAX_TABLE
###############################################################################

def PROT_tax_table(dict_of_arrays,
                   tax_level=1):
    """
    creates tax table
    """
    array_of_tax = FAST_filter(dict_of_arrays["Tax"][tax_level], dict_of_arrays["filters"])
    array_of_lengths = FAST_filter(dict_of_arrays["length"], dict_of_arrays["filters"])
    array_of_bac_SCGs = FAST_filter(dict_of_arrays["bacterial_SCGs"], dict_of_arrays["filters"])
    array_of_arc_SCGs = FAST_filter(dict_of_arrays["archaeal_SCGs"], dict_of_arrays["filters"])
    tt_tax_array = np.unique(array_of_tax)
    num_of_tax = np.shape(tt_tax_array)[0]
    tt_length_array = np.empty(num_of_tax, dtype="int")
    tt_bacterial_score_array = np.empty(num_of_tax, dtype="float")
    tt_archaeal_score_array = np.empty(num_of_tax, dtype="float")
    tt_bacterial_score_improvement_array = np.empty(num_of_tax, dtype="float")
    tt_archaeal_score_improvement_array = np.empty(num_of_tax, dtype="float")
    bacterial_score = FAST_dt_score(array_of_bac_SCGs)
    archaeal_score = FAST_dt_score(array_of_arc_SCGs)
    for i, item in np.ndenumerate(tt_tax_array):
        tt_length_array[i] = np.sum(array_of_lengths[array_of_tax == item])
        tt_bacterial_score_array[i] = FAST_dt_score(array_of_bac_SCGs[:, array_of_tax == item])
        tt_archaeal_score_array[i] = FAST_dt_score(array_of_arc_SCGs[:, array_of_tax == item])
        tt_bacterial_score_improvement_array[i] = FAST_dt_score(array_of_bac_SCGs[:, array_of_tax != item]) - bacterial_score
        tt_archaeal_score_improvement_array[i] = FAST_dt_score(array_of_arc_SCGs[:, array_of_tax != item]) - archaeal_score
    tax_table = {"tax": tt_tax_array,
                 "length": tt_length_array,
                 "bacterial_score": tt_bacterial_score_array,
                 "archaeal_score": tt_archaeal_score_array,
                 "bacterial_score_improvement": tt_bacterial_score_improvement_array,
                 "archaeal_score_improvement": tt_archaeal_score_improvement_array}
    return tax_table, bacterial_score, archaeal_score, False


###############################################################################
# PRINTING
###############################################################################

def iprint(*args, out=sys.stdout):
    if cl_args.log_to_stderr:
        out = sys.stderr
    if info:
        print(*args, file=out)


def lprint(*args, out=sys.stdout):
    if cl_args.log_to_stderr:
        out = sys.stderr
    if loud:
        print(*args, file=out)


###############################################################################
###############################################################################
###############################################################################

#   Arguments

###############################################################################

parser = argparse.ArgumentParser(prog="itBins",
                                 add_help=False,
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 allow_abbrev = False,
                                 description="╔══════════════════════════════════════════════════════════════════════════════╗\n"
                                             "║      ITBINS - ALGORITHMIC BIN CURATOR                                        ║\n"
                                             "╚══════════════════════════════════════════════════════════════════════════════╝\n"
                                             " Performs bin curation using GC content, coverage and taxonomy. Expect about 5\n"
                                             " seconds per bin, this varies slightly with the number of scaffolds in the bin\n"
                                             " and bins with a lot of diversity at taxonomic level 2 (usually phylum) can\n"
                                             " take way longer to curate, as the tax_drop task is expensive compared to the\n"
                                             " pick_expand_limit tasks of GC and cov.\n\n"

                                             "╔═█\033[7mGETTING STARTED\033[0m█════════════════════════════════════════════════════════════╗\n"
                                             "║ Installing:                                                                  ║\n"
                                             "║                                                                              ║\n"
                                             "║     conda create -n \"itBinsEnv\" python=3.10.11 pandas=1.4.2 numpy=1.21.5     ║\n"
                                             "║                                                                              ║\n"
                                             "║ mamba or micromamba may be better suited to resolve the environment          ║\n"
                                             "║                                                                              ║\n"
                                             "║                                                                              ║\n"
                                             "║ You just want to curate some bins:                                           ║\n"
                                             "║                                                                              ║\n"
                                             "║ (You will need the files produced by script XXX)                             ║\n"
                                             "║                                                                              ║\n"
                                             "║     python itBins.py --example-task-file                                     ║\n"
                                             "║     python itBins.py -ost -b overview.txt -g SCG.csv > itBins.log            ║\n"
                                             "║                                                                              ║\n"
                                             "║ The tasks.json file is itBins' config. You can modify it to tweak curation.  ║\n"
                                             "╚══════════════════════════════════════════════════════════════════════════════╝\n"
                                             ""

                                             " Copyright Julian M. Künkel. Licensed under the EUPL-1.2.\n"
                                             " Written at UDE-GEM\n\n"
                                             " Requires:\n"
                                             "\tpython 3.10.11\n"
                                             "\tpandas 1.4.2\n"
                                             "\tnumpy  1.21.5",
                                 epilog="If you have any problems or suggestions, send them through slack\n"
                                        "or contact me using my work email.\n")
parser.add_argument("-b", "--bins",
                    dest="ovw_path",
                    help="provide path to an overview file (.tsv, TAB-separated\n"
                         "with a header). This is the same file you would load\n"
                         "into uBin.\n\n",
                    metavar="BIN_OVERVIEW_FILE",
                    required=False,
                    type=str)
parser.add_argument("-g", "--genes",
                    dest="scg_path",
                    help="provide path to a singe copy gene file (.csv, comma-\n"
                         "separated with a header). This is the same file you\n"
                         "would load into uBin.\n\n",
                    metavar="SINGLE_COPY_GENE_FILE",
                    type=str)
parser.add_argument("-d", "--bin-DASTool",
                    default="Bin",
                    dest="bdt_col",
                    help="provide the name of the column in the overview file\n"
                         "containing the binnames of the uncurated bins,\n"
                         "defaults to \"Bin\"\n\n",
                    metavar="DASTOOL_BIN_NAME_COLUMN",
                    type=str)
parser.add_argument("-u", "--bin-uBin",
                    default="curated",
                    dest="bub_col",
                    help="provide the name of the column in the overview file\n"
                         "containing the binnames of the curated bins, defaults\n"
                         "to \"curated\"\n\n",
                    metavar="UBIN_BIN_NAME_COLUMN",
                    type=str)
parser.add_argument("-p", "--prefix",
                    default="",
                    dest="name_prefix",
                    help="provide a prefix to attach to the names of bins and\n"
                         "scaffolds. this is usefull, if multiple files are\n"
                         "processed and names are not unique between them.\n"
                         "Usually, this should be the name of the sample or\n"
                         "project, followed by \"_\"\n\n",
                    metavar="PREFIX",
                    type=str)
parser.add_argument("-o", "--output",
                    dest="out_path",
                    nargs="?",
                    const="./itBins_output.tsv",
                    help="provide file-path to save the modified table of scaffolds to.\n"
                         "If the path is omitted, it is saved as \"itBins_output.tsv\".\n"
                         "If flag is ommited, it is sent to stdout. If you want to\n"
                         "pipe the output, combine with --quiet to supress other\n"
                         "output being sent to stdout, or redirect it to stderr\n"
                         "using --log-to-stderr. If you don't want any output, use\n"
                         "--no-output.\n\n",
                    metavar="OUTPUT_FILE",
                    type=str)
parser.add_argument("-s", "--summary",
                    dest="summary_path",
                    nargs="?",
                    const="./itBins_summary.tsv",
                    help="Provide a file-path to save a summary file to. If the path\n"
                         "is omitted, it is saved as \"itBins_summary.tsv\"\n\n",
                    metavar="SUMMARY_FILE",
                    type=str)
parser.add_argument("-t", "--task-file",
                    dest="task_path",
                    nargs="?",
                    const="./tasks.json",
                    help="set into taskmode, provide path to a task file (.json).\n"
                         "The path may be omitted if \"./tasks.json\" exists. The\n"
                         "taskfile overwrites internal defaults.\n\n",
                    metavar="TASK_FILE",
                    type=str)
parser.add_argument("--example-task-file",
                    dest="task_example_path",
                    nargs="?",
                    const="./tasks.json",
                    help="provide path to create an example task file at, if a\n"
                         "path is omitted it defaults to \"./tasks.json\".\n"
                         "The example task file is a copy of the defaults.\n\n",
                    metavar="TASK_FILE",
                    type=str)
parser.add_argument("-i", "--info", "--verbose",
                    action="store_true",
                    dest="info",
                    help="set flag to put into info mode, will print progress\n"
                         "messages and histograms, overridden by quiet mode.\n\n")
parser.add_argument("-q", "--quiet",
                    action="store_true",
                    dest="quiet",
                    help="\nset flag to put into quiet mode, will not print progress\n"
                         "messages, overrides info mode.\n\n")
parser.add_argument("--log-to-stderr",
                    action="store_true",
                    dest="log_to_stderr",
                    help="\nset flag to redirect all progress and detail printing to\n"
                         "stderr.\n\n")
parser.add_argument("--no-output",
                    action="store_false",
                    dest="no_output",
                    help="\nset flag to supress writing an output file or printing\n"
                         "output to stdout.\n\n")
parser.add_argument("-h", "--help",
                    action="help",
                    help="\nshow this help message and exit.\n\n")
parser.add_argument("--version",
                    action="version",
                    version="%(prog)s 0.6.1",
                    help="\ndisplay version info\n\n")
parser.add_argument("--manual",
                    action="version",
                    version=("╔══════════════════════════════════════════════════════════════════════════════╗\n"
                             "║      ITBINS - MANUAL                                                         ║\n"
                             "╚══════════════════════════════════════════════════════════════════════════════╝\n\n"

                             "█\033[7mHOW THE PRORAM WORKS:\033[0m█\n\n"

                             "load input\n"
                             "prepare input\n"
                             "┌loop over bins\n"
                             "│   ┌loop over tasks\n"
                             "│   │   perform task\n"
                             "prepare output\n"
                             "save\n\n"

                             "█\033[7mUsing a task file\033[0m█\n\n"
                             "The task file allows to set parameters and which tasks are to be performed. To \n"
                             "create this configuration file, use the --example-task-file flag, which will \n"
                             "use the internal defaults.\\nn"
                             "Flags:\n"
                             "These override whatever you provide. Leave them out if you want to provide\n"
                             "them manually.\n\n"

                             "Parameters:\n"
                             "After collecting a list of the bins in the file, the program starts at the index\n"
                             "equal to min_runs and stops at the max_runs index. Setting max_runs to 0 sets it\n"
                             "to the number of bins internally.\n\n"

                             "Tasks:\n"
                             "You can replace \"task_XX\" with any key string you want, as long as they are\n"
                             "unique amoung each other. \"todo\" : \"task_name\" is used to determine which\n"
                             "task should be performed. Tasks are performed in the order given.\n\n"

                             #"---------1---------2---------3---------4---------5---------6---------7---------8\n"
                             "comp_cutoff\n"
                             "If the higher completeness between bacteria and archaea is below the cutoff,\n"
                             "the whole bin is rejected for low completeness. All filters are set to 0 for\n"
                             "each scaffold. Then, if \"break\" is set to True, the task loop breaks and the\n"
                             "algorithm procedes with the post task word and then the next bin. This task\n"
                             "should be applied at the start of the curation of a bin, to catch any bin not\n"
                             "worth working on in the first place and thus save time.\n\n"

                             "tax_prep\n"
                             "This task removes scaffolds based on top level taxonomy. remove_ab_low_comp\n"
                             "removes the lower completeness taxonomy between archaea and bacteria. Sets\n"
                             "Tax_Filter to 0 for relevant scaffolds.\n\n"

                             "GC_pick_expand_limit\n"
                             "This task creates a GC-histogram of the bin, detects maxima and minima and\n"
                             "then selects the largest area minimum-to-minimum peak. Then the peak is\n"
                             "expanded to the baseline on both sides. Finally the peak is limited to a\n"
                             "maximum radius from the maximum of:\n"
                             "    fixed_radius (in GC%) + linear_radius (in GC%) x peak position\n"
                             "Any scaffold with a coverage outside this radius gets its GC_Filter set to 0.\n"
                             "For GC content, the linear radius should be set to 0.\n\n"

                             "Cov_pick_expand_limit\n"
                             "This task creates a Cov-histogram of the bin, detects maxima and minima and\n"
                             "then selects the largest area minimum-to-minimum peak. Then the peak is\n"
                             "expanded to the baseline on both sides. Finally the peak is limited to a\n"
                             "maximum radius from the maximum of:\n"
                             "    fixed_radius + linear_radius x peak position\n"
                             "Any scaffold with a coverage outside this radius gets its Cov_Filter set to 0.\n\n"

                             "tax_drop\n"
                             "This task removes taxonomies, if doing so has a positive effect on the DASTool\n"
                             "binscore, by creating a list of taxonomies at the drop_level (1-6) and for\n"
                             "each taxonomie calculating the difference in DASTool binscore to the bin\n"
                             "without this taxonomy.\n"
                             "It removes the taxonomy with the greatest effect, if multiple\n"
                             "taxonomies are tied, it removes the one with the shortest sum of scaffold\n"
                             "lengths. The taxonomy is only removed, if the sum of the scaffoldlengths is\n"
                             "below the sum of the length of all scaffolds in the bin times max_length. If\n"
                             "drop_if_pos is set to \"true\", any length will be dropped, if the effect is\n"
                             "positive. Taxonomies with no impact will still only be dropped if they fit the\n"
                             "length requirement. Sets Tax_Filter to 0 for relevant scaffolds.\n\n"

                             "metrics\n"
                             "This task calculates a number of metrics and prints them to stdout. If the\n"
                             "summary is created, this task also creates an entry for the bin.\n\n"

                             "append_remainder\n"
                             "This task takes all the scaffolds removed from the bin so far, and assigns\n"
                             "them to a new bin. The bin gets the same name as the produced bin, with the\n"
                             "suffix attached. The new binname is appended to the list of bins and the max\n"
                             "number of runs is incremented by 1. After all bins have been processed, the\n"
                             "remainder bins will be processed using the same list of tasks. If a resonable\n"
                             "completeness cutoff is set as the first task, this represents only a small\n"
                             "increase in total time, as most remainders will be excluded based on low\n"
                             "completeness. \"reframe\" may be set to False, to prevent removing the\n"
                    "remainder scaffolds from the current frame, following tasks will still work on\n" 
                             "the full bin, but any scaffold assigned to the remainder will be worked on\n"
                             "again later, if \"recheck\" is set to True.\n\n"

                             "cont_cutoff\n"
                             "If the contamination of the higher completeness taxonomy between archaea and\n"
                             "bacteria is greater than the cutoff, the whole bin is rejected, setting all\n"
                             "filters to 0. Then, if \"break\" is set to True, the task loop breaks and the\n"
                             "algorithm procedes with the post task word and then the next bin. This task\n"
                             "should be applied only after the curation tasks, to catch bins that weren't\n"
                             "improved sufficiently.\n\n" 

                             "stop\n"
                             "This task breaks out of the task loop and continues with post task work on the\n"
                             "bin. Does not neccessarily jump to the next bin immediately.\n\n"
                             ""),
                    help="\ndisplay manual")
#t_list_of_old_times = []
#t_list_of_new_times = []
if len(sys.argv) < 2:
    parser.parse_args(["-h"])
    sys.exit(1)

#beginning = time.time()

if "-v" in sys.argv:
    print("To avoid ambiguity, the \"-v\" flag is not used.\n"
          "To set verbose, use \"-i\" or \"--info\" or \"--verbose\".\n"
          "To get version , use \"--version\"", file=sys.stderr)
    sys.exit(1)

try:
    cl_args = parser.parse_args()
except Exception:
    traceback.print_exc()

task_dict = {"flags": {"d": "Bin",
                       "u": "curated"},
             "parameters": {"min_runs": 0,
                            "max_runs": 0},
             "tasks": {"task_00": {"todo": "FAST_entry"},
                       "task_01": {"todo": "check_eukarya",
                                   "cutoff": 0.1},
                       "task_02": {"todo": "comp_cutoff",
                                   "cutoff": 0.7,
                                   "break": True},
                       "task_03": {"todo": "tax_prep",
                                   "remove_viruses": True,
                                   "remove_eukarya": True,
                                   "remove_bacteria": False,
                                   "remove_archaea": False,
                                   "remove_unclassified": False,
                                   "remove_other": False,
                                   "remove_ab_low_comp": True},
                       "task_04": {"todo": "peaklist",
                                   "type_of_peaklist": "GC"},
                       "task_05": {"todo": "pick_peak",
                                   "type_of_peaklist": "GC"},
                       "task_06": {"todo": "expand_to_baseline",
                                   "type_of_peaklist": "GC"},
                       "task_07": {"todo": "limit_peak",
                                   "fixed_radius": 7,
                                   "linear_position_radius": 0,
                                   "linear_bin_size_radius": 0,
                                   "type_of_peaklist": "GC"},
                       "task_08": {"todo": "peaklist",
                                   "type_of_peaklist": "Cov"},
                       "task_09": {"todo": "pick_peak",
                                   "type_of_peaklist": "Cov"},
                       "task_10": {"todo": "expand_to_baseline",
                                   "type_of_peaklist": "Cov"},
                       "task_11": {"todo": "limit_peak",
                                   "fixed_radius": 10,
                                   "linear_position_radius": 0.05,
                                   "linear_bin_size_radius": 0,
                                   "type_of_peaklist": "Cov"},
                       "task_12": {"todo": "tax_drop",
                                   "drop_level": 1,
                                   "max_length": 0.00000001,
                                   "drop_if_pos": True},
                       "task_13": {"todo": "FAST_exit"},
                       "task_14": {"todo": "append_remainder",
                                   "suffix": "_remainder",
                                   "reframe": True,
                                   "recheck": True},
                       "task_15": {"todo": "cont_cutoff",
                                   "cutoff": 0.1,
                                   "break": True},
                       "task_16": {"todo": "comp_cutoff",
                                   "cutoff": 0.7,
                                   "break": True},
                       "task_17": {"todo": "metrics"},
                       "task_18": {"todo": "stop"}}}

if cl_args.task_example_path is not None:
    with open(cl_args.task_example_path, "w") as filepath:
        try:
            json.dump(task_dict, filepath, indent='\t') #indent=12
        except Exception:
            print("\n\nWas unable to write task file, will exit", file=sys.stderr)
            traceback.print_exc()
    sys.exit()

curr_dict = copy.deepcopy(task_dict)
if cl_args.task_path is not None:  # in vars(args) :
    with open(cl_args.task_path) as file:
        try:
            new_dict = json.load(file)
            if "flags" in new_dict:
                if "d" in new_dict["flags"]:
                    cl_args.btd_col = new_dict["flags"]["d"]
                if "u" in new_dict["flags"]:
                    cl_args.bub_col = new_dict["flags"]["u"]
                if "c" in new_dict["flags"]:
                    cl_args.compare = new_dict["flags"]["c"]
            if "parameters" in new_dict:
                if "min_runs" in new_dict["parameters"]:
                    curr_dict["parameters"]["min_runs"] = new_dict["parameters"]["min_runs"]
                if "max_runs" in new_dict["parameters"]:
                    curr_dict["parameters"]["max_runs"] = new_dict["parameters"]["max_runs"]
            if "tasks" in new_dict:
                curr_dict["tasks"] = new_dict["tasks"]
        except Exception:
            print("\n\nWas unable to load task file, will exit", file=sys.stderr)
            sys.exit()

loud = 1
info = 0
if cl_args.info:
    info = 1
if cl_args.quiet:
    loud = 0
    info = 0

lprint("\n█\033[7mALGORITHMIC BIN CURATOR\033[0m█\n")

if cl_args.task_path is not None:
    lprint("Working from task file.")


scaffold_identifier = "scaffold"
GC_identifier = "GC"
Cov_identifier = "coverage"
Tax_identifier = "taxonomy"
length_identifier = "length"
Bin_dastool_identifier = "Bin"
Bin_uBin_identifier = "postuBin"

scaffold_dtype = "str"
GC_dtype = "float"
Cov_dtype = "float"
Tax_dtype = "str"
length_dtype = "int"
Bin_dastool_dtype = "str"
Bin_uBin_dtype = "str"


#lt0 = time.time()

try:
    ovw_loaded = pd.read_csv(cl_args.ovw_path, sep="\t")
    lprint("Loaded overview file.")
except Exception :
    print("\n\nWas unable to load overview file, will exit", file = sys.stderr)
    traceback.print_exc()
    sys.exit()

try :
    scg_loaded = pd.read_csv(cl_args.scg_path)
    lprint("Loaded SCG file.")
except Exception :
    print("\n\nWas unable to load SCG file, will exit", file = sys.stderr)
    traceback.print_exc()
    sys.exit()

#lt1 = time.time()

try :
    ovw_loaded[["taxL1", "taxL2", "taxL3", "taxL4", "taxL5", "taxL6"]] = ovw_loaded.taxonomy.str.split(";", n = 5, expand = True)
    #lt2 = time.time()
    ovw_loaded.loc[(ovw_loaded.taxL1 == "") | (ovw_loaded.taxL1 == "None"), "taxL1"] = "unclassified"
    ovw_loaded.loc[(ovw_loaded.taxL2 == "") | (ovw_loaded.taxL2 == "None"), "taxL2"] = "unclassified"
    ovw_loaded.loc[(ovw_loaded.taxL3 == "") | (ovw_loaded.taxL3 == "None"), "taxL3"] = "unclassified"
    ovw_loaded.loc[(ovw_loaded.taxL4 == "") | (ovw_loaded.taxL4 == "None"), "taxL4"] = "unclassified"
    ovw_loaded.loc[(ovw_loaded.taxL5 == "") | (ovw_loaded.taxL5 == "None"), "taxL5"] = "unclassified"
    ovw_loaded.loc[(ovw_loaded.taxL6 == "") | (ovw_loaded.taxL6 == "None"), "taxL6"] = "unclassified"
 
    ovw_loaded.drop(["taxonomy"], axis = 1, inplace = True)
    lprint("Separated tax-strings.")
except Exception :
    print("\n\nWas unable to separate tax-strings, will likely fail tax-tasks", file = sys.stderr)
    traceback.print_exc()
        
try :
    total = pd.merge(left = ovw_loaded, right = scg_loaded, left_on = "scaffold", right_on = "scaffolds")
    total.drop(["scaffolds"], axis = 1, inplace = True)
    lprint("Added SCGs to overview.")
except Exception :
    print("\n\nWas unable to add SCGs to overview, will exit.", file = sys.stderr)
    traceback.print_exc()
    sys.exit()
    

###   spalten für statistik und entscheidung des Algorithmus anfügen    
try :
    total.insert(total.shape[1], "algo", 1)
    total.insert(total.shape[1], "GC_Filter", 1)
    total.insert(total.shape[1], "Cov_Filter", 1)
    total.insert(total.shape[1], "Tax_Filter", 1)
    lprint("Prepared for statistics.")
except Exception :
    print("\n\nWas unable to add statistics columns, will exit.", file = sys.stderr)
    traceback.print_exc()
    sys.exit()

if cl_args.name_prefix != None :
    try :
        total.loc[:, "scaffold"] = cl_args.name_prefix + total.loc[:, "scaffold"]
        lprint("Added prefix to scaffold names.")
    except Exception :
        print("\n\nWas unable to add prefix to scaffold names.", file = sys.stderr)
        traceback.print_exc()
    try :
        total.loc[:, cl_args.bdt_col] = cl_args.name_prefix + total.loc[:, cl_args.bdt_col]
        lprint("Added prefix to DASTool bin names.")
    except Exception :
        print("\n\nWas unable to add prefix to DASTool bin names.", file = sys.stderr)
        traceback.print_exc()
#    if cl_args.compare :
#        try :
#            total.loc[:, cl_args.bub_col] = cl_args.name_prefix + total.loc[:, cl_args.bub_col]
#            lprint("Added prefix to uBin bin names.")
#        except Exception :
#            print("\n\nWas unable to add prefix to uBin bin names.", file = sys.stderr)
#            traceback.print_exc()

try :
    total["bdt"] = total.loc[:, cl_args.bdt_col]
    total["babc"] = total.bdt
    lprint("Found and copied uncurated bin column.")
except Exception as e :
    print(("\n\nWas unable to find uncurated bin column (" + str(cl_args.bdt_col) + "), will exit"), file = sys.stderr)
    traceback.print_exc()
    sys.exit()

try :
    total.babc.fillna("none", inplace = True)
    lprint("Replaced NaNs")
except Exception :
    print("\n\nWas unable to replace NaNs in bin_abc column.", file = sys.stderr)
    traceback.print_exc()
        
#if cl_args.compare :
#    try :
#        total["bub"] = total.loc[:, cl_args.bub_col]
#        lprint("Found curated bin column.")
#    except Exception :
#        print(("\n\nWas unable to find uncurated bin column (" + str(cl_args.bub_col) + "\nProgramm will continue without comparison.\n"), file = sys.stderr)
#        traceback.print_exc()
#        cl_args.compare = False
#if cl_args.compare :
#    try :
#        total.bub.fillna("none", inplace = True)
#        total.insert(total.shape[1], "decision", 0)
#        total.insert(total.shape[1], "evalu", 0)
#        total.insert(total.shape[1], "fdrop", 0)
#        total.insert(total.shape[1], "fkeep", 0)
#        total.insert(total.shape[1], "tdrop", 0)
#        total.insert(total.shape[1], "tkeep", 0)
#        total.decision = (total.bub != "none").astype("int")
#    except Exception :
#        print("\n\nWas unable to set up comparison statistics.\nProgramm will continue without comparison.\n", file = sys.stderr)
#        traceback.print_exc()
#        cl_args.compare = False

minruns = curr_dict["parameters"]["min_runs"]
totalruns = curr_dict["parameters"]["max_runs"]    
    
iprint("totalruns:", totalruns)
    
listOfBinsInFile = total.babc.drop_duplicates()
listOfBinsInFile = listOfBinsInFile.loc[listOfBinsInFile != "none"]
if info :
    iprint("Will work on these bins:")
if totalruns > 0 :
    iprint(listOfBinsInFile.reset_index(drop = True).iloc[minruns : totalruns])
    binList = listOfBinsInFile.reset_index(drop = True).iloc[minruns : totalruns].to_list()
else :
    iprint(listOfBinsInFile.reset_index(drop = True).iloc[minruns : ])
    binList = listOfBinsInFile.reset_index(drop = True).iloc[minruns : ].to_list()
    totalruns = listOfBinsInFile.shape[0]
    
numOfTasks = len(curr_dict["tasks"])
iprint("Num of tasks: ", numOfTasks)


#sys.exit()

###   Iteration über bins
lprint("Started working on bins\n")
start = time.time()
run = 0

if cl_args.summary_path != None :
    bin_summary_list_list = []
    
#if cl_args.compare :
#    zusammenfassung = pd.DataFrame(index = range(0, totalruns), columns = range(10))
#    zusammenfassung.columns = ["bin",
#                               "recall",
#                               "precision",
#                               "F1_score",
#                               "F1_ta",
#                               "binscore",
#                               "binscoreT",
#                               "gcdiffJ",
#                               "covdiffJ",
#                               "comp0"]

################################################################################
### START OF LOOP OVER BINS
################################################################################

#task_beginning = time.time()

for binIndex, listedBin in enumerate(binList) :
    curr_frame = total.loc[total["babc"] == listedBin]
    workingString1 = "\n┌\033[4mWorking on bin " + str(run) + ": \033[3m" + str(listedBin) + ".\033[0m"
    lprint(workingString1)
    bin_summary_list = [listedBin,
                        np.nan,
                        np.nan, 
                        np.nan, 
                        np.nan,
                        np.nan, 
                        np.nan, 
                        np.nan, 
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan,
                        np.nan]
    
    
################################################################################
### START OF LOOP OVER TASKS
################################################################################

    for key in curr_dict["tasks"] :
        lprint("│Performing Task: ", curr_dict["tasks"][key]["todo"])

        
################################################################################
### TASK: FAST_ENTRY
################################################################################

        if curr_dict["tasks"][key]["todo"] == "FAST_entry" :
            #t_0 = time.time()
            ### Enter the FAST lane
            try :
                curr_array = PROT_entry(curr_frame)
            except Exception :
                print("FAST_entry: Was unable to create the dictionary of np.arrays. Will continue with next bin.", file = sys.stderr)
                traceback.print_exc()
                break            
            #task_m1_timer += time.time() - t_0

################################################################################
### TASK: CHECK_EUKARYOTA
################################################################################            
            
        elif curr_dict["tasks"][key]["todo"] == "check_eukarya" :
            #t_0 = time.time()
            print("check euk")
            ### Get parameters from task dictionary
            try :
                cutoff = curr_dict["tasks"][key]["cutoff"]
            except Exception :
                print("Comp_cutoff: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"comp_cutoff\",\n"
                      "           \"cutoff\" : [0 ... 1],\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_0_timer += time.time() - t_0
                continue     
            
            #euk_entries = np.where(curr_array["Tax"][0, :] == "Eukaryota")[0]
            #euk_lengths = curr_array["length"][euk_entries]
            #print(euk_lengths)
            
            euk_length = np.sum(curr_array["length"][np.where(curr_array["Tax"][0, :] == "Eukaryota")[0]])
            total_length = np.sum(curr_array["length"])
            euk_portion = euk_length / total_length
            print("euk_length", euk_length)
            print("total_length", total_length)
            print("euk_portion", euk_portion)

            
            ### Check if the bin might be eukaryotic
            try :
                if euk_portion > cutoff :
                    if cl_args.summary_path != None :
                            bin_summary_list[17] = "suspected"
            except Exception :
                print("CHECK_EUKARYOTA: Unable to determine if bin is eukaryotic. Will continue with next bin.", file = sys.stderr)
                traceback.print_exc()
                break
                
            #task_m1_timer += time.time() - t_0
            
            # sum of lengths where taxL0 is euk
            # np.sum(curr_array["length", np.where(curr_array["Tax"][0, :] == "Eukaryota")])
            # 
            
            
################################################################################
### TASK: COMPLETENESS_CUTOFF
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "comp_cutoff" :
            #t_0 = time.time()
            ### Get parameters from task dictionary
            try :
                cutoff = curr_dict["tasks"][key]["cutoff"]
            except Exception :
                print("Comp_cutoff: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"comp_cutoff\",\n"
                      "           \"cutoff\" : [0 ... 1],\n"
                      "           \"break\" : [bool]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_0_timer += time.time() - t_0
                continue
            
            ### Perform completeness cutoff task
            try :
                curr_array, removed, bacterial_completeness, archaeal_completeness = PROT_completeness_cutoff (curr_array, cutoff)
            except Exception :
                print("Comp_cutoff: Was unable to perform completeness cutoff. Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_0_timer += time.time() - t_0
                continue
            
            ### Give feedback and break if completeness was to low
            higher_completeness = max(bacterial_completeness, archaeal_completeness)
            if removed :
                iprint("Bin removed due to completeness " + str(higher_completeness*100)[:6] + " < " + str(cutoff*100) + " %. Will continue with next bin.")
                if cl_args.summary_path != None :
                    bin_summary_list[16] = "removed due to low completeness"
                #task_0_timer += time.time() - t_0

                ### finish working on the bin if cutoff
                if curr_dict["tasks"][key]["break"] :
                    if cl_args.summary_path != None :
                        bin_summary_list_list.append(bin_summary_list)
                    break
            else :
                iprint("Bin completeness equal to or above threshold of " + str(cutoff*100) + " %.")
                #task_0_timer += time.time() - t_0

################################################################################
### TASK: TAX_PREP
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "tax_prep" :
            #t_1 = time.time()
            
            ### Get parameters from task dictionary
            try :
                viruses_treatment =      curr_dict["tasks"][key]["remove_viruses"]
                eukaryota_treatment =    curr_dict["tasks"][key]["remove_eukarya"]
                bacteria_treatment =     curr_dict["tasks"][key]["remove_bacteria"]
                archaea_treatment =      curr_dict["tasks"][key]["remove_archaea"]
                unclassified_treatment = curr_dict["tasks"][key]["remove_unclassified"]
                other_treatment =        curr_dict["tasks"][key]["remove_other"]
                ab_low_comp_treatment =  curr_dict["tasks"][key]["remove_ab_low_comp"]
            except :
                print("Tax_prep: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"tax_prep\",\n"
                      "           \"remove_viruses\" : [bool],\n"
                      "           \"remove_eukarya\" : [bool],\n"
                      "           \"remove_bacteria\" : [bool],\n"
                      "           \"remove_archaea\" : [bool],\n"
                      "           \"remove_unclassified\" : [bool],\n"
                      "           \"remove_other\" : [bool],\n"
                      "           \"remove_ab_low_comp\" : [bool]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_1_timer += time.time() - t_1
                continue
            
            ### Figure out what to remove between archaea and bacteria
            if ab_low_comp_treatment :
                bacterial_completeness = FAST_completeness(curr_array["bacterial_SCGs"])
                archaeal_completeness =  FAST_completeness(curr_array["archaeal_SCGs"])
                bacteria_treatment, archaea_treatment = (False, True) if bacterial_completeness > archaeal_completeness else (True, False)
            
            ### Remove taxonomies
            try :
                curr_array = PROT_remove_taxonomy(curr_array,
                                                  remove_viruses =      viruses_treatment,
                                                  remove_eukaryota =    eukaryota_treatment,
                                                  remove_bacteria =     bacteria_treatment,
                                                  remove_archaea =      archaea_treatment,
                                                  remove_unclassified = unclassified_treatment,
                                                  remove_other =        other_treatment)
                if viruses_treatment :      iprint("Viruses removed")
                if eukaryota_treatment :    iprint("Eukaryota removed")
                if bacteria_treatment :     iprint("Bacteria removed")
                if archaea_treatment :      iprint("Archaea removed")
                if unclassified_treatment : iprint("Unclassifieds removed")
                if other_treatment :        iprint("Others removed")
                #task_1_timer += time.time() - t_1
                
            except :
                print("Tax_prep: error while trying to remove taxonomies. Will continue with next task.")
                traceback.print_exc()
                #task_1_timer += time.time() - t_1
                continue
            
                
################################################################################
### TASK: PEAKLIST
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "peaklist" :
            #t_2 = time.time()
            
            ### Get parameters from task dictionary
            try :
                type_of_peaklist = curr_dict["tasks"][key]["type_of_peaklist"]
            except Exception :
                print("Peaklist: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"peaklist\",\n"
                      "           \"type_of_peaklist\" : [\"GC\", \"Cov\"]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_2_timer += time.time() - t_2
                continue
            
            if type_of_peaklist == "GC" :
                
                ### Create GC peaklist
                try :
                    peaklist_GC, histogram_GC = PROT_create_peaklist(curr_array, type_of_peaklist = "GC")
                    #task_2_timer += time.time() - t_2
                except Exception :
                    print("Peaklist: Was unable to create peaklist. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_2_timer += time.time() - t_2
                    continue
                
            elif type_of_peaklist == "Cov" :
                
                ### Create Cov peaklist
                try :
                    peaklist_Cov, histogram_Cov = PROT_create_peaklist(curr_array, type_of_peaklist = "Cov")
                    #task_2_timer += time.time() - t_2
                except Exception :
                    print("Peaklist: Was unable to create peaklist. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_2_timer += time.time() - t_2
                    continue
                
            else :
                
                ### Illegal argument
                exception_message = "Peaklist: Illegal argument \"" + str(type_of_peaklist) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\". Will continue with next task."
                raise Exception(exception_message)
                traceback.print_exc()
                #task_2_timer += time.time() - t_2
                continue
                

################################################################################
### TASK: PICK_PEAK
################################################################################
        
        elif curr_dict["tasks"][key]["todo"] == "pick_peak" :
            #t_3 = time.time()
            
            ### Get parameters from task dictionary
            try :
                type_of_peaklist = curr_dict["tasks"][key]["type_of_peaklist"]
            except Exception :
                print("Pick_peak: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"pick_peak\",\n"
                      "           \"type_of_peaklist\" : [\"GC\", \"Cov\"]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_3_timer += time.time() - t_3
                continue
            
            if type_of_peaklist == "GC" :
                
                ### Pick GC peak
                try :
                    peaklist_GC, histogram_GC = PROT_pick_peak(peaklist_GC, histogram_GC, type_of_peaklist = "GC")
                    #task_3_timer += time.time() - t_3
                except Exception :
                    print("Pick_peak: Was unable to pick peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_3_timer += time.time() - t_3
                    continue
                
            elif type_of_peaklist == "Cov" :
                
                ### Pick Cov peak
                try :
                    peaklist_Cov, histogram_Cov = PROT_pick_peak(peaklist_Cov, histogram_Cov, type_of_peaklist = "Cov")
                    #task_3_timer += time.time() - t_3
                except Exception :
                    print("Pick_peak: Was unable to pick peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_3_timer += time.time() - t_3
                    continue
                
            else :
                
                ### Illegal argument
                exception_message = "Pick_peak: Illegal argument \"" + str(type_of_peaklist) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\". Will continue with next task."
                raise Exception(exception_message)
                traceback.print_exc()
                #task_3_timer += time.time() - t_3
                continue


################################################################################
### TASK: EXPAND_TO_BASELINE
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "expand_to_baseline" :
            #t_4 = time.time()
            
            ### Get parameters from task dictionary
            try :
                type_of_peaklist = curr_dict["tasks"][key]["type_of_peaklist"]
            except Exception :
                print("Expand_to_baseline: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"pick_peak\",\n"
                      "           \"type_of_peaklist\" : [\"GC\", \"Cov\"]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_4_timer += time.time() - t_4
                continue
            
            if type_of_peaklist == "GC" :
                
                ### Expand GC peak
                try :
                    peaklist_GC, histogram_GC = PROT_expand_to_baseline(peaklist_GC, histogram_GC, type_of_peaklist = "GC")
                    #task_4_timer += time.time() - t_4
                except Exception :
                    print("Expand_to_baseline: Was unable to expand peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_4_timer += time.time() - t_4
                    continue
                
            elif type_of_peaklist == "Cov" :
                
                ### Expand Cov peak
                try :
                    peaklist_Cov, histogram_Cov = PROT_expand_to_baseline(peaklist_Cov, histogram_Cov, type_of_peaklist = "Cov")
                    #task_4_timer += time.time() - t_4
                except Exception :
                    print("Expand_to_baseline: Was unable to expand peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_4_timer += time.time() - t_4
                    continue
                
            else :
                
                ### Illegal argument
                exception_message = "Expand_to_baseline: Illegal argument \"" + str(type_of_peaklist) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\". Will continue with next task."
                raise Exception(exception_message)
                traceback.print_exc()
                #task_4_timer += time.time() - t_4
                continue
                
                
################################################################################
### TASK: LIMIT_PEAK
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "limit_peak" :
            #t_5 = time.time()
            
            ### Get parameters from task dictionary
            try :
                fixed_radius =           curr_dict["tasks"][key]["fixed_radius"]
                linear_position_radius = curr_dict["tasks"][key]["linear_position_radius"]
                linear_bin_size_radius = curr_dict["tasks"][key]["linear_bin_size_radius"]
                type_of_peaklist =       curr_dict["tasks"][key]["type_of_peaklist"]
            except Exception :
                print("Limit_peak: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"limit_peak\",\n"
                      "           \"fixed_radius\" : float,\n"
                      "           \"linear_position_radius\" : float[0 ... 1],\n"
                      "           \"linear_bin_size_radius\" : float[0 ... 1],\n"
                      "           \"type_of_peaklist\" : [\"GC\", \"Cov\"]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_5_timer += time.time() - t_5
                continue
            
            if type_of_peaklist == "GC" :
                
                ### Limit GC peak
                try :
                    curr_array = PROT_limit_peak(curr_array, peaklist_GC, histogram_GC, fixed_radius, linear_position_radius, linear_bin_size_radius, filter_type = "GC")
                    #task_5_timer += time.time() - t_5
                except Exception :
                    print("Expand_to_baseline: Was unable to expand peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_5_timer += time.time() - t_5
                    continue
                
            elif type_of_peaklist == "Cov" :
                
                ### Limit Cov peak
                try :
                    curr_array = PROT_limit_peak(curr_array, peaklist_Cov, histogram_Cov, fixed_radius, linear_position_radius, linear_bin_size_radius, filter_type = "Cov")
                    #task_5_timer += time.time() - t_5
                except Exception :
                    print("Limit_peak: Was unable to expand peak. Will continue with next task.", file = sys.stderr)
                    traceback.print_exc()
                    #task_5_timer += time.time() - t_5
                    continue
                
            else :
                
                ### Illegal argument
                exception_message = "Limit_peak: Illegal argument \"" + str(type_of_peaklist) + "\" passed for type_of_peaklist. Argument may be \"GC\" or \"Cov\". Will continue with next task."
                raise Exception(exception_message)
                traceback.print_exc()
                #task_5_timer += time.time() - t_5
                continue

################################################################################
### TASK: TAX_DROP
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "tax_drop" :
            #t_6 = time.time()
            #print("before getting params", curr_array["Tax"][1])
            ### Get parameters from task dictionary
            try :
                drop_level = curr_dict["tasks"][key]["drop_level"]
                max_length = curr_dict["tasks"][key]["max_length"]
                drop_if_pos = curr_dict["tasks"][key]["drop_if_pos"]
            except Exception :
                print("Limit_peak: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"limit_peak\",\n"
                      "           \"drop_level\" : int[0 ... 5],\n"
                      "           \"max_length\" : float[0 ... 1],\n"
                      "           \"drop_if_pos\" : bool,\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_6_timer += time.time() - t_6
                continue

            ### Loop over removing taxa
            did_something = True
            while did_something :
                try :
                    #print("before function call", curr_array["Tax"][drop_level])
                    curr_array, did_something = PROT_tax_drop(curr_array, drop_level = drop_level, max_length = max_length, always_drop_if_score_improves = drop_if_pos)
                except Exception :
                        print("Tax_drop: Was unable to drop taxonomy.", file = sys.stderr)
                        traceback.print_exc()
                        break
            #tax_table, bacterial_score, archaeal_score, useless = PROT_tax_table(curr_array)
            #print(pd.DataFrame(data = tax_table))
            #task_6_timer += time.time() - t_6


################################################################################
### TASK: PROT_METRICS
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "metrics" :
            #print('I love beer')      #this stays in (contributed by jonas hülsermann on fr, 20230707)
            #t_7 = time.time()
            if True :
                #algo = curr_array["filters"][0, :] & curr_array["filters"][1, :] & curr_array["filters"][2, :]
                
                algo_bacterial_completeness = FAST_completeness(FAST_filter(curr_array["bacterial_SCGs"], curr_array["filters"]))
                algo_archaeal_completeness = FAST_completeness(FAST_filter(curr_array["archaeal_SCGs"], curr_array["filters"]))
                algo_bacterial_contamination = FAST_contamination(FAST_filter(curr_array["bacterial_SCGs"], curr_array["filters"]))
                algo_archaeal_contamination = FAST_contamination(FAST_filter(curr_array["archaeal_SCGs"], curr_array["filters"]))
                algo_bacterial_score = FAST_dt_score(FAST_filter(curr_array["bacterial_SCGs"], curr_array["filters"]))
                algo_archaeal_score = FAST_dt_score(FAST_filter(curr_array["archaeal_SCGs"], curr_array["filters"]))
                algo_binsize = np.shape(curr_array["filters"][0][(curr_array["filters"][0] & curr_array["filters"][1] & curr_array["filters"][2]) == 1])[-1]
                #print(curr_array["filters"])
                #print(curr_array["GC"][(curr_array["filters"][0] & curr_array["filters"][1] & curr_array["filters"][2]) == 1])
                #original_array = PROT_entry(curr_frame)
                original_bacterial_completeness = FAST_completeness(curr_array["bacterial_SCGs"])
                original_archaeal_completeness = FAST_completeness(curr_array["archaeal_SCGs"])
                original_bacterial_contamination = FAST_contamination(curr_array["bacterial_SCGs"])
                original_archaeal_contamination = FAST_contamination(curr_array["archaeal_SCGs"])
                original_bacterial_score = FAST_dt_score(curr_array["bacterial_SCGs"])
                original_archaeal_score = FAST_dt_score(curr_array["archaeal_SCGs"])
                original_binsize = np.shape(curr_array["filters"])[-1]
                
                if cl_args.summary_path != None :
                    bin_summary_list[1:15] = [original_binsize,
                                              original_bacterial_completeness, 
                                        original_bacterial_contamination, 
                                        original_bacterial_score,
                                        original_archaeal_completeness, 
                                        original_archaeal_contamination, 
                                        original_archaeal_score, 
                                        algo_binsize,
                                        algo_bacterial_completeness,
                                        algo_bacterial_contamination,
                                        algo_bacterial_score,
                                        algo_archaeal_completeness,
                                        algo_archaeal_contamination,
                                        algo_archaeal_score,
                                        (max(algo_archaeal_score,algo_bacterial_score)-max(original_archaeal_score,original_bacterial_score))]
                    #bin_summary_list_list.append(bin_summary_list)
                

            lprint("│\n│Bin metrics:")
            lprint("│Size of bin: " + str(original_binsize))
            lprint("│Fraction kept by algorithm after filters: " + str(algo_binsize/original_binsize))
            lprint("│Bacterial completeness and contamination: " + str(original_bacterial_completeness) + ", " + str(original_bacterial_contamination))
            lprint("│Archaeal completeness and contamination: " + str(original_archaeal_completeness) + ", " + str(original_archaeal_contamination))
            lprint("│Bacterial and archaeal binscore: " + str(original_bacterial_score) + ", " + str(original_archaeal_score))
            lprint("│Bacterial completeness and contamination (algorithm): " + str(algo_bacterial_completeness) + ", " + str(algo_bacterial_contamination))
            lprint("│Archaeal completeness and contamination (algorithm): " + str(algo_archaeal_completeness) + ", " + str(algo_archaeal_contamination))
            lprint("│Bacterial and archaeal binscore (algorithm): " + str(algo_bacterial_score) + ", " + str(algo_archaeal_score))
            #lprint("│GC limits: " + str(gclower) + " : " + str(gcupper) + " (" + str(gcdiffJ) + ")")
            #lprint("│Cov limits: " + str(covlower) + " : " + str(covupper) + " (" + str(covdiffJ) + ")")

#            if cl_args.compare :
#                lprint(("│Comparison of the algorithm to the human:"))
#                lprint(("│\ttdrop: " + str(al_tn)))
#                lprint(("│\ttkeep: " + str(al_tp)))
#                lprint(("│\tfdrop: " + str(al_fn)))
#                lprint(("│\tfkeep: " + str(al_fp)))
#                lprint(("│\tPrecision: " + str(al_prec)))
#                lprint(("│\tRecall: " + str(al_reca)))
#                lprint(("│\tF1: " + str(al_f1)))
#                lprint(("│\tF1 (trivial acceptor): " + str(ta_f1)))
#                lprint(("│\tBacterial completeness and contamination (human): " + str(compBacT) + ", " + str(contBacT)))
#                lprint(("│\tArchaeal completeness and contamination (human): " + str(compArcT) + ", " + str(contArcT)))
#                lprint(("│\tBacterial and archaeal binscore (human): " + str(ScoreBacT) + ", " + str(ScoreArcT)))
#                lprint(("│\tGC limits (human): " + str(gclowertill) + " : " + str(gcuppertill) + " (" + str(gcdiffT) + ")"))
#                lprint(("│\tCov limits (human): " + str(covlowertill) + " : " + str(covuppertill) + " (" + str(covdiffT) + ")"))
            
            lprint("│")
            #task_7_timer += time.time() - t_7
            
            
################################################################################
### TASK: FAST_EXIT
################################################################################

        elif curr_dict["tasks"][key]["todo"] == "FAST_exit" :
            #t_71 = time.time()
            curr_frame.loc[:, "GC_Filter"] = pd.Series(data = curr_array["filters"][0], dtype = "int", index = curr_frame.index)
            curr_frame.loc[:, "Cov_Filter"] = pd.Series(data = curr_array["filters"][1], dtype = "int", index = curr_frame.index)
            curr_frame.loc[:, "Tax_Filter"] = pd.Series(data = curr_array["filters"][2], dtype = "int", index = curr_frame.index)
            #task_71_timer += time.time() - t_71
            print("\nCurrarray:\n", curr_array["filters"][0])
            print("\nCurrframe:\n", curr_frame.loc[:, "GC_Filter"])

            
################################################################################
### TASK: APPEND_REMAINDER
################################################################################
#jump
        elif curr_dict["tasks"][key]["todo"] == "append_remainder" :
            #t_8 = time.time()
            
            curr_size = curr_frame.shape[0]
            remainder_size = curr_frame.loc[~(curr_frame.Tax_Filter & curr_frame.GC_Filter & curr_frame.Cov_Filter).astype("bool"), "babc"].shape[0]
            if curr_size == remainder_size :
                print("\n\nWhole bin has been dropped, no remainder will be appended\n\n")
                continue
            
            
            try :
                curr_frame.loc[~(curr_frame.Tax_Filter & curr_frame.GC_Filter & curr_frame.Cov_Filter).astype("bool"), "babc"] = (listedBin + curr_dict["tasks"][key]["suffix"])
                #print("binnames:", curr_frame.babc.drop_duplicates())
                #total.loc[total["babc"] == listedBin, "babc"] = curr_frame.babc
                
                #curr_frame = total.loc[total["babc"] == listedBin]
                
                iprint("Assigning dropped scaffolds to new bin")
                #print("remainder_shape : ", curr_frame[curr_frame.babc == (listedBin + curr_dict["tasks"][key]["suffix"])].shape[0])
                iprint((curr_frame["babc"].drop_duplicates()))
                remainder_size = curr_frame.loc[~(curr_frame.Tax_Filter & curr_frame.GC_Filter & curr_frame.Cov_Filter).astype("bool"), "babc"].shape[0]
                iprint("remainder: ", remainder_size)
            except Exception :
                print("\n\nWhile append_remainder task: Was unable to update bin names\n", file = sys.stderr)
                traceback.print_exc()
            if curr_dict["tasks"][key]["recheck"] :
                if remainder_size > 0 :
                    try :
                        binList.insert(binIndex + 1, str(listedBin + curr_dict["tasks"][key]["suffix"]))
                        iprint("Appending remainder bin to bin list")
                    except Exception :
                        print("\n\nWhile append_remainder task: Was unable to append remainder to binList\n", file = sys.stderr)
                        traceback.print_exc()
                    totalruns += 1
                else :
                    iprint("No remainder")
                    
            iprint("Updating total frame")   
            total.loc[total["babc"] == listedBin, "babc"] = curr_frame.babc
            
            if curr_dict["tasks"][key]["reframe"] :
                try :
                    curr_frame = curr_frame.loc[curr_frame.babc == listedBin, :]
                    iprint("Removing remainder scaffolds from current frame")
                except Exception :
                    print("\n\nWhile append_remainder task: Was unable to update current frame\n", file = sys.stderr)
                    traceback.print_exc()
            #task_8_timer += time.time() - t_8
            #print("binnames:", curr_frame.babc.drop_duplicates())
            
            
################################################################################
### TASK: CONT_CUTOFF
################################################################################
        
        elif curr_dict["tasks"][key]["todo"] == "cont_cutoff" :
            #t_9 = time.time()
            
            ### Get parameters from task dictionary
            try :
                cutoff = curr_dict["tasks"][key]["cutoff"]
            except Exception :
                print("Cont_cutoff: error while trying to read task parameters. Make sure task is of form:\n"
                      "\"[key]\" : {\"todo\" : \"cont_cutoff\",\n"
                      "           \"cutoff\" : [0 ... 1],\n"
                      "           \"break\" : [bool]},\n"
                      "Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_9_timer += time.time() - t_9
                continue
            
            ### Perform completeness cutoff task
            try :
                curr_array, removed, relevant_contamination, contamination_type = PROT_contamination_cutoff(curr_array, cutoff)
                print("\n@@@@@\n", relevant_contamination, removed)
            except Exception :
                print("Cont_cutoff: Was unable to perform contamination cutoff. Will continue with next task.", file = sys.stderr)
                traceback.print_exc()
                #task_9_timer += time.time() - t_9
                continue
            
            ### Give feedback and break if completeness was to low
            higher_completeness = max(bacterial_completeness, archaeal_completeness)
            if removed :
                iprint("Bin removed due to contamination " + str(higher_completeness*100)[:6] + " > " + str(cutoff*100) + " %. Will continue with next bin.")
                #task_9_timer += time.time() - t_9

                ### finish working on the bin if cutoff
                if curr_dict["tasks"][key]["break"] :
                    if cl_args.summary_path != None :
                        bin_summary_list[16] = "removed due to high contamination"
                        #bin_summary_list_list.append(bin_summary_list)
                    break
            else :
                iprint("Bin completeness below or equal to threshold of " + str(cutoff*100) + " %.")
                #task_9_timer += time.time() - t_9


################################################################################
### TASK: STOP
################################################################################
        
        elif curr_dict["tasks"][key]["todo"] == "stop" :
            break
            
################################################################################
### UNKNOWN OR MISSSPELLED TASK
################################################################################

        else :
            print(("Unknown or missspelled task: \"" + str(curr_dict["tasks"][key]["todo"]) + "\". If using a task file, make sure that all tasks contain instructions at \"todo\" : \"instruction\", where instruction is the name of a task:\n\tcomp_cutoff\n\ttax_prep\n\tGC_pick_expand_limit\n\tCov_pick_expand_limit\n\ttax_drop\n\tmetrics\n\tappend_remainder\n\tstop\n\nWill continue with next task."), file = sys.stderr)
            raise Exception("well shit")
            traceback.print_exc()
            continue
            
################################################################################
### ACTUAL POST TASKS SECTION / END OF TASK LOOP
################################################################################

    if cl_args.summary_path != None :
        bin_summary_list_list.append(bin_summary_list)

    if cl_args.out_path != None :
        try :
            total.loc[total["babc"] == listedBin, "GC_Filter"] = curr_frame.GC_Filter
            total.loc[total["babc"] == listedBin, "Cov_Filter"] = curr_frame.Cov_Filter
            total.loc[total["babc"] == listedBin, "Tax_Filter"] = curr_frame.Tax_Filter
            total.loc[total["babc"] == listedBin, "babc"] = curr_frame.babc
            lprint("│Updating results.")
        except Exception :
            print("While updating results: failed to update total frame.", file = sys.stderr)
            traceback.print_exc()
            

    fertig = math.floor(time.time() - start)
    workingString2 = "│Bin " + str(run) + ": " + str(listedBin) + " processed after " + str(fertig) + " seconds.\n"
    workingStringLength2 = len(workingString2) -2
    lprint(workingString2 + "└" + "─" * workingStringLength2)
    run += 1

################################################################################
### POST BINS SECTION / END OF BIN LOOP
################################################################################

#task_done = time.time()

### Modify output
#~(fa & fb & fc).astype("bool")
#print(total.loc[(total.GC_Filter == 1) & (total.Cov_Filter == 0), ["GC_Filter", "Cov_Filter", "Tax_Filter"]])
#print(total.loc[87571:87580, ["scaffold", "bdt", "GC_Filter", "Cov_Filter", "Tax_Filter"]])
#print(total.loc[total.babc.isna(), ["scaffold", "bdt", "GC_Filter", "Cov_Filter", "Tax_Filter"]])


total.loc[~(total.GC_Filter & total.Cov_Filter & total.Tax_Filter).astype("bool"), "algo"] = 0
total.loc[total.algo == 0, "babc"] = "none"

if False:#cl_args.compare :
    total_out = total.loc[:,["scaffold", cl_args.bdt_col, cl_args.bub_col, "babc"]] #put inj flag
    total_out.columns = ["scaffold", "bin_DASTool", "bin_uBin", "bin_Algorithm"]
else:
    total_out = total.loc[:,["scaffold", cl_args.bdt_col, "babc"]]
    total_out.columns = ["scaffold", "bin_DASTool", "bin_Algorithm"]

total_out.fillna(value = "none", inplace = True)

################################################################################
### SAVE OUTPUT OR PRINT TO STDOUT
################################################################################

if cl_args.no_output : # setting flag sets False
    if cl_args.out_path != None :
        try :
            total_out.to_csv(cl_args.out_path, sep = "\t", index = False)
            lprint("\nSaving result to: " + str(cl_args.out_path))
        except Exception :
            print("\nWhile saving result: Was unable to save.\n", file = sys.stderr)
            traceback.print_exc()
    else :
        try :
            total_out.to_csv(sys.stdout, sep = "\t", index = False)
            lprint("\nSent result to: stdout")
        except Exception :
            print("\nWhile sending result to stdout: Was unable to print.\n", file = sys.stderr)
            traceback.print_exc()
else :
    lprint("\nResult discarded")

if cl_args.summary_path != None :
    print(bin_summary_list_list)
    print(np.array(bin_summary_list_list).shape)
    summary_frame = pd.DataFrame(data = bin_summary_list_list)
    summary_frame = summary_frame.iloc[:, 0:18]
    summary_frame.columns = ["bin",                                       
                             "#scaffolds(dt)", 
                             "bac_completeness(dt)", 
                             "bac_contamination(dt)", 
                             "bac_score(dt)", 
                             "arc_completeness(dt)", 
                             "arc_contamination(dt)", 
                             "arc_score(dt)",
                             "#scaffolds(abc)", 
                             "bac_completeness(abc)", 
                             "bac_contamination(abc)", 
                             "bac_score(abc)", 
                             "arc_completeness(dabct)", 
                             "arc_contamination(abc)", 
                             "arc_score(abc)", 
                             "score_improvement",
                             "comp/cont_flag",
                             "Eukaryote?"]
    print(summary_frame)
    summary_frame.loc[(summary_frame.score_improvement <= -0.1), "flag"] = "!RECHECK: drop in dt bin score"
    summary_frame.loc[(summary_frame.flag == ""), "flag"] = "fine"
    summary_frame.to_csv(cl_args.summary_path, sep = "\t", index = False)
    
################################################################################
### END OF PROGRAM
################################################################################

lprint("\n\033[7mALL DONE!\033[0m")       

