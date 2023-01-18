#!/usr/bin/env python3

import re
import json
import pandas as pd
from typing import List

def json_to_dataframes(json_file: str) -> List[pd.DataFrame]:

    with open(json_file, "r", encoding = "utf-8") as f:
        data = json.load(f)
        f.close()

    fragment = {'frag_seq':[],
                'frag_type1':[],
                'frag_type2':[],
                'position_frag_type1':[],
                'position_frag_type2':[],
                'frag_length':[],
                'frag_charge':[],
                'frag_intensity':[],
                'frag_mz':[],
                'perc_of_total_intensity':[],
                'prop_intensity_to_base_peak':[],
                'modification':[],
                'scan_number':[],
                'ambiguity':[]}

    spectrum = {'perc_internal':[],
                'perc_terminal':[],
                'perc_other':[],
                'total_int_internal':[],
                'total_int_terminal':[],
                'total_int_other':[],
                'top1_internal_ion_code':[],
                'top1_internal_seq':[],
                'top2_internal_ion_code':[],
                'top2_internal_seq':[],
                'top3_internal_ion_code':[],
                'top3_internal_seq':[],
                'intensity_explained_aa':[],
                'intensity_explained_precursor':[],
                'scan_number':[],
                'score':[],
                'peptide_seq':[],
                'peptide_length':[]}

    for entry in data:
        pep_seq = entry["sequence"]
        pep_len = len(pep_seq)

        for i, code in enumerate(entry["annotation"]["theoretical_code"]):
            # Fragment

            # skip non-annonated fragments
            if code is None:
                continue

            start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(code)
            frag_seq = pep_seq[start-1:end]
            fragment['frag_seq'].append(frag_seq)
            fragment['frag_type1'].append(ion_cap_start)
            fragment['frag_type2'].append(ion_cap_end)
            fragment['position_frag_type1'].append(start)
            fragment['position_frag_type2'].append(end)
            fragment['frag_length'].append(end-start+1)
            fragment['frag_charge'].append(charge)
            fragment['frag_intensity'].append(entry["annotation"]["intensity"][i])
            fragment['frag_mz'].append(entry["annotation"]["mz"][i])
            total_intensity = sum(entry["annotation"]["intensity"])
            fragment['perc_of_total_intensity'].append(entry["annotation"]["intensity"][i]/total_intensity)
            intensity_base_peak = max(entry["annotation"]["intensity"])
            fragment['prop_intensity_to_base_peak'].append(entry["annotation"]["intensity"][i]/intensity_base_peak)
            fragment['modification'].append(parse_modfication(entry["proforma"],start,end))
            fragment['scan_number'].append(0) # Change later
            fragment['ambiguity'].append(0) # Change later

        # Spectrum
        percentages_and_total_intensities = calculate_internal_terminal_non_annotated_ions(entry["annotation"]["theoretical_code"],entry["annotation"]["intensity"])
        spectrum['perc_internal'].append(percentages_and_total_intensities['internal'])
        spectrum['perc_terminal'].append(percentages_and_total_intensities['terminal'])
        spectrum['perc_other'].append(percentages_and_total_intensities['non_annotated'])
        spectrum['total_int_internal'].append(percentages_and_total_intensities['total_int_internal'])
        spectrum['total_int_terminal'].append(percentages_and_total_intensities['total_int_terminal'])
        spectrum['total_int_other'].append(percentages_and_total_intensities['total_int_non'])
        top_3 = find_top3_most_intense_internal_ions(entry["annotation"]["theoretical_code"],entry["annotation"]["intensity"],pep_seq)
        spectrum['top1_internal_ion_code'].append(top_3[0][0])
        spectrum['top1_internal_seq'].append(top_3[1][0])
        spectrum['top2_internal_ion_code'].append(top_3[0][1])
        spectrum['top2_internal_seq'].append(top_3[1][1])
        spectrum['top3_internal_ion_code'].append(top_3[0][2])
        spectrum['top3_internal_seq'].append(top_3[1][2])
        spectrum['intensity_explained_aa'].append(0)
        spectrum['intensity_explained_precursor'].append(0)
        spectrum['scan_number'].append(0)
        spectrum['score'].append(0)
        spectrum['peptide_seq'].append(pep_seq)
        spectrum['peptide_length'].append(len(pep_seq))

    return [pd.DataFrame(fragment),pd.DataFrame(spectrum)]

def find_top3_most_intense_internal_ions(fragments,intensities,pep_seq):
    mapping = dict()
    for i, intensity in enumerate(intensities):
        # Add only the internal ions
        if fragments[i] is not None and "t" not in fragments[i]:
            mapping[intensity] = fragments[i]

    top_3_intensities = sorted(list(mapping.keys()), reverse = True)[:3]
    top_3_codes = [mapping[intensity] for intensity in top_3_intensities]

    top_3_sequences = []

    for code in top_3_codes:
        start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(code)
        top_3_sequences.append(pep_seq[start-1:end])

    return [top_3_codes,top_3_sequences]

def calculate_internal_terminal_non_annotated_ions(fragments,intensities):
    # Number of different ion types
    internal = 0
    terminal = 0
    non_annotated = 0

    # Total intensities of ion types
    total_int_internal = 0
    total_int_terminal = 0
    total_int_non = 0

    for i,fragment in enumerate(fragments):
        if fragment is not None:
            if 't' in fragment:
                terminal += 1
                total_int_terminal += intensities[i]
            else:
                internal += 1
                total_int_internal += intensities[i]
        else:
            non_annotated += 1
            total_int_non += intensities[i]

    return {'internal':internal/len(fragments),
            'terminal':terminal/len(fragments),
            'non_annotated':non_annotated/len(fragments),
            'total_int_internal':total_int_internal,
            'total_int_terminal':total_int_terminal,
            'total_int_non':total_int_non}

def parse_modfication(proforma, start, end):

    # Finds the modifications in the proforma string
    pattern = r"\[([A-Za-z0-9_]+)\]"
    mods = re.findall(pattern, proforma)

    positions = parse_modification_positions(proforma)

    modifications = ""

    # Finds the positions of the modifications
    for i, position in enumerate(positions):
        if position >= start and position <= end:
            modifications += mods[i]

    return modifications

def parse_modification_positions(seq):

    i = 1
    is_modification = False
    annotation = []
    for aa in seq:
        if aa == "[":
            is_modification = True
            annotation.append(0)
        elif aa == "]":
            is_modification = False
        else:
            if not is_modification:
                annotation.append(i)
                i += 1

    positions = []

    for i, anno in enumerate(annotation):
        if anno == 0:
            positions.append(annotation[i-1])

    return positions

def parse_fragment_code(fragment_code: str):

    # test if fragment code format is valid*
    fragment_code_pattern = re.compile(".+(:).+(@)[0-9]+(:)[0-9]+(\()(\+|\-)[0-9](\))(\[(.*?)\])?")
    if bool(fragment_code_pattern.match(fragment_code)) == False:
        raise RuntimeError("Incorrect fragment code format: {0}".format(fragment_code))

    ## Parse fragment code

    start, end = [
        int(i) for i in re.search("(?<=\@)(.*?)(?=\()", fragment_code).group(1).split(":")
    ]  # Get start and end amino acid indexes
    ion_cap_start, ion_cap_end = [
        str(i) for i in re.search("^(.*?)(?=\@)", fragment_code).group(1).split(":")
    ]  # Get start and end ion caps name
    charge = int(re.search("(?<=\()(.*?)(?=\))", fragment_code).group(1))  # get charge state
    formula = re.search("(?<=\[)(.*?)(?=\])", fragment_code)
    if formula == None:
        formula = ""
    else:
        formula = str(re.search("(?<=\[)(.*?)(?=\])", fragment_code).group(1))

    return start, end, ion_cap_start, ion_cap_end, charge, formula