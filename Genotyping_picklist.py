# James Jensen 2019/2020
# Modified by Chris J, 5/2020
# Last modification 12/2020

# prior to use "credentials.json" needs to be put in the same folder as the script
# to do this, visit:
# https://developers.google.com/sheets/api/quickstart/python
# follow the directions for the first 2 steps, download the credentials.json file and
#  move it to the same directory as this script

# usage examples

# minimum info, includes excel file and sheet of file with data. output files are named automatically
# python platemap2genotyping_setup.py CPS13-B.source_plate_info.xlsx source_plate1

# similar to above, but a custom prefix is used to name output files
# python platemap2genotyping_setup.py CPS13-B.source_plate_info.xlsx source_plate1 --output_prefix crazy_ed

# Similar to first example, but iterates through spefically named sheets in a file
# for i in {1..4}; do python platemap2genotyping_setup.py CPS13-B.source_plate_info.xlsx source_plate${i}; done

# --mix_plates -restrict potential primer mix source plates to these, in format: '10','14','7'
# --seq_plates -restrict potential sequencing primer mix source plates to these, in format: '10','14','7'

import argparse
import warnings
import pickle
import os.path
import sys
from string import ascii_uppercase
import itertools
import traceback
import re

import pandas as pd
import numpy as np

from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request


def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    # traceback.print_stack()
    return str(msg) + '\n'


warnings.formatwarning = custom_formatwarning

# If modifying these scopes, delete the file token.pickle.
SCOPES = ['https://www.googleapis.com/auth/spreadsheets.readonly']


def get_creds():
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)#/Users/chris.johnson/.config/gcloud/application_default_credentials.json
            creds = flow.run_local_server()
        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
    return creds


def fetch_data(spreadsheet_id, sheet, data_range):
    creds = get_creds()
    service = build('sheets', 'v4', credentials=creds)
    sheet_range_str = f'{sheet}!{data_range}'
    sheet = service.spreadsheets()
    result = sheet.values().get(spreadsheetId=spreadsheet_id,
                                range=sheet_range_str).execute()
    values = result.get('values', [])
    return values


#PRIMER_MASTERLIST_SPREADSHEET_ID = '1PP3Nv1u8jhiFCJXvcTdZ2YCoJSmLsJhbimS7Sx2exIo'
#PRIMERS_SHEETNAME = 'Primers'
#PRIMER_MIXES_SHEETNAME = 'Primer Mixes'
PRIMER_MASTERLIST_SPREADSHEET_ID = '1l2F7wasrVrYSHFIfX0sN2SFkzFz46CCHJfIBXb5R0n4'
PRIMERS_SHEETNAME = 'Primers'
PRIMER_MIXES_SHEETNAME = 'PM and sample name'

def primer_masterlist_sheet2df(sheetname, data_range, index_name):
    data = fetch_data(PRIMER_MASTERLIST_SPREADSHEET_ID, sheetname, data_range)
    df = pd.DataFrame.from_records(data[1:], columns=data[0], index=index_name)
    return df


def fetch_primers_data():
    return primer_masterlist_sheet2df(PRIMERS_SHEETNAME, 'A1:K', 'CO#')


def fetch_primer_mixes_data():
    primer_mixes_df = primer_masterlist_sheet2df(PRIMER_MIXES_SHEETNAME,
                                                 'A1:O',
                                                 'GRNA REC_ID').dropna(how='all')
    return primer_mixes_df


def sgrna2primer_mix_info(sgrna_id_str, info_type, selected_primer_mixes_df):
    #print(sgrna_id_str, info_type)
    items = sgrna_id_str.split('-')
    sgrna_id = items[-1]
    if (sgrna_id_str != sgrna_id_str) or (len(items) < 3): return None # New
    try:
        matching_row = selected_primer_mixes_df.loc[sgrna_id_str]
        if isinstance(matching_row, pd.DataFrame): matching_row = matching_row.iloc[0]

        if info_type == 'primer_mix':
            info = matching_row['Primer Mix #'] #matching_row.name
        elif info_type == 'seq_primer':
            try:
                seq_primer = matching_row['sequecing primer (default)']
                if seq_primer is not None and 'for' in seq_primer:
                    try:
                        parts = seq_primer.split('; ')
                        found = False
                        for part in parts:
                            primer, target_sgrnas_str = part.strip().split(' for ')
                            target_sgrnas = target_sgrnas_str.split(', ')
                            for target_sgrna in target_sgrnas:
                                if target_sgrna == sgrna_id:
                                    seq_primer = primer
                                    found = True
                        if not found:
                            seq_primer = None
                    except Exception as e:
                        seq_primer = None
                if seq_primer is None or len(seq_primer) == 0:
                    warnings.warn(f'{sgrna_id}: No seq_primer found')
                    print(f'{sgrna_id}: No seq_primer found') # new
                    seq_primer = None
            except Exception as e:
                warnings.warn(f'{sgrna_id}: No seq_primer found; error: ' + e)
                print(f'{sgrna_id}: No seq_primer found; error: ' + e)
                seq_primer = None
            """
            to_fro_primers = matching_row.loc[['Forward Primer', 'Reverse Primer']]
            try:
                seq_primer = to_fro_primers[to_fro_primers.str.contains('*', regex=False)][0].replace('*', '')
            except IndexError as e:
                warnings.warn(f'{sgrna_id}: No selected seq_primer found')
                seq_primer = None
            """
            info = seq_primer
        elif info_type == 'tm':
            try:
                #tm = matching_row['Tm (60 = plat II)'].replace('C', '')
                #info = tm
                info = 60 # Dummy Value
            except IndexError as e:
                warnings.warn(f'{sgrna_id}: No tm found')
                print(f'{sgrna_id}: No tm found')
                info = None
        elif info_type == 'size':
            try:
                #size = matching_row.Size
                #info = size
                info = 100 # Dummy Value
            except IndexError as e:
                warnings.warn(f'{sgrna_id}: No size found')
                print(f'{sgrna_id}: No size found')
                info = None
    except IndexError as e:
        warnings.warn(f'{sgrna_id}: No primer_mix found')
        print(f'{sgrna_id}: No primer_mix found')
        info = None
    return info


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


ham_picklist_cols = ['Position', 'Labware_ID', 'Dest_ID', 'Dest_Pos']


def remove_zero_pad_in_wellname_series(ws):
    new_vals = []
    for wellname in ws:
        try:
            row = wellname[0]
            col = int(wellname[1:])
            new_val = f'{row}{col}'
        except:
            new_val = wellname
        new_vals.append(new_val)
    return pd.Series(new_vals, index=ws.index, name=ws.name)


def fetch_reagent(platemap, address):
    row = address[0]
    col = int(address[1:])
    return platemap.loc[row, col]


def add_pos_and_labware_to_ham_primer_mix_df(primer_mix_df, selected_primer_mixes_df):
    plate_nums = []
    primer_wells = []
    for primer_mix in primer_mix_df.primer_or_source_well:
        row = selected_primer_mixes_df[selected_primer_mixes_df['Primer Mix #'] == primer_mix].reset_index(drop=True)
        address = row['primer plate well'][0]
        PM_plate = row['PM plate #'][0]
        primer_wells.append(address)
        plate_nums.append(PM_plate)
    primer_mix_df['Position'] = primer_wells # Get source plate address for the primer
    # primer_mix_df.drop('primer_or_source_well', axis=1, inplace=True)
    primer_mix_df['Source_plate'] = plate_nums # Figure out how to set multiple things here.
    primer_mix_df['Labware_ID'] = ''
    primer_mix_df = primer_mix_df.loc[:, ['Position', 'Labware_ID', 'Dest_Pos', 'primer_or_source_well','Source_plate']]
    return primer_mix_df


def assign_destination(df):
    'Stamping protocol only'
    df['Dest_Pos'] = df['Position']
    return df


def assign_sample_labware(stacked, column):
    '''sample_position = []
    POSITIONS = {}
    st_idx = 0
    for i in range(len(stacked)):
        row = stacked.iloc[i]
        pos = row[column]
        if pos not in POSITIONS:
            st_idx += 1
            POSITIONS[pos] = st_idx
        this_idx = POSITIONS[pos]
        sample_position.append(column+str(this_idx))
    stacked['Labware_ID'] = sample_position'''
    stacked['Labware_ID'] = stacked[column]
    return stacked


def map_to_list(platemap):
    out = []
    rows = platemap.index.tolist()
    columns = platemap.columns.tolist()
    for r in rows:
        for c in columns:
            out.append([r+str(c),platemap.loc[r,c]])
    out_df = pd.DataFrame(out, columns = ['Position','primer_or_source_well']).dropna()
    return out_df


def assign_primer_mixes(df_in, selected_primer_mixes_df):
    df = df_in.copy()
    df.rename(columns = {'primer_or_source_well':'rec_id'}, inplace = True)
    df['primer_or_source_well'] = ''
    for idx, row in df.iterrows():
        rec_id = row['rec_id']
        primer_mix = selected_primer_mixes_df.loc[rec_id,'Primer Mix #']
        if isinstance(primer_mix, pd.Series): primer_mix = primer_mix.iloc[0]
        df.loc[idx,'primer_or_source_well'] = primer_mix
    df = df[['primer_or_source_well','Position']]
    df.rename(columns = {'Position':'Dest_Pos'}, inplace = True)
    return df


def make_hamilton_picklist2_df(consolidated_primer_mix_dfs, consolidated_seq_primer_dfs, platemap, selected_primer_mixes_df, gsheet_primers):

    ham_dna_df = map_to_list(platemap)
    ham_dna_df = ham_dna_df[ham_dna_df['primer_or_source_well'] != 'nan'] # drop Nan that have been typed as string
    ham_primer_mix_df = assign_primer_mixes(ham_dna_df, selected_primer_mixes_df)

    ham_dna_df = assign_destination(ham_dna_df)
    ham_dna_df['Labware_ID'] = 'Plate1' #static assignment for now
    ham_dna_df['Dest_ID'] = 'Destination1' #static assignment for now
    ham_dna_df = sort_dest_col_then_row(ham_dna_df)

    ham_primer_mix_df = add_pos_and_labware_to_ham_primer_mix_df(ham_primer_mix_df, selected_primer_mixes_df)
    ham_primer_mix_df = assign_sample_labware(ham_primer_mix_df, 'Source_plate')
    ham_primer_mix_df['Dest_ID'] = 'Destination1' # static assignment for now
    ham_primer_mix_df['Position'] = remove_zero_pad_in_wellname_series(ham_primer_mix_df['Position'])
    ham_primer_mix_df = sort_dest_col_then_row(ham_primer_mix_df)

    seq_primer_df = assign_seq_primer(ham_primer_mix_df, selected_primer_mixes_df)
    seq_primer_df = get_seq_primer_address(seq_primer_df, gsheet_primers)
    seq_primer_df = assign_sample_labware(seq_primer_df, 'seq_primer_plate')
    seq_primer_df = seq_primer_df[['seq_primer_address','Labware_ID','Dest_ID','Dest_Pos','sequencing_primer','seq_primer_plate']]
    seq_primer_df.rename(columns = {'seq_primer_address':'Position'},inplace = True)
    seq_primer_df['Position'] = remove_zero_pad_in_wellname_series(seq_primer_df['Position'])

    ham_df = pd.concat([ham_dna_df,
                        ham_primer_mix_df], axis=0).loc[:, ham_picklist_cols + ['primer_or_source_well','Source_plate']]
    ham_df.rename(columns={'primer_or_source_well': 'Additional_information'}, inplace=True)
    return ham_df, seq_primer_df


def assign_seq_primer(df, selected_primer_mixes_df):
    out = df.copy()
    seq_primers = []
    for idx in range(len(out)):
        row = out.iloc[idx]#['primer_or_source_well']
        mix_id = row['primer_or_source_well']
        try:
            row = selected_primer_mixes_df[selected_primer_mixes_df['Primer Mix #'] == mix_id].reset_index(drop=True)
            seq_primer = row['sequecing primer (default)'][0]
        except:
            seq_primer = 'error'
        seq_primers.append(seq_primer)
    out['sequencing_primer'] = seq_primers
    return out


def get_seq_primer_address(df, gsheet_primers):
    out = df.copy()
    primer_addresses = []
    primer_plates = []
    for idx in range(len(out)):
        row = out.iloc[idx]
        seq_primer = row['sequencing_primer']
        try:
            primer_line = gsheet_primers.loc[seq_primer]
            if isinstance(primer_line, pd.DataFrame): primer_line = primer_line.iloc[0]
            address = primer_line['Well']
            plate = primer_line['Primer plate/ Seq plate #']
        except:
            address = plate = 'error'
        primer_addresses.append(address)
        primer_plates.append(plate)
    out['seq_primer_address'] = primer_addresses
    out['seq_primer_plate'] = primer_plates
    return out


def sort_dest_col_then_row(df):
    dest_plates = sorted(df['Dest_ID'].unique())
    sorted_dfs = []
    for dest_plate in dest_plates:
        subdf = df.query('Dest_ID == @dest_plate').copy()
        subdf['dest_row'] = [dest_pos[0] for dest_pos in subdf.Dest_Pos]
        subdf['dest_col'] = [''.join(dest_pos[1:]) for dest_pos in subdf.Dest_Pos]
        subdf['dest_col'] = subdf['dest_col'].astype(int)
        sorted_df = subdf.sort_values(by=['dest_col', 'dest_row'])
        sorted_df.drop(['dest_row', 'dest_col'], axis=1, inplace=True)
        sorted_dfs.append(sorted_df)
    return pd.concat(sorted_dfs)


def FormatWell(well):
    try:
        match = re.match(r"([a-z]+)([0-9]+)", well, re.I)
        row,col = match.groups()
        col = '0'+col
        col = col[-2:]
        return row+col
    except:
        return well

# """



def make_picklist(platemap_file, sheet_name, prefix = '', mix_plates = False, seq_plates = False, group_by_temp = False):
    if prefix == '':
        prefix = '.'.join(platemap_file.split('.')[:-1])+'_'+sheet_name
    #"""
    error_file = f'{prefix}.errors.txt'
    error_file_handle = open(error_file, 'w')
    sys.stderr = error_file_handle
    #"""
    # Will work only with REC_ID as input

    source_platemap = pd.read_excel(platemap_file, index_col=0, sheet_name=sheet_name).astype(str)

    selected_primer_mixes_df = fetch_primer_mixes_data()
    if mix_plates:
        first = True
        mix_plates = mix_plates.split(',')
        for i in mix_plates:
            temp = selected_primer_mixes_df[selected_primer_mixes_df['PM plate #'] == i]
            if first:
                final = temp.copy()
                first = False
            else:
                final = pd.concat([final,temp])
        selected_primer_mixes_df = final.copy()

    gsheet_primers = fetch_primers_data()
    if seq_plates:
        first = True
        seq_plates = seq_plates.split(',')
        for i in seq_plates:
            temp = gsheet_primers[gsheet_primers['Primer plate/ Seq plate #'] == i]
            if first:
                final = temp.copy()
                first = False
            else:
                final = pd.concat([final,temp])
        gsheet_primers = final.copy()

    info_types = 'primer_mix seq_primer tm size'.split()
    info_dfs = [source_platemap.applymap(lambda x: sgrna2primer_mix_info(x, info_type, selected_primer_mixes_df)) for info_type in info_types]
    primer_mix_df, seq_primer_df, tm_df, size_df = info_dfs
    if not group_by_temp:
        tm_df.fillna(value=np.nan, inplace=True)
        tm_df[tm_df.notnull()] = '60'
    orig_index = primer_mix_df.index
    orig_columns = primer_mix_df.columns
    tms = tm_df.values.ravel().astype(str)
    tms = tms[tms != 'None']
    tms = tms[tms != 'nan']
    unique_tms = sorted(set(tms))
    ns_primers_per_tm = {}
    for tm in unique_tms:
        ns_primers_per_tm[tm] = (tm_df == tm).sum().sum()
    out_index = list(ascii_uppercase)[:8]
    out_columns = range(1, 13)

    out_df_template = pd.DataFrame(index=out_index, columns=out_columns)

    tms_for_index = []
    out_column_idx = 0
    sorted_tms_and_n_primers = sorted(ns_primers_per_tm.items(), key=lambda x: int(x[0]),
                                      reverse=False)  # [1] would be # guides per temp
    last_tm = sorted_tms_and_n_primers[-1][0]
    all_groups = []
    for tm, n_primers in sorted_tms_and_n_primers:
        is_last_tm = tm == last_tm
        this_tm_mask = tm_df == tm
        rows_cols = np.where(this_tm_mask)
        primer_mixes = primer_mix_df.values[rows_cols]
        seq_primers = seq_primer_df.values[rows_cols]
        sizes = size_df.values[rows_cols]
        source_addresses = []
        for row, col in zip(*rows_cols):
            row_letter = orig_index.tolist()[row]
            col_number = orig_columns.tolist()[col]
            filler = '0' if int(col_number) < 10 else ''
            source_address = '{}{}{}'.format(row_letter, filler, col_number)
            source_addresses.append(source_address)
        zipped = sorted(list(zip(primer_mixes, source_addresses, seq_primers, sizes, [tm] * n_primers)),
                        key=lambda x: int(x[0].replace('PM', '')))
        groups = list(grouper(zipped, 8))
        all_groups.extend(groups)
    groups_of_groups = list(grouper(all_groups, 6))

    consolidated_primer_mix_dfs = [out_df_template.copy() for i in range(len(groups_of_groups))]
    consolidated_seq_primer_dfs = [out_df_template.copy() for i in range(len(groups_of_groups))]
    for pmdf, spdf, gg in zip(consolidated_primer_mix_dfs, consolidated_seq_primer_dfs, groups_of_groups):
        tms_for_index = []
        for col, g in enumerate(gg):
            if g is not None:
                for row, info in enumerate(g):
                    if info is not None:
                        primer_mix, address, seq_primer, size, tm = info
                        if row == 0:
                            tms_for_index.extend([tm] * 2)
                        pmdf.iloc[row, col * 2] = primer_mix
                        pmdf.iloc[row, col * 2 + 1] = address
                        spdf.iloc[row, col * 2] = seq_primer
                        spdf.iloc[row, col * 2 + 1] = size
        n_empty_columns = len(out_columns) - len(tms_for_index)
        tms_for_index.extend(['EMPTY'] * n_empty_columns)
        multi_index_parts = [tms_for_index, out_columns]
        out_multi_index = pd.MultiIndex.from_arrays(multi_index_parts, names=['tm', 'col'])
        pmdf.columns = out_multi_index
        spdf.columns = out_multi_index
    ham_df,seq_df = make_hamilton_picklist2_df(consolidated_primer_mix_dfs, consolidated_seq_primer_dfs,
                                               source_platemap, selected_primer_mixes_df, gsheet_primers)
    ham_df.to_csv(f'{prefix}_hamilton_PCR_picklist.csv', index=False)
    seq_df.to_csv(f'{prefix}_hamilton_seq_picklist.csv', index=False)
    # """
    #error_file_handle.close()
    error_cols = 'primer_mix seq_primer tm size'.split()
    error_df = pd.DataFrame(index=error_cols)
    with open(error_file, 'r') as f:
        for line in f.readlines():
            guide, missing_info = line.strip().split(': ')
            missing_datum = missing_info.split()[1]
            if guide not in error_df.columns:
                error_df[guide] = 0
            error_df.loc[missing_datum, guide] = 1
    error_df = error_df.T

    engine = 'xlsxwriter'  # must have installed
    out_file = f'{prefix}.consolidated_primer_info.xlsx'
    zipped_dfs = zip(consolidated_primer_mix_dfs, consolidated_seq_primer_dfs)
    with pd.ExcelWriter(out_file, engine=engine) as writer:
        for i, (consolidated_primer_mix_df, consolidated_seq_primer_df) in enumerate(zipped_dfs):
            seq_plate_prefix = 'seq_plate{}'.format(i + 1)
            consolidated_primer_mix_dfs[i].to_excel(writer, sheet_name=f'{seq_plate_prefix}.primer_mixes')
            consolidated_seq_primer_dfs[i].to_excel(writer, sheet_name=f'{seq_plate_prefix}.seq_primers')
        error_df.to_excel(writer, sheet_name='errors')
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get primer mixes, sequencing primers, and Tms for a platemap")
    parser.add_argument("platemap_file", help="Excel file containing platemap", type=str)
    parser.add_argument("sheet_name", help="Name of sheet in platemap Excel file", type=str)
    parser.add_argument("--output_prefix", help="Desired prefix for output files", type=str)
    parser.add_argument("--group_by_temp", help="Group guides by their PCR temperature", action='store_true')
    parser.add_argument("--mix_plates", help="Specify a list of plates to use for primer mixes", type=str)
    parser.add_argument("--seq_plates", help="Specify a list of plates to use for sequencing primers", type=str)
    args = parser.parse_args()

    if args.output_prefix:
        prefix = args.output_prefix
    else:
        prefix = '.'.join(args.platemap_file.split('.')[:-1])+'_'+args.sheet_name
    platemap_file = args.platemap_file
    sheet_name = args.sheet_name
    if args.mix_plates: mix_plates = args.mix_plates
    else: mix_plates = False
    if args.seq_plates: seq_plates = args.seq_plates
    else: seq_plates = False
    if args.group_by_temp: group_by_temp = args.group_by_temp
    else: group_by_temp = False
    make_picklist(platemap_file, sheet_name, prefix, mix_plates = mix_plates, seq_plates = seq_plates, group_by_temp = group_by_temp)
