# James Jensen 2019/2020
# Modified by Chris J, 5/2020

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


def get_sheet_colors(service, wbId: str, ranges: list):
    params = {'spreadsheetId': wbId,
              'ranges': ranges,
              'fields': 'sheets(data(rowData(values(effectiveFormat/textFormat.bold))))'}
    sheet = service.spreadsheets()
    return sheet.get(**params).execute()


def get_primer_mixes_rows_boldness():
    sheet = PRIMER_MIXES_SHEETNAME
    data_range = 'A2:A'  # exclude header row
    sheet_range_str = f'{sheet}!{data_range}'
    desiredA1NotationRanges = [sheet_range_str]
    creds = get_creds()
    service = build('sheets', 'v4', credentials=creds)
    all_data = get_sheet_colors(service, PRIMER_MASTERLIST_SPREADSHEET_ID, desiredA1NotationRanges)
    rowdata = all_data['sheets'][0]['data'][0]['rowData']
    # added if bool(row) to test againts empty dict objects that were being included in rowdata
    rows_are_bold = [row['values'][0]['effectiveFormat']['textFormat']['bold'] for row in rowdata if bool(row)]
    return rows_are_bold


PRIMER_MASTERLIST_SPREADSHEET_ID = '1PP3Nv1u8jhiFCJXvcTdZ2YCoJSmLsJhbimS7Sx2exIo'
PRIMERS_SHEETNAME = 'Primers'
PRIMER_MIXES_SHEETNAME = 'Primer Mixes'


def primer_masterlist_sheet2df(sheetname, data_range, index_name):
    data = fetch_data(PRIMER_MASTERLIST_SPREADSHEET_ID, sheetname, data_range)
    df = pd.DataFrame.from_records(data[1:], columns=data[0], index=index_name)
    return df


def fetch_primers_data():
    return primer_masterlist_sheet2df(PRIMERS_SHEETNAME, 'A1:H', 'Name')


def fetch_primer_mixes_data():
    primer_mixes_df = primer_masterlist_sheet2df(PRIMER_MIXES_SHEETNAME,
                                                 'A1:G',
                                                 'Primer Mix #').dropna(how='all')
    rows_are_bold = get_primer_mixes_rows_boldness()
    primer_mixes_df['bold'] = rows_are_bold[:primer_mixes_df.shape[0]]
    return primer_mixes_df


def sgrna2primer_mix_info(sgrna_id_str, info_type):
    items = sgrna_id_str.split('-')
    if len(items) > 1:
        sgrna_id = items[1]
        if sgrna_id.startswith('C'):
            sgrna_id_n = sgrna_id.replace('C', '')
            regex_str = f'[C|/]-?{sgrna_id_n}(?:\D|$)'
            mix_name_matches_sgrna = selected_primer_mixes_df.Name.str.contains(regex_str, regex=True)
            try:
                matching_row = selected_primer_mixes_df[mix_name_matches_sgrna].iloc[0]
                # print(matching_row.index)
                if info_type == 'primer_mix':
                    info = matching_row.name
                elif info_type == 'seq_primer':
                    try:
                        seq_primer = matching_row['Sequencing Primer']
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
                            seq_primer = None
                    except Exception as e:
                        warnings.warn(f'{sgrna_id}: No seq_primer found; error: ' + e)
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
                        tm = matching_row['Tm (60 = plat II)'].replace('C', '')
                        info = tm
                    except IndexError as e:
                        warnings.warn(f'{sgrna_id}: No tm found')
                        info = None
                elif info_type == 'size':
                    try:
                        size = matching_row.Size
                        info = size
                    except IndexError as e:
                        warnings.warn(f'{sgrna_id}: No size found')
                        info = None
            except IndexError as e:
                warnings.warn(f'{sgrna_id}: No primer_mix found')
                info = None
            return info

    return None


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)


def unique_everseen(items):
    seen = set()
    return [x for x in items if x not in seen and not seen.add(x)]


ham_picklist_cols = ['Position', 'Labware_ID', 'Dest_ID', 'Dest_Pos']


def make_ham_stacked_df(orig_df):
    df = orig_df.copy()
    #col2tm = dict([reversed(tup) for tup in df.columns.tolist()])
    df.columns = df.columns.droplevel(0)
    #tms = unique_everseen([col2tm[col] for col in df.columns])
    df.index.name = 'row'
    stacked = pd.DataFrame(df.stack(), columns=['primer_or_source_well'])
    stacked.reset_index(inplace=True)
    stacked['Dest_Pos'] = stacked['row'].astype(str) + stacked['col'].astype(str)
    return stacked


def remove_zero_pad_in_wellname_series(ws):
    new_vals = []
    for wellname in ws:
        row = wellname[0]
        col = int(wellname[1:])
        new_val = f'{row}{col}'
        new_vals.append(new_val)
    return pd.Series(new_vals, index=ws.index, name=ws.name)


def fetch_reagent(platemap, address):
    row = address[0]
    col = int(address[1:])
    return platemap.loc[row, col]


def make_ham_dna_df(stacked, platemap):
    dna_df = stacked[stacked.col % 2 == 0].copy().sort_values(by=['col', 'row'])
    dna_df.drop(['row', 'col'], axis=1, inplace=True)
    dna_df.rename(columns={'primer_or_source_well': 'Position'}, inplace=True)
    dna_df['Position'] = remove_zero_pad_in_wellname_series(dna_df['Position'])
    dna_df['Labware_ID'] = 'Sample_Source_1'
    dna_df['primer_or_source_well'] = [fetch_reagent(platemap, address) for address in dna_df.Position]
    dna_df = dna_df.loc[:, ['Position', 'Labware_ID', 'Dest_Pos', 'primer_or_source_well']]
    return dna_df


def make_ham_wt_df(stacked):
    wt_df = stacked[stacked.col % 2 == 1].copy().sort_values(by=['col', 'row'])
    wt_df.drop(['row', 'col'], axis=1, inplace=True)
    wt_df['Labware_ID'] = 'WT_DNA'
    pos_cycle = itertools.cycle([f'{letter}1' for letter in ascii_uppercase[:8]])
    wt_df['Position'] = [next(pos_cycle) for _ in range(wt_df.shape[0])]
    wt_df = wt_df.loc[:, ['Position', 'Labware_ID', 'Dest_Pos', 'primer_or_source_well']]
    return wt_df


def make_ham_primer_mix_df(stacked, dna_df, wt_df):
    primer_mix_df = stacked.copy()
    primer_mix_df.loc[dna_df.index, 'primer_or_source_well'] = primer_mix_df.loc[
        wt_df.index, 'primer_or_source_well'].values
    primer_mix_df['col_pair'] = (primer_mix_df.col - 1) // 2
    primer_mix_df.sort_values(by=['col_pair', 'row', 'col'], inplace=True)
    primer_mix_df.drop(['row', 'col', 'col_pair'], axis=1, inplace=True)
    return primer_mix_df


def add_pos_and_sample_tube_to_ham_primer_mix_df(primer_mix_df):
    position_idxs = []
    primers_mixes_seen = set()
    sample_tubes = []
    sample_tube = 1
    position_idx = 1
    for primer_mix in primer_mix_df.primer_or_source_well:
        if primer_mix not in primers_mixes_seen:
            primers_mixes_seen.add(primer_mix)
            position_idx += 1
            if position_idx == 32:
                position_idx = 2
                sample_tube += 1
        position_idxs.append(position_idx)
        sample_tubes.append(f'Sample_Tube_{sample_tube}')

    primer_mix_df['Position'] = position_idxs
    # primer_mix_df.drop('primer_or_source_well', axis=1, inplace=True)
    primer_mix_df['Labware_ID'] = sample_tubes
    primer_mix_df = primer_mix_df.loc[:, ['Position', 'Labware_ID', 'Dest_Pos', 'primer_or_source_well']]
    return primer_mix_df


def assign_destination(stacked, col_type):
    dest = []
    dest_idx = 0
    for idx, row in stacked.iterrows():
        dest_pos = row['Dest_Pos']
        if dest_pos == ('A2' if col_type == 'even' else 'A1'):
            dest_idx += 1
        dest.append(f'Destination{dest_idx}')
    stacked['Dest_ID'] = dest
    return stacked


def assign_sample_tubes(stacked):
    sample_tubes = []
    st_idx = 0
    for i in range(len(stacked)):
        row = stacked.iloc[i]
        pos = row['Position']
        if pos == 2 and (i == 0 or stacked.iloc[i - 1]['Position'] == 31):
            st_idx += 1
        sample_tubes.append(f'Sample_Tube_{st_idx}')
    stacked['Labware_ID'] = sample_tubes
    return stacked


def make_hamilton_picklist_df(consolidated_primer_mix_dfs, platemap):
    ham_dna_dfs = []
    ham_wt_dfs = []
    ham_primer_mix_dfs = []
    for i, orig_df in enumerate(consolidated_primer_mix_dfs):
        stacked = make_ham_stacked_df(orig_df)
        ham_dna_df = make_ham_dna_df(stacked, platemap)
        ham_dna_dfs.append(ham_dna_df)
        ham_wt_df = make_ham_wt_df(stacked)
        ham_wt_dfs.append(ham_wt_df)
        ham_primer_mix_df = make_ham_primer_mix_df(stacked, ham_dna_df, ham_wt_df)
        ham_primer_mix_dfs.append(ham_primer_mix_df)
    ham_dna_df = assign_destination(pd.concat(ham_dna_dfs), 'even')
    ham_wt_df = assign_destination(pd.concat(ham_wt_dfs), 'odd')
    ham_primer_mix_df = pd.concat(ham_primer_mix_dfs)
    ham_primer_mix_df = add_pos_and_sample_tube_to_ham_primer_mix_df(ham_primer_mix_df)
    ham_primer_mix_df = assign_destination(ham_primer_mix_df, 'odd')
    ham_primer_mix_df = assign_sample_tubes(ham_primer_mix_df)
    ham_primer_mix_df = sort_dest_col_then_row(ham_primer_mix_df)

    ham_df = pd.concat([ham_dna_df,
                        ham_wt_df,
                        ham_primer_mix_df], axis=0).loc[:, ham_picklist_cols + ['primer_or_source_well']]
    ham_df.rename(columns={'primer_or_source_well': 'primer_mix'}, inplace=True)
    return ham_df


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
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get primer mixes, sequencing primers, and Tms for a platemap")
    parser.add_argument("platemap_file", help="Excel file containing platemap", type=str)
    parser.add_argument("sheet_name", help="Name of sheet in platemap Excel file", type=str)
    parser.add_argument("--output_prefix", help="Desired prefix for output files", type=str)
    parser.add_argument("--group_by_temp", help="Group guides by their PCR temperature", action='store_true')
    args = parser.parse_args()
    """;

    class FakeArgs:

        def __init__(self, platemap_file, sheet_name, output_prefix):
            self.output_prefix = output_prefix
            self.platemap_file = platemap_file
            self.sheet_name = sheet_name

    args = FakeArgs('phenosprint1.q4.2019.source_and_dest_plate_info.xlsx',
                    'genoscreen1', 'q4.2019.phenosprint1.genoscreen1')
    #"""
    if args.output_prefix:
        prefix = args.output_prefix
    else:
        #prefix = os.path.splitext(args.platemap_file)[0]+'.'+args.sheet_name
        prefix = args.platemap_file.split('.')[0]+'.'+args.sheet_name
    platemap_file = args.platemap_file

    #"""
    error_file = f'{prefix}.errors.txt'
    error_file_handle = open(error_file, 'w')
    sys.stderr = error_file_handle
    #"""
    source_platemap = pd.read_excel(platemap_file, index_col=0, sheet_name=args.sheet_name).astype(str)
    source_platemap = source_platemap.applymap(
        lambda x: x if len(x.split('-')) < 3 else '-'.join([x.split('-')[0], x.split('-')[2]]))
    gsheets_primer_mixes_df = fetch_primer_mixes_data()
    selected_primer_mixes_df = gsheets_primer_mixes_df[gsheets_primer_mixes_df.bold]
    info_types = 'primer_mix seq_primer tm size'.split()
    info_dfs = [source_platemap.applymap(lambda x: sgrna2primer_mix_info(x, info_type)) for info_type in info_types]
    primer_mix_df, seq_primer_df, tm_df, size_df = info_dfs
    if not args.group_by_temp:
        tm_df[tm_df.notnull()] = '60'
    orig_index = primer_mix_df.index
    orig_columns = primer_mix_df.columns
    tms = tm_df.values.ravel().astype(str)
    tms = tms[tms != 'None']
    unique_tms = sorted(set(tms))
    ns_primers_per_tm = {}
    for tm in unique_tms:
        ns_primers_per_tm[tm] = (tm_df == tm).sum().sum()
    out_index = list(ascii_uppercase)[:8]
    out_columns = range(1, 13)

    out_df_template = pd.DataFrame(index=out_index, columns=out_columns)

    consolidated_primer_mix_dfs = []
    consolidated_seq_primer_dfs = []
    consolidated_primer_mix_df = out_df_template.copy()
    consolidated_seq_primer_df = out_df_template.copy()
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
        zipped = sorted(list(zip(primer_mixes, source_addresses, seq_primers, sizes, [tm] * n_primers)), key=lambda x: int(x[0].replace('PM', '')))
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

    """
    test_primer_mix_df = source_platemap = pd.read_excel('Hamilton picklist setup from DS script.xlsx', skiprows=12, index_col=0).iloc[1:9, :12]
    tpmdf_t = test_primer_mix_df.T
    tpmdf_t['tm'] = [65] * 8 + [62] * 4
    tpmdf_t.reset_index(inplace=True)
    tpmdf_t.set_index('tm', inplace=True)
    tpmdf_t.rename(columns={'index':'col'}, inplace=True)
    tpmdf_t.set_index('col', inplace=True, append=True)
    tpmdf = tpmdf_t.T
    consolidated_primer_mix_dfs = [tpmdf]
    """
    ham_df = make_hamilton_picklist_df(consolidated_primer_mix_dfs, source_platemap)
    '''
    # reformat well identifiers to have 2 digit column descriptions: A03 instead of A3
    for i in range(len(ham_df)):
        ham_df.iloc[i,0] = FormatWell(ham_df.iloc[i,0])
        ham_df.iloc[i,3] = FormatWell(ham_df.iloc[i,3])
    '''
    """
    with pd.option_context('display.max_rows', None):
        display(ham_df)
    """
    ham_df.to_excel(f'{prefix}.hamilton_picklist.xls', index=False)
    # """
    error_file_handle.close()
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
    # """
    engine = 'xlsxwriter'  # must have installed
    out_file = f'{prefix}.consolidated_primer_info.xlsx'
    zipped_dfs = zip(consolidated_primer_mix_dfs, consolidated_seq_primer_dfs)
    with pd.ExcelWriter(out_file, engine=engine) as writer:
        for i, (consolidated_primer_mix_df, consolidated_seq_primer_df) in enumerate(zipped_dfs):
            seq_plate_prefix = 'seq_plate{}'.format(i + 1)
            consolidated_primer_mix_dfs[i].to_excel(writer, sheet_name=f'{seq_plate_prefix}.primer_mixes')
            consolidated_seq_primer_dfs[i].to_excel(writer, sheet_name=f'{seq_plate_prefix}.seq_primers')
        error_df.to_excel(writer, sheet_name='errors')
