# coding : utf-8
import pandas as pd
import glob

# get the columns associated to one HLA and return the subsets columns associated
def find_subsets(netmhc_file):
    # read the netmhcpan output table
    netmhc_table = pd.read_csv(netmhc_file,sep="\t")
    # empty list for all the subsets
    subsets = []
    # empty list for one subset
    sub = []
    # initiate bool to check if met HLA is the first one met
    first = True
    # get column names from the table (mostly unnamed ones)
    column_names = list(netmhc_table.columns)
    # loop through columns
    for col_i in range(len(column_names)):
        # check if a HLA column has not yet been encountered, the columns before the first HLA column are not part of the subsets
        if first :
            # check if the column is a HLA column
            if "HLA" in column_names[col_i]:
                # set first to false
                first = False
                # add the column name to the sub list
                sub.append(column_names[col_i])
        # the last 2 columns are general columns and are not part of the subsets
        elif col_i >= len(column_names)-2:
            # append the subset to the subsets list
            subsets.append(sub)
            # break the loop
            break
        else:
            # test if HLA in column name
            if "HLA" in column_names[col_i]:
                # append subset to subsets list (beginning of new subset means end of previous one)
                subsets.append(sub)
                # reset subset list
                sub = []
            # append column name to subset
            sub.append(column_names[col_i])
    return(subsets,netmhc_table)


def format_netMHCpan(netmhc_table,subsets):
    # empty list that will contain dataframes to concatenate
    to_concat = []
    # get the rows before the first HLA column and the 2 last ones
    common_columns = netmhc_table.iloc[:,:3].copy().join(netmhc_table.iloc[:,-2:].copy())
    # change the column names for the ones in the first row (actual column names in that netmhc output format)
    common_columns.columns = common_columns.iloc[0]
    # drop the first row after the column names switch
    common_columns = common_columns.drop(0)
    print(common_columns)
    # loop through subsets columns
    for sub in subsets:
        # create dataframe containing associated columns
        sub_df = netmhc_table[sub].copy()
        # rename the columns using the first row
        sub_df.columns = sub_df.iloc[0]
        # drop the first row
        sub_df = sub_df.drop(0)
        # create a column containing the HLA (one subset = one HLA so only one HLA for all the rows)
        sub_df["HLA"] = sub[0]
        # join the subset to the common columns
        to_add = common_columns.join(sub_df)
        sub_df = sub_df.drop(["core","icore"],axis=1)
        # add dataframe to the list of dataframes to concatenate
        to_concat.append(to_add)
    # concatenate all dataframes (different HLAs)
    df = pd.concat(to_concat)
    # reset the index and drop the former one
    df = df.reset_index(drop=True)
    return(df)

def filter_netMHC_table(netmhc_table,elr_thr,netmhc_file,NB_min = 0,ELS_thr = 0):
    # convert column types to desired ones
    netmhc_table["NB"] = netmhc_table["NB"].astype(int)
    netmhc_table["EL_Rank"] = netmhc_table["EL_Rank"].astype(float)
    netmhc_table["EL-score"] = netmhc_table["EL-score"].astype(float)
    # filter the NB column (number of WB/SB of peptide accross all HLA)
    # filter the EL_Rank values (after seeing plots, 2 might be a good value)
    # filter the EL-score values (no real filter atm, must just not be 0)
    netmhc_table = netmhc_table[((netmhc_table["NB"]>NB_min)&(netmhc_table["EL_Rank"]<=elr_thr)&(netmhc_table["EL-score"]>ELS_thr))]
    # drop possible duplicates, not taking the netmhcpan index column into account
    netmhc_table = netmhc_table.drop_duplicates(list(netmhc_table.columns)[1:])
    # save to csv
    netmhc_table.to_csv(netmhc_file.split(".out")[0]+f"_netmhc_ELR_{elr_thr}.csv",index=False)
    return netmhc_table


def handle_netMHCpan(netmhc_file,args):
    subsets,netmhc_table = find_subsets(netmhc_file)
    netmhc_df = format_netMHCpan(netmhc_table,subsets)
    netmhc_table = filter_netMHC_table(netmhc_df,args.el_rank,netmhc_file)
    return netmhc_table
