import pandas as pd
import math

# --- Global Paramters --- 
# Helps keep tpfs together enough to minimise error
netsize = 3

# How many texels smaller than the smallest contig we will search for 
lowcutoff = 1.5

# Sex Chromosome IDS
sex = ["X","Y","Z","W"]
prefix = 'R'
borderlen = 80

agp_file = '/nfs/treeoflife-01/teams/tola/users/dp24/rapid-pretext/input-data/idSyrVitr1_1.pretext.agp'
tpf_file = '/nfs/treeoflife-01/teams/tola/users/dp24/rapid-pretext/input-data/idSyrVitr1.20221125.decontaminated.fa.tpf'


def read_agp (agp):
    """
    Reads the input agp, skips the 3*# lines, filters for lines that arn't gap scaffs, then drops the end columns.
    """
    agp_df = pd.read_csv(agp,sep = '\t',skiprows=3, names=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'])
    agp_df = agp_df[agp_df['e'] == 'W']
    #agp_df.drop(inplace=True, columns=['k', 'l'])

    return agp_df


def read_tpf (tpf):
    """
    Reads the input tpf file
    """
    tpf_df = pd.read_csv(tpf,sep = '\t', names=['a', 'b', 'c', 'd'])
    tpf_df = tpf_df[tpf_df['a'] == '?']
    tpf_df.drop(inplace=True, columns=['a'])

    return tpf_df


def reformat_tpf_df (df):
    """
    Reformat the tpf dataframe so this: `scaffold_1:1-2795873` becomes `scaffold_1  1  2795873 {length of scaff}`
    e.g. 1 column becomes 4
    """
    df[['scaffold', 'coords']] = df.b.str.split(":", expand = True)
    df[['tpf_start', 'tpf_end']] = df.coords.str.split("-", expand = True)
    df[['tpf_start','tpf_end']].convert_dtypes().astype(int)

    # Temp df required for scaff length calc
    temp_df = df[['b']]
    temp_df['length'] = (
        df['tpf_end'].astype('int32') - df['tpf_start'].astype('int32')
    )

    df = pd.concat([df, temp_df], axis=1, join='outer')
    df.drop(inplace=True, columns=['b', 'c', 'coords'])
    temp_df = [] # Deletes temp_df to save mem

    return df


def append_dict(k,v,d):
    if k not in d:
        d[k]=[v]
    else:
        d[k].append(v)
    return d


def dividers (agp_df):
    dividers = {}
    return [append_dict(v.f, v.h, dividers) for i, v in agp_df.iterrows()]


def texel_calc (tpf):
    """
    Count the size of genome based on scaffold lengths
    Calculates the number of texels in pretextmap
    """
    size = tpf['length'].sum()
    
    return size, int(round(size/32768,0))


def filter_keepers (agp_df, tpf_df, texel_size):
    """
    Discard all small fragments
    keep all scaffs > 10 texels in length calculate whether to keep
    """
    delete_list = []
    haplotig_list = []
    sex_list = []
    scaff_dict = {}
    for i, v in agp_df.iterrows():
        tpf = tpf_df[tpf_df['scaffold'] == v.f]
        # if superscaff end - start > 10 * texel_size
        if (v.c - v.b) > (10 * texel_size):
            append_dict(v.a,[v.f,v.g,v.h,v.i,v.b,v.c],scaff_dict)
        # TODO: should this just be against the largest (first) scaff for a scaff?
        elif v.c - v.b > [vv.length for ii, vv in tpf.iterrows()][0] - lowcutoff*texel_size:
            # Showcases a need to rename columns
            append_dict(v.a,[v.f,v.g,v.h,v.i,v.b,v.c],scaff_dict)
        else:
            # everything else isn't needed
            pass
        
        # No haplotig data in dataset
        if not pd.isnull(v.j):
            if 'HAPLOTIG' in v.j:
                haplotig.append(v.j)
            for xxx in sex_list:
                if xxx in v.j:
                    print(f"It's a {xxx}!!")
                    sex_list.append(xxx)

    return delete_list, haplotig_list, sex_list, scaff_dict


def agp_dividers(agpdict):

    tdivs={}
    result={}
    for k,v in agpdict.items():
        [append_dict(i[0],int(i[2]),tdivs) for i in v]

    for k, v in tdivs.items():
        result[k]=sorted(v)
    return result


def double_check (tpf_end, agp_end, scaff_id, dic):
    if abs(tpf_end - agp_end) < abs(tpf_end-dic[scaff_id]):
        dic[scaff_id] = tpf_end
    return dic

def get_nearest (tpf, dividers, texel_size):
    closest = {}
    for k, v in dividers.items():
        for scaff_end in v:
            tpf['tpf_end'] = tpf['tpf_end'].astype(int)
            for index, line in tpf.iterrows():
                ids = f'{k}:{str(scaff_end)}'
                if abs(line.tpf_end - scaff_end) < netsize*texel_size:
                    if ids not in closest:
                        closest[ids]=line.tpf_end
                    else:
                        double_check(line.tpf_end, scaff_end, ids, closest)
    return closest

def main():
    agp_df = read_agp(agp_file)
    tpf_df = read_tpf(tpf_file)
    reformatted_tpf = reformat_tpf_df(tpf_df)
    genome_size, texel_size = texel_calc(reformatted_tpf)
    delete_list, hap_list, sex_list, scaff_dict = filter_keepers(agp_df, reformatted_tpf, texel_size)
    dividers = agp_dividers(scaff_dict)
    agp2tpf = get_nearest(reformatted_tpf, dividers, texel_size)


