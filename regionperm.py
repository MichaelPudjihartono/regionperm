import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm
import argparse
import os
import sys
import time
import datetime


def parse_args():
    parser = argparse.ArgumentParser(description="""regionperm: Genomic region association analysis with permutation tests.
        This program allows you to assess the association between a set of genomic regions and other genomic features using permutation tests.""")
    parser.add_argument('-A', required=True, help='''File containing the set of regions to randomize in BED format. 
        File A is subset of a finite set of all valid regions in the universe file.
        For example, a small set of genes (file A) as a subset of all genes in the genome (universe file).''')
    parser.add_argument('-B', required=True, help='File containing the set of genomic features that will be analyzed for their association with the regions in file A in BED format.')
    parser.add_argument('-U', '--universe', required=True, help='File containing the total valid regions from which subsequent random iterations will be sampled from in BED format. File A is a subset of the universe file')
    parser.add_argument('-n', '--num-iterations', type=int, default=1000, help='Number of iterations for permutation test.')
    parser.add_argument('-m', '--match-by', choices=['count', 'length'], default='count', help="For each iteration, do you want to randomize regions by matching the count or length of the original A region set?. Default: 'count'")
    parser.add_argument('-o', '--output-dir', required=True, help='Directory to write results.')

    return parser.parse_args()

class Logger:
    def __init__(self, output_dir):
        self.console = sys.stdout
        self.log = open(os.path.join(output_dir, 'regionperm.log'), 'w')

        time = datetime.datetime.now()
        self.write((f'{time.strftime("%A")},'
                    f' {time.strftime("%b")}' 
                    f' {time.strftime("%d")},'
                    f' {time.strftime("%Y")}' 
                    f' {time.strftime("%I")}:'
                    f'{time.strftime("%M")}'
                    f' {time.strftime("%p")}'))
        
    def write(self, message):
        self.console.write(message+'\n')
        self.log.write(message+'\n')
        self.log.flush()


def read_input(file_path):
    df = pd.read_csv(file_path, sep="\t")
    df = df.iloc[:, :3]
    return df

def add_fragment_length_column(df):
    df['fragment_length'] = df.apply(lambda row: row.iloc[2]-row.iloc[1], axis=1)
    return df

def add_feature_id_column(df):
    df['feature_id'] = df.apply(lambda row: str(row.iloc[0])+':'+str(row.iloc[1])+'-'+str(row.iloc[2]), axis=1)
    return df

def get_universe_intersect(universe_region_set, B_region_set):
    universe_region_set_BEDtool = pybedtools.BedTool.from_dataframe(universe_region_set)
    B_region_set_BEDtool = pybedtools.BedTool.from_dataframe(B_region_set)

    universe_intersect = universe_region_set_BEDtool.intersect(B_region_set_BEDtool, wa=True, wb=True)
    universe_intersect = universe_intersect.to_dataframe()
    universe_test_stat = len(set(universe_intersect.iloc[:,-1])) #The last column correspond to 'feature_id'
    #Thus, test statistic counts the total number of unique B region set that intersects the universe region set

    return universe_intersect, universe_test_stat

def get_intersect(region_set, universe_intersect):
    left_on = list(region_set.columns[0:3])
    right_on = list(universe_intersect.columns[0:3])

    intersect = pd.merge(region_set, universe_intersect, how='inner', left_on=left_on, right_on=right_on)
    test_stat = len(set(intersect.iloc[:,-1])) #The last column correspond to 'feature_id'
    #Thus, test statistic counts the total number of unique B region set that intersects the specified region set

    return test_stat

def create_complement_region_set(A_region_set, universe_region_set): #This creates the complement of A_region_set
    left_on = list(A_region_set.columns[0:3])
    right_on = list(universe_region_set.columns[0:3])

    complement_region_set = pd.merge(A_region_set, universe_region_set, how='outer', left_on=left_on, right_on=right_on)
    complement_region_set = complement_region_set[complement_region_set.isnull().any(axis=1)]

    return complement_region_set

def permutation_simulation(A_region_set, universe_region_set, universe_intersect, 
                           original_test_stat, num_iterations, match_by):
    
    perm_test_stats = [original_test_stat]
    more_than_original_count = 1
    less_than_original_count = 1

    bar_format = '{desc}: {percentage:3.0f}% |{bar}| {n_fmt}/{total_fmt} {unit}'

    for iteration in tqdm(
        range(num_iterations), 
        desc='Permutation simulations in progress', 
        unit='simulation', 
        bar_format = bar_format, 
        ncols=100):

        iteration_region_set = universe_region_set.sample(frac=1).reset_index().drop('index', axis=1)
        iteration_region_set['fragment_length_cumsum'] = iteration_region_set['fragment_length'].cumsum()

        if match_by == 'count':
            total_count = A_region_set.shape[0]
            up_to_index = total_count
        elif match_by == 'length':
            total_length = A_region_set['fragment_length'].sum()
            up_to_index = iteration_region_set[iteration_region_set['fragment_length_cumsum'] >= total_length].index[0] + 1
        
        iteration_region_set = iteration_region_set.iloc[:up_to_index]
        iteration_test_stat = get_intersect(iteration_region_set, universe_intersect)

        perm_test_stats.append(iteration_test_stat)

        if iteration_test_stat >= original_test_stat:
            more_than_original_count += 1
        elif iteration_test_stat <= original_test_stat:
            less_than_original_count += 1

    return perm_test_stats, more_than_original_count, less_than_original_count

def generate_perm_test_result(perm_test_stats, original_test_stat, complement_test_stat, 
                              universe_test_stat, more_than_original_count, less_than_original_count, args):
    
    alternative = 'greater' if more_than_original_count<less_than_original_count else 'less'

    if alternative == 'greater':
        p_val = more_than_original_count/len(perm_test_stats)
    elif alternative == 'less':
        p_val = less_than_original_count/len(perm_test_stats)
    
    #Calculate z-score
    mean = sum(perm_test_stats) / len(perm_test_stats)
    std_dev = (sum((x - mean) ** 2 for x in perm_test_stats) / len(perm_test_stats)) ** 0.5
    z_score = (original_test_stat - mean) / std_dev

    #Generate perm_test_result df
    perm_test_result = pd.DataFrame({'A_region_set':[f'{os.path.basename(args.A)}'], 
                                     'B_region_set':[f'{os.path.basename(args.B)}'],
                                     'match':[f'{args.match_by}'],
                                     'p_value':[p_val],
                                     'z_score':[z_score],
                                     'n_iterations':[f'{args.num_iterations}'],
                                     'alternative':[alternative],
                                     'original_evaluation':[original_test_stat],
                                     'complement_evaluation':[complement_test_stat],
                                     'universe_evaluation':[universe_test_stat]})
    
    perm_test_result.to_csv(os.path.join(args.output_dir, 'perm_test_result.txt'), sep="\t", index=False)

    return perm_test_result

def generate_perm_test_raw_file(perm_test_stats, num_iterations, output_dir):
    identity = [f'iteration_{i+1}' for i in range(num_iterations)]
    identity = ['original_evaluation'] + identity
    perm_test_raw_data = pd.DataFrame({'Test_stat_evaluation':perm_test_stats,
                                      'Identity':identity})
    
    perm_test_raw_data.to_csv(os.path.join(output_dir, 'perm_test_raw_data.txt'), sep="\t", index=False)

def generate_perm_test_figure(perm_test_stats, perm_test_result, output_dir):
    fig, ax = plt.subplots(figsize=(8, 6))

    #Plot histogram
    sns.histplot(data=perm_test_stats, bins=100, kde=False, color='#a5aab0', stat='density', alpha=0.5, ax=ax)

    #Set axis labels and texts
    ax.set_xlabel('Number of overlaps', weight='bold')
    ax.set_ylabel('Density', weight='bold')

    text = f"p-value: {perm_test_result.loc[0,'p_value']}\nz-score: {perm_test_result.loc[0,'z_score']}\nn-iterations: {perm_test_result.loc[0, 'n_iterations']}"
    ax.text(x=0.5, y=1.05, s=text, transform=ax.transAxes, ha='center')

    #Add vertical line to indicate significance_cutoff and observed value
    percentile = 95 if perm_test_result.loc[0, 'alternative'] == 'greater' else 5
    significance_cutoff = np.percentile(perm_test_stats, percentile)
    ax.axvline(x=significance_cutoff, ymax=0.7, color= 'red', linestyle='--')
    ax.text(x=significance_cutoff, y=ax.get_ylim()[1]*0.72, s=f'{int(significance_cutoff)}',
           color='red', va='center', ha='center', fontsize=8, weight='bold')
    
    ax.axvline(x=perm_test_result.loc[0, 'original_evaluation'], ymax=0.7, color='green', linestyle='--')
    ax.text(x=perm_test_result.loc[0, 'original_evaluation'], y=ax.get_ylim()[1]*0.72, s=f"{int(perm_test_result.loc[0, 'original_evaluation'])}",
           color='green', va='center', ha='center', fontsize=8, weight='bold')
    
    ax.legend(['\u03B1 = 0.05', 'Observed'])

    #Plot KDE line
    sns.kdeplot(data=perm_test_stats, color='black', ax=ax)
    ax.axvline(x=significance_cutoff, ymax=0.7, color= 'red', linestyle='--')
    ax.axvline(x=perm_test_result.loc[0, 'original_evaluation'], ymax=0.7, color='green', linestyle='--')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'perm_test_figure.png'), dpi=300, facecolor='white')


if __name__=='__main__':
    start_time = time.time()
    pd.options.mode.chained_assignment = None

    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    logger = Logger(args.output_dir)

    logger.write('SETTINGS\n========')
    for arg in vars(args):
        logger.write(f'{arg}:\t {vars(args)[arg]}')
    logger.write('\n')

    logger.write('Processing inputs...')
    A_region_set = read_input(args.A)
    A_region_set = add_fragment_length_column(A_region_set)

    B_region_set = read_input(args.B)
    B_region_set = add_feature_id_column(B_region_set)

    universe_region_set = read_input(args.universe)
    universe_region_set = add_fragment_length_column(universe_region_set)

    universe_intersect, universe_test_stat = get_universe_intersect(universe_region_set, B_region_set)

    original_test_stat = get_intersect(A_region_set, universe_intersect)

    complement_region_set = create_complement_region_set(A_region_set, universe_region_set)
    complement_test_stat = get_intersect(complement_region_set, universe_intersect)


    logger.write('Starting simulations...')
    perm_test_stats, more_than_original_count, less_than_original_count = permutation_simulation(A_region_set, universe_region_set, universe_intersect, 
                                                                                                 original_test_stat, args.num_iterations, args.match_by)
    
    logger.write('Writing results...\n')
    perm_test_result = generate_perm_test_result(perm_test_stats, original_test_stat, complement_test_stat, 
                                                 universe_test_stat, more_than_original_count, less_than_original_count, args)
    
    
    generate_perm_test_raw_file(perm_test_stats, args.num_iterations, args.output_dir)
    generate_perm_test_figure(perm_test_stats, perm_test_result, args.output_dir)

    logger.write(f'Original evaluation: {original_test_stat}')
    logger.write(f"Alternative: {perm_test_result.loc[0,'alternative']}")
    logger.write(f"P-value: {perm_test_result.loc[0,'p_value']}")
    logger.write(f"Z-score: {perm_test_result.loc[0,'z_score']}\n")

    logger.write('Done.')
    logger.write(f'Time elapsed: {(time.time() - start_time) / 60:.2f} minutes.')
