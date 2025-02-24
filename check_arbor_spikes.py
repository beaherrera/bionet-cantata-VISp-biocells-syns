import argparse
import pandas as pd


pd.set_option('display.max_rows', None)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('spikes_csv')

    args = parser.parse_args()
    
    spikes_file = args.spikes_csv
    spikes_df = pd.read_csv(spikes_file, index_col=0)
    print(spikes_df[spikes_df['population'] == 'internal'])
