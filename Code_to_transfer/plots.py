import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def consesus_genome_coverage(consensus_genome, coverage):
    df = pd.DataFrame.from_records(coverage)
    df['0'].astype(str) + df['1']

    #print(df)