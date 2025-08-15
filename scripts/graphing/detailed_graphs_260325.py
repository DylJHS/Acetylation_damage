# Load the necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import datetime as dt


# Read the csv file
df = pd.read_csv('/home/djhs/wsl_projects/Dros_H3K9ac_bulkChIC_Analysis/data/tables/Detailed_Alignment_and_BLAST_Stats_with_mapping.csv')

# Today's date
today = dt.datetime.today().strftime('%Y-%m-%d')
print(
    '\n\t\t\t Today\'s date: \n',
    today,
    '\n'
)

# Print the details about the dataframe
print(
    '\n\t\t\t Info on the df: \n',
      '\n',
      df.describe(),
      '\n',
      df.head(),
        '\n'
)

# Set the sample as the index
df.set_index('Sample', inplace=True)

# Remove the 004 sample
df.drop('SCC-bulkChIC-UMC-JAN-004', inplace=True)

# Remove the columns that contain "Gap Openings" or "Mismatches"
df = df[[col for col in df.columns if 'Gap Openings' not in col and 'Mismatches' not in col]]

# Break the dataframe into two dataframes, the fly one and the human one based on columns having "Human" or "Fly" in the name
df_fly = df[[col for col in df.columns if 'Fly' in col]].reset_index()
df_human = df[[col for col in df.columns if 'Human' in col]].reset_index()

# Remove the word "Fly" or "Human" from the column names
df_fly.columns = [col.replace('Fly', '').strip("_").strip() for col in df_fly.columns]
df_human.columns = [col.replace('Human', '').strip("_").strip() for col in df_human.columns]
print(
    '\n\t\t\t Columns in the two dataframes: \n',
    df_fly.columns,
      '\n', 
      df_human.columns,
        '\n'
)

# Add a new column to the dataframes that contains the organism name
df_fly['Reference Organism'] = 'Drosophila'
df_human['Reference Organism'] = 'Human'

# Check that the columns in the two dataframes are the same
print(
    '\n\t\t\t Check for column equality',
      df_fly.columns == df_human.columns,
        '\n'
)

# Concatenate the two dataframes
df_combined = pd.concat([df_fly, df_human])
print(
    '\n\t\t\t Combined dataframe: \n',
      df_combined.head(),
        '\n'
)

# Export the combined dataframe to a csv file
# df_combined.to_csv('/home/djhs/wsl_projects/Dros_H3K9ac_bulkChIC_Analysis/data/combined_alignment_stats.csv', index=False)

# Modify the column names
df_combined['Sample'] = df_combined['Sample'].apply(lambda x: ('-').join(x.split('-')[-2:]))
df_combined['% Reads Mapped'] = df_combined['% Mapped']

# Log transform the Read Count columns and total hits columns
df_combined['R1 Total Hits (log2+1)'] = np.log2(df_combined['R1 Total Hits'] + 1)
df_combined['R2 Total Hits (log2+1)'] = np.log2(df_combined['R2 Total Hits'] + 1)
# Drop the original columns
df_combined.drop([
    'R1 Total Hits', 
    'R2 Total Hits', 
    'Properly Paired (%)', 
    'R1 Avg Bit Score', 
    'R2 Avg Bit Score',
    '% Mapped'
    ], axis=1, inplace=True
  )
# Desired order
df_combined = df_combined[[
    'Sample', 
    'Reference Organism',
    'QC-passed Reads',
    '% Reads Mapped',
    'R1 Total Hits (log2+1)',
    'R2 Total Hits (log2+1)',
    'R1 Alignment Rate (%)',
    'R2 Alignment Rate (%)',
    'R1 Avg % Identity',
    'R2 Avg % Identity'
    ]]

print(
    '\n\t\t\t Log transformed columns: \n',
    df_combined.head(),
    '\n'
)


# Melting the dataframe
df_melt = pd.melt(df_combined, id_vars=['Sample', 'Reference Organism'], var_name='Alignment Stat', value_name='Value')
print(
    '\n\t\t\t Melted dataframe: \n',
    df_melt.head(),
    '\n'
)

g = sns.catplot(
    data=df_melt,
    x='Sample',
    y='Value',
    hue='Reference Organism',
    col='Alignment Stat',              
    kind='bar',
    col_wrap=2,               
    height=4,
    aspect=1.2,
    palette='Set2',
    sharex=False,
    sharey=False,
    legend = True
)

# Improve formatting
g.set_titles(col_template="{col_name}", size=12)
g.set_axis_labels("Sample", "Value", fontsize=12)
g.set_xticklabels(fontsize=10)
g.fig.subplots_adjust(top=0.9, wspace=0.1, hspace=0.3)
g.fig.suptitle("Alignment Stats by Sample", fontsize=18, fontweight='bold')
plt.show()

# Save the full figure
g.savefig(f'/home/djhs/wsl_projects/Dros_H3K9ac_bulkChIC_Analysis/figures/alignment_full_stats_faceted_{today}.png')
