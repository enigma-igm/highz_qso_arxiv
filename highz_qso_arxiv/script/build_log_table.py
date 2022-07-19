import numpy as np
import pandas as pd

df = pd.read_csv('../arxiv/LRIS_2203/log/LRIS_2203.csv')
df = df[["Unnamed: 0", "My Notes", "Notes", "ra", "dec", "mag"]]
latex_table = df.to_latex(index=False)
print(latex_table)