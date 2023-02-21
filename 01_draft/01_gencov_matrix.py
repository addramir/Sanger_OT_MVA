import pandas as pd

# Read the input data from a CSV file
data = pd.read_csv('heart_trait_GenCor.csv')

# Reshape the data from long format to wide format using pivot_table
wide_data = pd.pivot_table(
    data,
    index=['p1'],
    columns=['p2'],
    values=['rg']
)

wide_data.columns = [col[1] for col in wide_data.columns.values]

wide_data[wide_data<0.9]=None

df=wide_data.loc[~wide_data["FINNGEN_R6_FG_CVD"].isna()]

df=df[df.index.values]
df=df.iloc[[0,1,4,5],[0,1,4,5]]
df[df>1]=1

df=df.reset_index()
df.to_csv('gcor.csv', index=False)

h2={"study_id":['FINNGEN_R6_FG_CVD', 'FINNGEN_R6_I9_CVD', 'SAIGE_459', 'SAIGE_459_9'],
    "h2":[0.027,0.0233,0.0091,0.0091],
    "intercept":[1.1275,1.1274,1.0167,1.018]}
h2=pd.DataFrame(h2)
h2.to_csv('h2.csv', index=False)