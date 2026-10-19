import pandas as pd

file = r"D:\Major Project\2nd Step\tp53_line1_pipeline\tp53_line1_pipeline\data\GSE33000_raw_data_Extract.txt"

df = pd.read_csv(file, sep="\t", nrows=5)

print("Columns:", df.columns[:8].tolist())
print(df.iloc[:, :5].to_string(index=False))