import pandas as pd

src = r"D:\Bio\inspect_real\score_map_template.csv"
df = pd.read_csv(src)

lower_pathogenic = {
    "EVE", "FATHMM", "LRT", "PoET", "Provean", "SIFT", "SIFT4G", "TranceptEVE_L"
}
df["orientation"] = df["model"].apply(
    lambda m: "lower_pathogenic" if m in lower_pathogenic else "higher_pathogenic"
)
df.to_csv(r"D:\Bio\score_map.csv", index=False)
print(df[["model", "orientation"]].to_string())
