import streamlit as st
import pandas as pd
import plotly.express as px

st.title(" 5MW Electrolyzer Criticality Analysis")

costWeight = st.slider("Cost Weight", 0.0, 1.0, 0.5, 0.05)
reliabilityWeight = 1 - costWeight
st.write(f"Reliability Weight: **{reliabilityWeight:.2f}**")

data = [
    ["MV/LV Transformer", 450000, 120000, 9.5, 9.0],
    ["Converter Stack (IGBTs)", 380000, 180000, 9.0, 8.5],
    ["MV Switchgear", 280000, 85000, 7.5, 8.0],
    ["Cooling System (PE + Transformer)", 265000, 225000, 7.0, 8],  
    ["DC Bus Capacitors", 95000, 140000, 6.0, 7.0],
    ["Gate Drivers & Protection", 85000, 60000, 5.5, 8.5],
    ["Harmonic Filters", 180000, 95000, 6.5, 6.0]
]

df = pd.DataFrame(data, columns=["Component", "CAPEX", "OPEX15yr", "CostScore", "ReliabilityScore"])

df["TotalCost"] = df["CAPEX"] + df["OPEX15yr"]
df["OPEX_Ratio"] = df["OPEX15yr"] / df["TotalCost"]

df["CostImpact"] = (df["CostScore"] * 0.6) + (df["OPEX_Ratio"] * 10 * 0.4)
df["CI"] = (costWeight * df["CostImpact"]) + (reliabilityWeight * df["ReliabilityScore"])
df = df.sort_values(by="CI", ascending=False)
df["Rank"] = range(1, len(df) + 1)

st.subheader("📊 Criticality Ranking")
st.dataframe(df[['Rank','Component','CI','CAPEX','OPEX15yr','CostImpact','ReliabilityScore']])

st.subheader("🔥 Criticality Index Bar Chart")
fig_bar = px.bar(df, x="CI", y="Component", orientation="h", color="CI", color_continuous_scale="RdYlGn_r")
fig_bar.update_layout(xaxis_title="Criticality Index", yaxis_title="")
st.plotly_chart(fig_bar, width="stretch")

st.subheader("⚖️ Cost vs Reliability Matrix")
fig_scatter = px.scatter(df, x="CostImpact", y="ReliabilityScore", size="CI", color="Rank",
                         text="Component", color_continuous_scale="RdYlGn_r")
fig_scatter.update_traces(textposition="top center")
fig_scatter.update_layout(xaxis_title="Cost Impact Score", yaxis_title="Reliability Score")
st.plotly_chart(fig_scatter, width="stretch")
