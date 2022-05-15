# REACT dates
# BA.1 first on 3rd December
# BA.1.1 first on 9th December
# BA.2 first on 7th January
# Delta first on 2021-04-27
# REACT -- defined as 30 days prior to first appearance in REACT


# UKSHA dates
# Delta 18 Mar 2021
# Omicron 6 Nov 2021
# BA.2 22 Dec 2021
# UKSHA -- dates defined as first day with cumulative 5 cases of variant


## to look into reinfection by proxy
## definition: get dates and covida < 1 month


import pandas as pd
import numpy as np
from scipy import stats

df18 = pd.read_csv("E:/dt20/data/rep17_18_lineage.csv",
                        parse_dates = ["agprev5_new", "d_comb"])

df13 = pd.read_csv("E:/dt20/data/rep13_lineage.csv",
                        parse_dates = ["agprev5_new", "d_comb"],
                   encoding='unicode_escape')



def processDfOmicron(df, timedelta = 1, reactDates=False):

    #df["timedelta_swabs"] = (df["d_comb"] - df["agprev5_new"]).dt.days

    if reactDates == False:
        reinfectionCriteria = \
            [
                (df["react_lineage"] == "BA.1") & \
                ((pd.to_datetime("2021-11-06") - df["agprev5_new"]).dt.days >= timedelta),
                (df["react_lineage"] == "BA.1.1") & \
                ((pd.to_datetime("2021-11-06") - df["agprev5_new"]).dt.days>= timedelta),
                (df["react_lineage"] == "BA.2") & \
                ((pd.to_datetime("2021-12-22") - df["agprev5_new"]).dt.days>= timedelta)
            ]
    else:
        reinfectionCriteria = \
            [
                (df["react_lineage"] == "BA.1") & \
                ((pd.to_datetime("2021-12-03") - df["agprev5_new"]).dt.days>= timedelta),
                (df["react_lineage"] == "BA.1.1") & \
                ((pd.to_datetime("2021-12-09") - df["agprev5_new"]).dt.days>= timedelta),
                (df["react_lineage"] == "BA.2") & \
                ((pd.to_datetime("2022-01-07") - df["agprev5_new"]).dt.days>= timedelta)
            ]

    reinfectionChoices = [True, True, True]
            
    reportedDate = pd.notnull(df["agprev5_new"])

    df["reinfection"] = np.select(reinfectionCriteria,
                                         reinfectionChoices,
                                         default=False)

    df["pastSwabDate"] = np.where(reportedDate, True, False)
    
    return df

def processDfDelta(df, timedelta = 1, reactDates=False):

    #df["timedelta_swabs"] = (df["d_comb"] - df["agprev5_new"]).dt.days

    if reactDates == False:
        reinfectionCriteria = \
            [
                (df["react_lineage"] == "BA.1.617.2") & \
                ((df["agprev5_new"] - pd.to_datetime("2021-03-18")).dt.days >= timedelta)
            ]
    else:
        reinfectionCriteria = \
            [
                (df["react_lineage"] == "BA.1.617.2") & \
                ((df["agprev5_new"] - pd.to_datetime("2021-04-27")).dt.days >= timedelta)
            ]

    reinfectionChoices = [True]
            
    reportedDate = pd.notnull(df["agprev5_new"])

    df["reinfection"] = np.select(reinfectionCriteria,
                                         reinfectionChoices,
                                         default=False)

    df["pastSwabDate"] = np.where(reportedDate, True, False)
    
    return df

def getLineage(df):
    return pd.concat(
    [(df[["reinfection", "pastSwabDate", "react_lineage"]]
        .groupby("react_lineage")
        .sum()),
    (df[["estbinres", "react_lineage"]]
         .groupby("react_lineage")
         .count()        
         .rename(columns={"estbinres":"total"}))],
     axis = 1)

def getChiSq(df):
    df = getLineage(df)
    p_val = stats.chi2_contingency(pd.concat(
    [df["reinfection"], df["total"]-df["reinfection"]], axis=1))[1]

    out = pd.DataFrame([p_val], columns=["chi2 p-val"])            
    return out

df18E = processDfOmicron(df18)
df13E = processDfDelta(df13)

df18R = processDfOmicron(df18, reactDates=True)
df13R = processDfDelta(df13, reactDates=True)

df18R2 = processDfOmicron(df18, reactDates=True, timedelta=30)
df13R2 = processDfDelta(df13, reactDates=True, timedelta=30)


writer = pd.ExcelWriter("E:/custom_marc/tables/R1718reinfections.xlsx",
                        engine="xlsxwriter")

getLineage(df18E).to_excel(writer)
getChiSq(df18E).to_excel(writer, startrow = 0, startcol=4, index=False)

getLineage(df18R).to_excel(writer, startcol = 6)
getChiSq(df18R).to_excel(writer, startrow = 0, startcol=10, index=False)

getLineage(df18R2).to_excel(writer, startcol = 12)
getChiSq(df18R2).to_excel(writer, startrow = 0, startcol=16, index=False)


##getLineage(df13E).to_excel(writer, startrow = 6)
##getLineage(df13R).to_excel(writer, startrow = 6, startcol = 6)
##getLineage(df13R).to_excel(writer, startrow = 6, startcol = 12)
# getChiSq(df13).to_excel(writer, startrow = 6, startcol=4, index=False)

writer.save()


