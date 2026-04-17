# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import datetime
import os
import re
import glob
import os

#%%

absolute_path = os.path.dirname(__file__)
#%%
profile_df = pd.read_csv(absolute_path + '\\twitter_data\\Users_profile_data.csv')
profile_df = profile_df.dropna()

#%%
file = absolute_path + '\\twitter_data\\Tweets_of_4776users_CSVs\\'
csv_files = glob.glob(os.path.join(file, "*.csv"))
df_tweets = pd.DataFrame()

for f in csv_files:
    temp = pd.read_csv(f, encoding='latin-1')
    df_tweets = df_tweets.append(temp, ignore_index=True)

df_tweets = df_tweets.add_suffix('_tweeted_data')
print(df_tweets)
#%%
df_tweets['UserID_tweeted_data'] = df_tweets['UserID_tweeted_data'].str.strip('@')
#%%
def parser(fn):
    txt = open(fn).read()
    preparse = re.findall('"\d+",\d+,".+?","[^\t\n\r\f\v]+?","[^\t\n\r\f\v]+?"', txt, re.DOTALL)
    parsedRows = []
    for line in preparse:
        columns = line.split(',')
        output = {}

        output['Index'] = columns.pop(0)
        output['Tweet_ID'] = columns.pop(0)
        output['UserID'] = columns.pop()
        output['Company'] = columns.pop()
        output['Text'] = ','.join(columns)

        parsedRows.append(output)

    return parsedRows
            
    #parsed = [t.split(',') for t in preparse]
    #print(parsed)

#%%

data = parser(absolute_path + r'\twitter_data\EDITED - Reference_tweet_data.txt')
reference_df = pd.DataFrame(data)

reference_df['Index'] = reference_df['Index'].str.strip('"')
reference_df['Index'] = reference_df['Index'].astype(int)
reference_df['UserID'] = reference_df['UserID'].str.strip('"')
reference_df =reference_df.add_suffix('_ref_tweets')
print ("\nUnique values :  \n",reference_df.nunique())

duplicated = reference_df[reference_df.duplicated(subset='UserID_ref_tweets', keep=False)]
reference_without_duplicates = reference_df.drop_duplicates(
    subset = ['UserID_ref_tweets', 'Company_ref_tweets'],
    keep = 'last').reset_index(drop=True)
#%%
df = pd.merge(reference_without_duplicates, profile_df, left_on='UserID_ref_tweets', right_on='user_id')
print ("\nFeatures : \n" ,df.columns.tolist())
print ("\nMissing values :  ", df.isnull().sum().values.sum())
print ("\nUnique values :  \n",df.nunique())

#%%
finalDF = pd.merge(df_tweets,df, left_on='UserID_tweeted_data', right_on='UserID_ref_tweets')
print ("\nFeatures : \n" ,finalDF.columns.tolist())
print ("\nMissing values :  ", finalDF.isnull().sum().values.sum())
print ("\nUnique values :  \n",finalDF.nunique())

#%%
finalDF['text_tweeted_data'] = finalDF['text_tweeted_data'].astype(str)
finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].astype(str)
finalDF['Company_ref_tweets'] = finalDF['Company_ref_tweets'].str.strip('"')
finalDF['mention'] = finalDF.apply(lambda x: x.Company_ref_tweets in x.text_tweeted_data, axis=1)
print(finalDF['mention'].value_counts())

#%%
finalDF['mention_type'] = 2
finalDF.loc[finalDF['mention']==False,'mention_type'] = 1
print(finalDF['mention_type'].value_counts())

#%%
finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: pd.Timestamp(x.DateTime2_tweeted_data), axis=1)
finalDF['DateTime_tweeted_data'] = finalDF.apply(lambda x: x.DateTime_tweeted_data.timestamp(), axis=1)
#%%
"""WRITE FORMATTED DATA TO CSV"""
finalDF.to_csv(absolute_path + r'\final_output.csv')

#%%
"""IMPORT FORMATTED DATA"""
finalDF = pd.read_csv(r'C:\Users\Rob\OneDrive\NCSU PhD\twitter\final_output.csv')

print(finalDF.columns)

#%%
"""REMOVE EXTRA COLUMNS"""
final_column_list = ['DateTime2_tweeted_data', 'Index_ref_tweets', 'mention_type']
finalDF = finalDF[finalDF.columns.intersection(final_column_list)]
finalDF['mention_type'] = finalDF['mention_type'].astype(int)

#%%
"""lAST MONTH OF ACTIVITY"""

finalDF['DateTime2_tweeted_data'] = pd.to_datetime(finalDF['DateTime2_tweeted_data'])
month_ago = finalDF['DateTime2_tweeted_data'].max() - datetime.timedelta(days=30)

last_30_days_DF = finalDF.query('DateTime2_tweeted_data > @month_ago')
print(last_30_days_DF.head())

#%%

def round_time(dt: datetime, unit=20):
    #seconds = dt - dt.date()
    #unit_seconds = unit.total_seconds()
    #rounded_seconds = seconds - (seconds % unit_seconds)
    #return dt.date() + rounded_seconds
    return dt.replace(second=0, microsecond=0, minute=dt.minute-(dt.minute%unit), hour=dt.hour)
    
"""Bin Each User"""
each_user = last_30_days_DF.groupby('Index_ref_tweets')
output = pd.DataFrame()
counter = 0
start_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].min())
end_time=round_time(dt=last_30_days_DF['DateTime2_tweeted_data'].max())
print(start_time)
print(end_time)

t_index = pd.DatetimeIndex(pd.date_range(start=start_time, end=end_time, freq="20min"))

for user_name, user_data in each_user:
    print(f"Read {counter} out of {len(each_user)}")
    #df_bin=df_subset.set_index('DateTime2_tweeted_data',drop=False)
    binuser=user_data.groupby(pd.Grouper(key='DateTime2_tweeted_data', freq='20min')).max()
    binuser = binuser.reindex(t_index)
    binuser['mention_type'] = binuser['mention_type'].fillna(5)
    binuser['Index_ref_tweets'] = binuser['Index_ref_tweets'].fillna(user_name)
    binuser = binuser.transpose()
    print(binuser)
    output = output.append(binuser, ignore_index=True)
    counter += 1
    
#%%
"""SPLIT INTO THREE DATAFRAMES"""
"""Mentions"""
mentions = output.replace(to_replace=[1, 5, 2], value=[0, 0, 1])
mentions.to_csv(absolute_path + r'\mentions_dataset.csv')

no_mentions = output.replace(to_replace=[1, 5, 2], value=[1, 0, 0])
no_mentions.to_csv(absolute_path + r'\no_mentions_dataset.csv')

no_tweets = output.replace(to_replace=[1, 5, 2], value=[0, 1, 0])
no_tweets.to_csv(absolute_path + r'\no_tweet_dataset.csv')