import pandas as pd

df = pd.read_csv('timer.csv')

timer = df['time']

print(sum(timer)/1000)
