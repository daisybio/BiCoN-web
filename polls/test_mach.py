import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

# Load the diabetes dataset
diabetes_data = datasets.load_diabetes()

# Print all keys and number of raw and columns
print(diabetes_data.keys,diabetes_data.data.shape)

# Print all the feature_names in dataset
print(diabetes_data.feature_names)

di = pd.DataFrame(diabetes_data.data)
di.columns = diabetes_data.feature_names

di['target'] = diabetes_data.target

x=di.drop('target',axis=1)

# Create linear regression object
rm = linear_model.LinearRegression()

rm.fit(x,di.target)

print(rm.intercept_)
print(rm.coef_)

print(rm.predict(x)[:10])

plt.scatter(di.target,rm.predict(x))
plt.xlabel('old data')
plt.ylabel('predicted data')
plt.show()
#mse = np.mean(di.target — rm(x)*2)
#print(mse)
