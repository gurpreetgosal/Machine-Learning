# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:39:43 2015

@author: Gurpreet Gosal
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 23:01:08 2015

@author: ggosals
"""
# Gradient Tree Boosting Implementation using xgboost package

from nltk.stem import WordNetLemmatizer, porter
# Need to download wordnet package. Type nltk.download() in Python
from sklearn.feature_extraction.text import TfidfVectorizer, TfidfTransformer
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm, grid_search
import pandas as pd
import re
import numpy as np
import xgboost as xgb
from sklearn import preprocessing

train = pd.read_json("./train.json")
y = train['cuisine']
lencode = preprocessing.LabelEncoder()
lencode.fit(y)
list(lencode.classes_)
y=lencode.transform(y)



train['ingred_l_noun'] = [' '.join([' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', word)) for word in line.split()]).strip() for line in lists]).strip() for lists in train['ingredients']]

tv = TfidfVectorizer(stop_words='english', ngram_range=(1 , 1), analyzer="word", token_pattern=r'\w+')
x = tv.fit_transform(train['ingred_l_noun'])
dtrain = xgb.DMatrix(x, y)

# Ensemble technique: XGBoost

param = {}
# use softmax multi-class classification
param['objective'] = 'multi:softmax'
param['eta'] = 0.6
param['max_depth'] = 6
param['silent'] = 0
param['num_class'] = 20
param['nthread'] =1
param['eval_metric'] = 'auc'
param['max_delta_step'] =1

num_round = 280
bst = xgb.train(param, dtrain, num_round )

test = pd.read_json("./test.json")
#test['ingred_l'] = [' '.join([line for line in lists]).strip() for lists in test['ingredients']]
test['ingred_l_noun'] = [' '.join([' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', word)) for word in line.split()]).strip() for line in lists]).strip() for lists in test['ingredients']]

x_test = tv.transform(test['ingred_l_noun'])
dtest = xgb.DMatrix(x_test)

pred = bst.predict(dtest);
pred =lencode.inverse_transform(pred.astype(np.int64))
test['cuisine'] = pred

test[['id', 'cuisine']].to_csv("./result131.csv", index=False)
#%% Write Results to File

