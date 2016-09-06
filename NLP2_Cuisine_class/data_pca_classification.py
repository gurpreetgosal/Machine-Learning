# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 23:01:08 2015

@author: z6liao, nbraniff,ggosal
"""
import numpy as np
from nltk.stem import WordNetLemmatizer
# Need to download wordnet package. Type nltk.download() in Python
from sklearn.feature_extraction.text import TfidfVectorizer, TfidfTransformer
from sklearn.ensemble import AdaBoostClassifier
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm, grid_search
from sklearn.decomposition import PCA
from sklearn.cross_validation import train_test_split
from sklearn.decomposition import RandomizedPCA
from sklearn import cross_validation
import pandas as pd
import re
import scipy.io as sio
train = pd.read_json("train.json")

y = train['cuisine']
train['ingred_l_noun'] = [' '.join([' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', word)) for word in line.split()]).strip() for line in lists]).strip() for lists in train['ingredients']]
tv = TfidfVectorizer(stop_words='english', ngram_range=(1 , 1), analyzer="word", token_pattern=r'\w+')
x = tv.fit_transform(train['ingred_l_noun'])

colnames=tv.get_feature_names()

colnms=tuple([a.encode('UTF8') for a in colnames])

#%% Export the data to Matlab for SPCA

st=list(set(y))
ynum=[st.index(a) for a in y]
sio.savemat("dat.mat",{'x': x})
sio.savemat("ynum.mat",{'ynum': ynum})
sio.savemat("cuisines",{'cuisines':st})
sio.savemat("colnms",{'colnms':colnms})

#%% Load the SPCA results from Matlab

spca = sio.loadmat('spca_2.mat')

Y_tr=spca['Y_tr']
Y_ts=spca['Y_ts']
Ztr_SPCA=spca['Ztr_SPCA']
Zts_SPCA=spca['Zts_SPCA']

Y_tr=Y_tr.flatten()
Y_ts=Y_ts.flatten()

#%% Build SPCA SVM, score on test set (plit happens in matlab)

lsvm = svm.LinearSVC()
lsvm.fit(Ztr_SPCA, Y_tr)

#mboost = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1), n_estimators=5000, learning_rate=1, algorithm='SAMME.R', random_state=True)
#mboost.fit(X_trn_pca, y_train)

res=lsvm.score(Zts_SPCA,Y_ts)


#%% Apply PCA to whole dataset to and test on Kaggle test set for submission to Kaggle
temp=x.todense()
#pca = PCA(n_components=2000)
#X_pca = pca.fit(temp).transform(temp)

pca = RandomizedPCA(n_components=2000)
X_pca = pca.fit(temp).transform(temp)

lsvm = svm.LinearSVC()
lsvm.fit(X_pca, y)

#model = lsvm
test = pd.read_json("test.json")
#test['ingred_l'] = [' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', line)) for line in lists]).strip() for lists in test['ingredients']]
test['ingred_l_noun'] = [' '.join([' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', word)) for word in line.split()]).strip() for line in lists]).strip() for lists in test['ingredients']]

x_test = tv.transform(test['ingred_l_noun'])
X_tst_pca=pca.transform(x_test.todense())
##test['cuisine'] = model.predict(x_test)
test['cuisine'] = lsvm.predict(X_tst_pca)
test[['id', 'cuisine']].to_csv("result_pca.csv", index=False)

#%%Eval linear SVM (or Adaboost which is commented out) on 25% test set

X_train, X_test, y_train, y_test = train_test_split(train['ingred_l_noun'], y, test_size=0.25, random_state=41)

trnV = TfidfVectorizer(stop_words='english', ngram_range=(1 , 1), analyzer="word", token_pattern=r'\w+')
X_tf_trn = trnV.fit_transform(X_train)

#temp=X_tf_trn.todense()
#pca = RandomizedPCA(n_components=250)
#X_trn_pca = pca.fit(temp).transform(temp)

#lsvm = svm.LinearSVC()
#lsvm.fit(X_trn_pca, y_train)

#lsvm = svm.SVC(probability=True,kernel='linear')
#lsvm.fit(X_trn_pca, y_train)

lsvm = svm.SVC(probability=True,kernel='linear')
lsvm.fit(X_tf_trn, y_train)

#mboost = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1), n_estimators=5000, learning_rate=1, algorithm='SAMME.R', random_state=True)
#mboost.fit(X_trn_pca, y_train)

X_tf_tst = trnV.transform(X_test)
#X_tst_pca=pca.transform(X_tf_tst.todense())

#res=lsvm.score(X_tst_pca,y_test)
#pred=lsvm.predict(X_tst_pca)

res=lsvm.score(X_tf_tst,y_test)
pred=lsvm.predict(X_tf_tst)

#res=mboost.score(X_tst_pca,y_test)
#pred=mboost.predict(X_tst_pca)

wrong=pd.value_counts(y_test[y_test!=pred])
tot=pd.value_counts(y_test)

frac=wrong/tot[wrong.axes[0]]

prec=tot[wrong.axes[0]]/sum(tot)

#%% Generate error distribution across classes
n_groups = 20
plt.figure(1)
fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
error_config = {'ecolor': '0.3'}

rects1 = plt.bar(index, tot, bar_width,
                 alpha=opacity,
                 color='b',
                 error_kw=error_config,
                 label='# Recipes')

rects2 = plt.bar(index + bar_width, wrong[tot.axes[0]], bar_width,
                 alpha=opacity,
                 color='r',
                 error_kw=error_config,
                 label='Incorrect')

nms=tuple([a.encode('UTF8') for a in tot.axes[0]])
#plt.xlabel('Group')    tv.get_feature_names()
#plt.ylabel('Scores')
#$plt.title('Scores by group and gender')
locs, labels=plt.xticks(index + bar_width, nms)
plt.legend()
plt.setp(labels, rotation=90)
plt.gcf().subplots_adjust(bottom=0.25)




#%%PCA plot, used for plotting PCA results but mainly used R for this

#reds = y_train == -1
#blues = y_train == 1
#for i in range(5):
#    plt.figure()
#    plt.plot(X_trn_pca[reds, i], X_trn_pca[reds, i+1], "ro")
#    plt.plot(X_trn_pca[blues, i], X_trn_pca[blues, i+1], "bo")
#    plt.xlabel(i+1)
#    plt.ylabel(i+2)

#%%

#BELOW: mosty old code I used for testing various models on the full dataset for submission to Kaggle

#%% CV LSVM

#clf = svm.LinearSVC()
#scores = cross_validation.cross_val_score(clf, x, y, cv=5)

#%% SVM
#
#lsvm = svm.LinearSVC()
#lsvm.fit(X_pca, y)

#%% SVM
#
#lsvm = svm.SVC(probability=True,kernel='linear')
#lsvm.fit(X_pca, y)

#%% Multiclass Adaboost

#mboost = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1), n_estimators=10000, learning_rate=1, algorithm='SAMME.R', random_state=True)
#mboost.fit(x, y)

#%% Multiclass Adaboost SVM

#mboost = AdaBoostClassifier(svm.SVC(probability=True,kernel='linear'), n_estimators=10, learning_rate=1, algorithm='SAMME', random_state=True)
#mboost.fit(x, y)

#%% Write Results to File
#
#model = lsvm
#test = pd.read_json("test.json")
#test['ingred_l'] = [' '.join([WordNetLemmatizer().lemmatize(re.sub('[^A-Za-z]', ' ', line)) for line in lists]).strip() for lists in test['ingredients']]
#x_test = tv.transform(test['ingred_l'])
#X_tst_pca=pca.transform(x_test.todense())
##test['cuisine'] = model.predict(x_test)
#test['cuisine'] = model.predict(X_tst_pca)
#test[['id', 'cuisine']].to_csv("result_pca.csv", index=False)

#%%

f=open('f1.txt','w') 
for ele in nms: 
    f.write(ele+'\n') 
f.close()
