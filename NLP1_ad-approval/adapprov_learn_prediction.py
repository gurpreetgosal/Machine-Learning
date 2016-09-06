
import numpy as np
import scipy as sp
import scipy.io as sio
import matplotlib as plt
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.naive_bayes import BernoulliNB, MultinomialNB
from sklearn.pipeline import Pipeline
from sklearn import cross_validation
from sklearn.grid_search import GridSearchCV
from sklearn import svm
from sklearn.svm import LinearSVC
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model import Perceptron
from sklearn.linear_model import PassiveAggressiveClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.extmath import density
from sklearn import metrics
import pandas as pd


# Load .mat file which contains document-term dictionary
mat_contents = sio.loadmat('ads-data.mat')
print '\n number of advertisements to be classified= ', np.size(mat_contents['train_classes'])

labels = mat_contents['train_classes']
labels = labels.ravel()
data = mat_contents['train_mat']
test_data = mat_contents['test_mat']

# We further divide the training data into train data, validation data and test data

#train_data =
# Build a Pipeline for vectorizer -> tfidf-transformer ->  classifier
'''text_clf = Pipeline([('vect', CountVectorizer()), ('tfidf', TfidfTransformer()),
                     ('clf', MultinomialNB())]) '''

tfidf_transformer = TfidfTransformer()
data_tfidf = tfidf_transformer.fit_transform(data)
test_data_tfidf = tfidf_transformer.transform(test_data)

# Train and validate Classifier
#predicted = clf.predict(test_data_tfidf)
#clf = MultinomialNB().fit(train_data_tfidf, train_labels)
#clf = MultinomialNB(fit_prior=  'boolean')
clf = svm.SVC(C=1, cache_size=2000, class_weight= 'balanced', coef0=8.3,
  decision_function_shape=None, degree=77, gamma=1.0, kernel='linear',
  max_iter=-1, probability=False, random_state=None, shrinking=True,
  tol=0.004, verbose=False).fit(data_tfidf , labels)

#clf = MultinomialNB()
#clf = MultinomialNB().fit(data_tfidf , labels)

# now test the classifier performance on test data
predicted = clf.predict(test_data_tfidf)
ids = np.empty([len(predicted), 1])
ids[:,0] = range(1, len(predicted )+1)
predicted = predicted.reshape((len(predicted),1))

matrix1 = np.append(ids, predicted, axis =1)
columns = ['Id', 'prediction']

df = pd.DataFrame(matrix1, columns=columns)
df.to_csv('ads_pred_opt2.csv')

gg=2



