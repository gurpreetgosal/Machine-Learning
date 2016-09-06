
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
from sklearn.ensemble import AdaBoostClassifier


# Load .mat file which contains document-term dictionary
mat_contents = sio.loadmat('ads-data.mat')
print '\n number of advertisements to be classified= ', np.size(mat_contents['train_classes'])

labels = mat_contents['train_classes']
labels = labels.ravel()
data = mat_contents['train_mat']

# We further divide the training data in to train data, validation data and test data

train_share = 0.7
valid_share = (1-train_share)/2
test_share  = valid_share

#train_data =
# Build a Pipeline for vectorizer -> tfidf-transformer ->  classifier
'''text_clf = Pipeline([('vect', CountVectorizer()), ('tfidf', TfidfTransformer()),
                     ('clf', MultinomialNB())]) '''

tfidf_transformer = TfidfTransformer()
data_tfidf = tfidf_transformer.fit_transform(data)


# Train and validate Classifier
#predicted = clf.predict(test_data_tfidf)
#clf = MultinomialNB().fit(train_data_tfidf, train_labels)
#clf = MultinomialNB(alpha = 0.00001)
'''clf = SGDClassifier(loss='modified_huber', penalty='l2', class_weight = 'balanced',
                                         alpha=1e-4, n_iter=20, random_state=42)'''
clf = svm.SVC(C=1, cache_size=2000, class_weight= 'balanced', coef0=8.3,
  decision_function_shape=None, degree=77, gamma=1.0, kernel='linear',
  max_iter=-1, probability=False, random_state=None, shrinking=True,
  tol=0.004, verbose=False)

abSVM = AdaBoostClassifier(clf,
                         algorithm="SAMME",
                         n_estimators=200)

scores = cross_validation.cross_val_score(abSVM, data_tfidf, labels, cv=5, scoring='f1_weighted')

'''
parameters = {'vect__ngram_range': [(1, 1), (1, 2)],
              'tfidf__use_idf': (True, False),
              'clf__alpha': (1e-2, 1e-3)}
              '''
parameters = {
    'vec__max_df': (0.5, 0.625, 0.75, 0.875, 1.0),
    'vec__max_features': (None, 5000, 10000, 20000),
    'vec__min_df': (1, 5, 10, 20, 50),
    'tfidf__use_idf': (True, False),
    'tfidf__sublinear_tf': (True, False),
    'vec__binary': (True, False),
    'tfidf__norm': ('l1', 'l2'),
    'clf__alpha': (1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
    }

gs_clf = GridSearchCV(clf, parameters)
gs_clf = gs_clf.fit(data, labels)
best_parameters, score, _ = max(gs_clf.grid_scores_, key=lambda x: x[1])

for param_name in sorted(parameters.keys()):
    print("%s: %r" % (param_name, best_parameters[param_name]))



# now test the classifier performance on test data
test_data = mat_contents['test_mat']
test_data_tfidf = tfidf_transformer.fit_transform(data)

gg=2



