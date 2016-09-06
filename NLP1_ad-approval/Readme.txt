For this classification task, we are given 2500 labelled instances of web advertisements that appeared on a dozen websites that were related to farms, which have been labelled to indicate if they were approved or rejected by the website owners. Our job is to create a classifier that can try to predict if the owners of these websites would approve of new advertisements based on this information about their preferences.


 Implemented tf-idf vectorizer and trained on kernelled SVM as shown in:

adapprov_learn.py  (parameter selection for SVM using training data ads-data-train.txt)
adapprov_learn_prediction.py  (for final predictions on test data ads-data-test.txt)

Predictions are stored in ads-submission-sample.csv