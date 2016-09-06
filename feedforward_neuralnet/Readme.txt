Multi-layered Feedforward Neural Network neural network to recognize hand-drawn digits. The input is a vector of features that represents a 2-D binary pattern, and the output will be a classification of the input as one of the ten possible digits,’0’,’1’,…’9’. 

For each input sample binary pattern is divided into cells of of 4x4 with total of 64 4x4 cells. The pixels within a cell which are on are counted as 1 whereas off pixels are counted as 0. Finally for each cells all ‘on’ pixels are added. This way we have feature vector of dimension 64.And we have 64 input neurons, hidden neurons are selected by experimentation and 10 output neurons.

nn.c contains the code which implements Backpropagation algorithm to train the weight parameters of neural network.

To run this code:

Compile:
gcc nn.c -lm

Run:
nn a k h log_file out_file < digits_train.txt

h = number of hidden neutrons, 
a = learning rate, 
k = is a constant used in sigmoid function that determines its slope and spread
log_file = file to output the training error after each epoch (log_network.txt)
out_file = the final weight parameters (network.h)
digits_train.txt = training data