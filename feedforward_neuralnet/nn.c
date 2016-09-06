//
//  nn_model.c
//  Neural_Network
//
//  Created by Gurpreet Gosal on 2016-07-15.
//  Copyright © 2016 uWaterloo. All rights reserved.
//

/*
 *  Multilayer network using the backpropagation learning rule.
 *
 *  Compile using: gcc nn.c -lm
 */
//#include "nn_model.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/*
 *  When determining classification errors during training,
 *  if the difference between actual output and desired output
 *  is less than or equal to this, consider it close enough to accept.
 *  YOU MAY WISH TO CHANGE THIS.
 *  Valid range: 0.00, ..., 1.00
 */
#define ACCEPT_TRAIN 0.2

/*
 *  When determining classification errors during testing,
 *  if the difference between actual output and desired output
 *  is less than or equal to this, consider it close enough to accept.
 *  NO NEED TO CHANGE THIS.
 */
#define ACCEPT_TEST   0.49999

/*
 *  Max number of epochs in training.
 *  YOU MAY WISH TO CHANGE THIS.
 *  Valid range: 1, ...
 */
#define MAX_EPOCH  10000

/*
 *  Min errors to terminate training.
 *  YOU MAY WISH TO CHANGE THIS.
 *  Valid range: 0, ...
 */
#define MIN_ERRORS  5

/*
 *  Fraction of the data to use for training.
 *  The rest is used for testing.
 *  YOU MAY WISH TO CHANGE THIS.
 *  Valid range: 0.00, ..., 1.00
 */
#define TRAINING_PROPORTION  0.7
#define VALIDATION_PROPORTION 0.0
#define TEST_PROPORTION 0.3

/*
 *  Max number of units in input, hidden, and output layers
 *  NO NEED TO CHANGE THESE.
 */
#define MAX_INPUT    70
#define MAX_HIDDEN   70
#define MAX_OUTPUT   16

/*
 *  Max number of instances/examples.
 *  HOPEFULLY NO NEED TO CHANGE THIS.
 */
#define N_EXAMPLES 2823


typedef struct {
    int    n_input;
    int    n_hidden;
    int    n_output;
    double h[MAX_HIDDEN];
    double o[MAX_OUTPUT];
    double W1[MAX_INPUT][MAX_HIDDEN];
    double W2[MAX_HIDDEN][MAX_OUTPUT];
    double learning_rate;
    double k;
} nn_type;


typedef struct {
    int    n;
    double x[N_EXAMPLES][MAX_INPUT];
    double y[N_EXAMPLES][MAX_OUTPUT];
} dd_type;


void    readData(  int n_input, int n_output, dd_type *all_data );
void    splitData( int n_input, int n_output, dd_type *all_data,
                  dd_type *training_data, dd_type *test_data );
void    trainNetwork( FILE *log_file, nn_type *nn, dd_type *training_data );
void    testNetwork(  FILE *out_file, nn_type *nn, dd_type *test_data );
void    printNetwork( FILE *out_file, nn_type *nn );

dd_type all_data;
dd_type training_data;
dd_type validation_data;
dd_type test_data;

int
main( int argc, char *argv[] )
{
    nn_type nn;
    double  atof();
    FILE    *log_file;
    FILE    *out_file;
    
    if( argc != 6 ) {
        fprintf( stderr, "Usage: nn learning_rate k hidden log_file out_file < digits_train.txt\n" );
        fprintf( stderr, "       log_file - file to record progress of training\n");
        fprintf( stderr, "       out_file - file to record final network\n");
        exit(0);
    }
    
    nn.learning_rate = atof( argv[1] );
    nn.k             = atof( argv[2] );
    nn.n_hidden      = atoi( argv[3] );
    
    if( (log_file = fopen( argv[4], "w" )) == NULL ) {
        fprintf( stderr, "Could not open file %s\n", argv[4] );
        exit( 0 );
    }
    
    if( (out_file = fopen( argv[5], "w" )) == NULL ) {
        fprintf( stderr, "Could not open file %s\n", argv[5] );
        exit( 0 );
    }
    
    fprintf( log_file, "learning rate: %0.2f\n", nn.learning_rate );
    fprintf( log_file, "multiplicative constant (k): %0.1f\n", nn.k );
    fprintf( log_file, "hidden units: %d\n", nn.n_hidden );
    
    /*
     *  Number of input lines.
     *  NO NEED TO CHANGE THIS.
     */
    nn.n_input = 64;
    
    /*
     *  Number of output lines.
     *  NO NEED TO CHANGE THIS.
     */
    nn.n_output = 10;
    
    /*
     *  Total amount of data.
     *  YOU MAY WISH TO CHANGE THIS.
     */
    all_data.n = N_EXAMPLES;  /* total amount of data */
    if( all_data.n > N_EXAMPLES ) {
        fprintf( stderr, "Too many examples; increase N_EXAMPLES\n" );
        exit( 0 );
    }
    
    readData( nn.n_input, nn.n_output, &all_data );
    splitData( nn.n_input, nn.n_output,
              &all_data, &training_data, &test_data );
    trainNetwork( log_file, &nn, &training_data );
    testNetwork(  log_file, &nn, &test_data );
    printNetwork( out_file, &nn );
    
    fclose( log_file );
    fclose( out_file );
}


/*
 *  Read in the data.
 */
void
readData( int n_input, int n_output, dd_type *data )
{
    int  i, j;
    
    /*
     *  Read in the training examples.
     0,1,6,15,12,1,0,0,0,7,16,6,6,10,0,0,0,8,16,2,0,11,2,0,0,5,16,3,0,5,7,0,0,7,13,3,0,8,7,0,0,4,12,0,1,13,5,0,0,0,14,9,15,9,0,0,0,0,6,14,7,1,0,0,0
     */
    for( i = 0; i < data->n; i++ ) {
        
        /*
         *  Read in feature values of training example.
         *  Normalize to the range [-1, 1].
         */
        int value;
        for( j = 1; j <= n_input; j++ ) {
            scanf( "%d,", &value );
            data->x[i][j] = ((double)value/8.0) - 1.0;
        }
        
        /*
         *  Read in class value of training example and
         *  set desired outputs using a unary encoding.
         */
        char c;
        scanf( "%c\n", &c );
        data->y[i][1] = 0;
        data->y[i][2] = 0;
        data->y[i][3] = 0;
        data->y[i][4] = 0;
        data->y[i][5] = 0;
        data->y[i][6] = 0;
        data->y[i][7] = 0;
        data->y[i][8] = 0;
        data->y[i][9] = 0;
        data->y[i][10] = 0;
        switch (c) {
            case '0':
                data->y[i][1] = 1;
                break;
            case '1':
                data->y[i][2] = 1;
                break;
            case '2':
                data->y[i][3] = 1;
                break;
            case '3':
                data->y[i][4] = 1;
                break;
            case '4':
                data->y[i][5] = 1;
                break;
            case '5':
                data->y[i][6] = 1;
                break;
            case '6':
                data->y[i][7] = 1;
                break;
            case '7':
                data->y[i][8] = 1;
                break;
            case '8':
                data->y[i][9] = 1;
                break;
            case '9':
                data->y[i][10] = 1;
                break;
            default:
                fprintf( stderr, "error in input\n" );
                exit( 6 );
                break;
        }
        //printf("\n x[1][1]=%d\ty[1]=%15.10f", data->x[0][3],data->y[0][1]);
        //getchar();
    }

}


/*
 *  Split the data into training set and test set.
 */
void
splitData( int n_input,
          int n_output,
          dd_type *all_data,
          dd_type *training_data,
          dd_type *test_data )
{
    int  i, j, k, l, remove;
    
    srandom( time(0) );
    
    /*
     *  Use % of the data for training.
     */
    training_data->n = TRAINING_PROPORTION * all_data->n;
    remove = all_data->n - training_data->n;
    test_data->n = remove;
    for( i = 1; i <= remove;  ) {
        k = random() % all_data->n;
        /* mark k^th example as not in training data if it's not
         already marked */
        if( training_data->x[k][0] != -1.0 ) {
            training_data->x[k][0] = -1.0;
            i++;
        }
    }
    
    k = 0;
    l = 0;
    for( i = 0; i < all_data->n; i++ ) {
        if( training_data->x[i][0] != -1.0 ) {
            for( j = 0; j <= n_input; j++ ) {
                training_data->x[k][j] = all_data->x[i][j];
            }
            for( j = 1; j <= n_output; j++ ) {
                training_data->y[k][j] = all_data->y[i][j];
            }
            k++;
        }
        else {
            for( j = 0; j <= n_input; j++ ) {
                test_data->x[l][j] = all_data->x[i][j];
            }
            for( j = 1; j <= n_output; j++ ) {
                test_data->y[l][j] = all_data->y[i][j];
            }
            l++;
        }
    }
}


/*
 *  Sigmoid function.
 */
double
f( double k, double x )
{
    return( 1.0 / (1.0 + exp(-k * x)) );
    
}


/*
 *  Determine the actual output of the network
 *  for the given example x.
 */
void
determineOutput( nn_type *nn, double x[] )
{
    double  t, f();
    int     i, j;
    
    /*
     *  Determine the activation levels of the units in
     *  the hidden layer.
     */
    x[0] = 1.0;
    for( j = 1; j <= nn->n_hidden; j++ ) {
        t = 0.0;
        for( i = 0; i <= nn->n_input; i++ ) {
            t += nn->W1[i][j] * x[i];
        }
        nn->h[j] = f( nn->k, t );
    }
    
    /*
     *  Determine the activation levels of the units in
     *  the output layer.
     */
    nn->h[0] = 1.0;
    for( j = 1; j <= nn->n_output; j++ ) {
        t = 0.0;
        for( i = 0; i <= nn->n_hidden; i++ )
            t += nn->W2[i][j] * nn->h[i];
        nn->o[j] = f( nn->k, t );
    }
}


/*
 *  Determine how to adjust weights to reduce sum
 *  of squared error for this training example.
 */
void
adjustWeights( nn_type *nn, double x[], double y[] )
{
    int     i, j;
    double  t;
    double  E1[MAX_HIDDEN], E2[MAX_OUTPUT];
    
    /*
     *  Determine how to adjust weights between hidden and
     *  output layer to reduce sum of squared error for
     *  this training example.
     */
    for( j = 1; j <= nn->n_output; j++ ) {
        E2[j] = nn->k * nn->o[j]*(1.0-nn->o[j])*(y[j]-nn->o[j]);
    }
    
    /*
     *  Determine how to adjust weights between input and
     *  hidden layer to reduce sum of squared error for
     *  this training example.
     */
    for( j = 1; j <= nn->n_hidden; j++ ) {
        t = 0.0;
        for( i = 1; i <= nn->n_output; i++ )
            t += E2[i] * nn->W2[j][i];
        E1[j] = nn->k * nn->h[j]*(1.0-nn->h[j])*t;
    }
    
    /*
     *  Adjust weights between hidden and output layer.
     */
    for( i = 0; i <= nn->n_hidden; i++ )
        for( j = 1; j <= nn->n_output; j++ ) {
            nn->W2[i][j] = nn->W2[i][j] + nn->learning_rate*E2[j]*nn->h[i];
        }
    
    /*
     *  Adjust weights between input and hidden layer.
     */
    for( i = 0; i <= nn->n_input; i++ )
        for( j = 1; j <= nn->n_hidden; j++ ) {
            nn->W1[i][j] = nn->W1[i][j] + nn->learning_rate*E1[j]*x[i];
        }
}


void
trainNetwork( FILE *log_file, nn_type *nn, dd_type *data )
{
    int     i, j, p;
    int     epoch, errors;
    double  drand48();
    void    determineOutput();
    void    adjustWeights();
    
    srand48( time(0) );
    
    /*
     *  Initialize weights and thresholds to small random values in
     *  the range (-0.5, 0.5). I use the convention that -W1[0][j]
     *  and -W2[0][j] are the thresholds for the hidden units and
     *  the output units, respectively.
     */
    for( i = 0; i <= nn->n_input; i++ )
        for( j = 1; j <= nn->n_hidden; j++ ) {
            nn->W1[i][j] = drand48() - 0.5;
        }
    
    for( i = 0; i <= nn->n_hidden; i++ )
        for( j = 1; j <= nn->n_output; j++ ) {
            nn->W2[i][j] = drand48() - 0.5;
        }
    
    nn->h[0] = 1.0;
    
    epoch = 0;
    do {
        epoch++;
        errors = 0;
        
        /*
         *  Loop through training data (one epoch).
         */
        for( p = 0; p < data->n; p++ ) {
            
            /*
             *  Determine the actual output of the network
             *  for the given example x[p].
             */
            determineOutput( nn, data->x[p] );
            
            /*
             *  Determine the number of classification errors by
             *  comparing the desired output with the actual output.
             *  YOU MAY WISH TO CHANGE THIS.
             */
            for( j = 1; j <= nn->n_output; j++ ) {
                if( fabs( nn->o[j] - data->y[p][j] ) > ACCEPT_TRAIN ) {
                    errors++;
                    break;
                }
            }
            
            /*
             *  Determine how to adjust weights to reduce sum
             *  of squared error for the given example x[p], y[p].
             */
            if( errors ) {
                adjustWeights( nn, data->x[p], data->y[p] );
            }
        }
        
        fprintf( log_file, "Epoch: %3d,  ", epoch );
        fprintf( log_file, "Number of classification errors: %d (/%d)\n",
                errors, data->n );
        fflush( log_file );
        
        /*
         *  Stopping criteria.
         *  YOU MAY WISH TO CHANGE THIS.
         */
    } while ((errors > MIN_ERRORS) && (epoch < MAX_EPOCH));
    
    fprintf( log_file, "Training epochs = %5d\n", epoch );
    fprintf( log_file, "Classification errors on training data = %d (/%d)\n",
            errors, data->n );
    fflush( log_file );
}

void
testNetwork( FILE *log_file, nn_type *nn, dd_type *data )
{
    int     j, p;
    int     errors;
    void    determineOutput();
    
    errors = 0;
    for( p = 0; p < data->n; p++ ) {
        
        determineOutput( nn, data->x[p] );
        
        /*
         *  Determine the number of errors by comparing
         *  the desired output with the actual output.
         */
        for( j = 1; j <= nn->n_output; j++ ) {
            if( fabs( nn->o[j] - data->y[p][j] ) > ACCEPT_TEST ) {
                errors++;
                break;
            }
        }
    }
    
    fprintf( log_file, "Classification errors on test data = %d (/%d)\n",
            errors, data->n );
}


void
printNetwork( FILE *out_file, nn_type *nn )
{
    int     i, j;
    
    fprintf( out_file, "nn.n_input = %d;\n", nn->n_input );
    fprintf( out_file, "nn.n_hidden = %d;\n", nn->n_hidden );
    fprintf( out_file, "nn.n_output = %d;\n", nn->n_output );
    fprintf( out_file, "nn.learning_rate = %0.3f;\n", nn->learning_rate );
    fprintf( out_file, "nn.k = %0.3f;\n", nn->k );
    
    /*
     *  Print out the final weights and threshold values.
     */
    fprintf( out_file, "// Hidden layer:\n" );
    for( j = 1; j <= nn->n_hidden; j++ )
        for( i = 0; i <= nn->n_input; i++ ) {
            fprintf( out_file, "\tnn.W1[%d][%d] = %6.4f;\n", i, j, nn->W1[i][j] );
        }
    
    fprintf( out_file, "// Output layer:\n" );
    for( j = 1; j <= nn->n_output; j++ )
        for( i = 0; i <= nn->n_hidden; i++ ) {
            fprintf( out_file, "\tnn.W2[%d][%d] = %6.4f;\n", i, j, nn->W2[i][j] );
        }
}

