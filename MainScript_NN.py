#Step 1: Import necessary libraries
#%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import tensorflow as tf
from tensorflow import keras
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt 
import random 


import csv

#Step 1.1: Path to saving and loading FF_NN model data.
LiftFlexible = True
PatuFlexible = False
LoadSaved_NN = False

if LiftFlexible:
    # Path to NN model
    keras_model_path  = 'ExData/NN_Data/LiftBoom_NN'    
    SimData           = np.genfromtxt('ExData/LiftFlexible/LiftBoomData.csv', 
                                      delimiter=",",dtype=np.float32)
    
elif PatuFlexible:
    # Path to NN model
    keras_model_path  = 'ExData/NN_Data/PatuFlexible_NN'
    SimData           = np.genfromtxt('ExData/LiftFlexible/PatuFlexible.csv', 
                                      delimiter=",",dtype=np.float32)    
    
else:
    # Path to NN model
    keras_model_path  = 'ExData/NN_Data/PatuRigid_NN'
    SimData           = np.genfromtxt('ExData/LiftFlexible/PatuRigid.csv', 
                                      delimiter=",",dtype=np.float32)      
      
    
# Path to NN model


if LoadSaved_NN:
    #Load previous keras model
    restored_keras_model = tf.keras.models.load_model(keras_model_path)
else:
    #Step 2: Loading and preparing the data
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # 2.1: Import dataset and split into train and test data
    n_samples, n_features = SimData.shape
    n_features           -= 1
    X                     = SimData[:,0:n_features]
    Y                     = SimData[:,n_features]
    
    
    print('Shape of Lift Boom data =',SimData.shape)
    print('Number of samples =',n_samples)
    print('Number of features =',n_features)
    print(X.shape, Y.shape)
    
    from sklearn.model_selection import train_test_split
    #X = SimData.iloc[:, :6]  # Use the first six columns as input features
    #Y = SimData.iloc[:, -1]  # Use the last column as the target variable
    
    # Normalize the data
    #from sklearn.preprocessing import MinMaxScaler
    #scaler              = MinMaxScaler()
    #X_scaled            = scaler.fit_transform(X)

    # Splitting data into training and testing sets
    from sklearn.model_selection import train_test_split

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=6)
    Y_train = Y_train.reshape(-1, 1)  # Reshape Y_train to have shape (n_samples, 1)
    Y_test = Y_test.reshape(-1, 1)  # Reshape Y_test to have shape (n_samples, 1) 
    
    # 2.2: Length of testing and training data
    print('Length of training dataset=',len(X_train))
    print('Length of testing dataset=',len(X_test))
    print('Shape of testing dataset=',X_test.shape)

    # 2.3: Data printing and data visualization
    #print(x_train[0])
    #plt.matshow(x_train[1])
    #plt.show()

    #Step 3: Create the model (Define the network architecture using keras).
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    input_dim  = n_features
    output_dim = 1
    
    model   = tf.keras.Sequential([tf.keras.layers.Dense(6, activation='relu', input_shape=(input_dim,)),
                                   tf.keras.layers.Dense(128, activation='relu'),
                                   tf.keras.layers.Dense(64, activation='relu'),
                                   tf.keras.layers.Dense(output_dim) # Use 'softmax' for multi-class classification
                            ])
    model.summary()

    #Step 4: Compiling and training the model (Train the model using sgd).  
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.06),
              loss='mean_squared_error',
              metrics=['accuracy'])

    history=model.fit(X_train,Y_train,validation_data=(X_test,Y_test), epochs=10)

    #Step 5: Evaluate the network and make predictions (Testing Performance).
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    test_loss, test_acc = model.evaluate(X_test,Y_test)
    print("Loss=%.3f" %test_loss)
    print("Accuracy=%.3f" %test_acc)


    n = random.randint(0,9999)
    plt.imshow(X_test[n])
    plt.show()

    predicted_value=model.predict(X_test)
    print("Handwritten number in the image is= %d" %np.argmax(predicted_value[n]))

    #Step 6: Training loss and accuracy.  
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    history.history.keys()

    #Figure 1
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.title('Model Accuracy')
    plt.ylabel('Accuracy')
    plt.xlabel('Epoch')
    plt.legend(['Train','Validation'],loc='upper left')
    plt.show()


    #Figure 2
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('Model Loss')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.legend(['Train','Validation'],loc='upper left')
    plt.show()


    #Figure 3
    plt.plot(history.history['accuracy'])
    plt.plot(history.history['val_accuracy'])
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('Training Loss And Accuracy')
    plt.ylabel('Accuracy/Loss')
    plt.xlabel('Epoch')
    plt.legend(['Accuracy','Val Accuracy','Loss','Val Loss'],loc='upper left')
    plt.show()

    #Step 7: Save model  
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    model.save(keras_model_path)
