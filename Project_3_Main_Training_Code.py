
# coding: utf-8

# In[14]:


import numpy as np
import cv2
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from itertools import product, combinations
import random
import scipy
import tensorflow as tf
import keras


from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Dropout, Flatten, Conv3D, MaxPool3D, BatchNormalization, Input, Concatenate
from tensorflow.keras.optimizers import RMSprop, Adam
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.callbacks import ReduceLROnPlateau, TensorBoard
from tensorflow.keras.callbacks import ModelCheckpoint

#import seaborn as sns
#sns.set_style('white')
from sklearn.metrics import confusion_matrix, accuracy_score


def Conv(filters=16, kernel_size=(3,3,3), activation='relu', input_shape=None):
    if input_shape:
        return Conv3D(filters=filters, kernel_size=kernel_size, padding='Same', activation=activation, input_shape=input_shape)
    else:
        return Conv3D(filters=filters, kernel_size=kernel_size, padding='Same', activation=activation)

# Define Model
def get_siamese_model(input_dim):
    
    # Define the tensors for the two input images
    left_input = Input(input_dim)
    right_input = Input(input_dim)
    
    model = Sequential()
    

    model.add(Conv(32, (3,3,3), input_shape=input_dim))
    model.add(Conv(32, (3,3,3)))
    
    # model.add(BatchNormalization())
    model.add(MaxPool3D())
    # model.add(Dropout(0.25))

    model.add(Conv(64, (3,3,3)))
    model.add(Conv(64, (3,3,3)))
    
    #model.add(BatchNormalization())
    model.add(MaxPool3D())
    #model.add(Dropout(0.25))
    
    model.add(Conv(128, (3,3,3)))
    model.add(Conv(128, (3,3,3)))

    model.add(Flatten())

    model.add(Dense(512, activation='relu'))
    
    #model.add(Dropout(0.5))

    #model.add(Dense(1024, activation='relu'))
    #model.add(Dropout(0.5))
    
    encoded_l = model(left_input)
    encoded_r = model(right_input)
    
    concatenated_layer= Concatenate()([encoded_l, encoded_r])
    

    final_layer = Dense(512,activation='relu')(concatenated_layer) 
    final_layer = Dense(512,activation='relu')(final_layer) 
    final_layer = Dense(2, activation = 'sigmoid')(final_layer)

    siamese_net = Model(inputs=[left_input,right_input],outputs=final_layer)
          
    return siamese_net


# In[18]:


model = get_siamese_model((20,20,20,4))
model.summary()
optimizer = RMSprop(lr = 0.0001)
model.compile(loss="binary_crossentropy",optimizer=optimizer,  metrics=['accuracy'])


# In[20]:


NNN = 160374


# In[21]:


targets = np.tile(np.array([1,0]), NNN//2)


# In[22]:


labels = dict(zip(  range(0,len(targets))  , targets.astype(int)))


# In[23]:


partition = {}


# In[27]:


partition['train'] = list(range(0,int(0.8*len(targets))))
partition['val'] = list(range(int(0.8*len(targets)),int(0.9*len(targets))))
partition['test'] = list(range(int(0.9*len(targets)),len(targets)))


# In[28]:


import numpy as np
import tensorflow as tf

class DataGenerator(tf.keras.utils.Sequence):
    'Generates data for Keras'
    def __init__(self, list_IDs, labels, batch_size=32, dim=(20,20,20), n_channels=4,
                 n_classes=2, shuffle=True):
        'Initialization'
        self.dim = dim
        self.batch_size = batch_size
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index*self.batch_size:(index+1)*self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        [X1,X2], y = self.__data_generation(list_IDs_temp)

        return [X1,X2], y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples' # X : (n_samples, *dim, n_channels)
        # Initialization
        X1 = np.empty((self.batch_size, *self.dim, self.n_channels))
        X2 = np.empty((self.batch_size, *self.dim, self.n_channels))
        y = np.empty((self.batch_size), dtype=int)

        # Generate data
        for i, ID in enumerate(list_IDs_temp):
            # Store sample
            X1[i,] = np.load('/home/gxs372/protein/Stanford_Interface/data/p1/' + repr(ID) + '.npy')
            X2[i,] = np.load('/home/gxs372/protein/Stanford_Interface/data/p2/' + repr(ID) + '.npy')

            # Store class
            y[i] = self.labels[ID]

        return [X1,X2], keras.utils.to_categorical(y, num_classes=self.n_classes)


# In[29]:



#from my_classes import DataGenerator

# Parameters
params = {'dim': (20,20,20),
          'batch_size': 32,
          'n_classes': 2,
          'n_channels': 4,
          'shuffle': True}

# Datasets
#partition = # IDs
#labels = # Labels

# Generators
training_generator = DataGenerator(partition['train'], labels, **params)
validation_generator = DataGenerator(partition['val'], labels, **params)
testing_generator = DataGenerator(partition['test'], labels, **params)


# In[ ]:


filepath="weights.best.hdf5"
checkpoint = ModelCheckpoint(filepath, monitor='val_accuracy', verbose=1, save_best_only=True, mode='max')
callbacks_list = [checkpoint]
#history = model.fit(train_data, train_label, validation_data=(val_data,val_label) ,epochs=10, batch_size=32, shuffle=True,verbose=1, callbacks=callbacks_list)


# Train model on dataset
history = model.fit_generator(generator=training_generator, validation_data=validation_generator, epochs=25, callbacks=callbacks_list, verbose = 0)


# In[14]:

res_test = model.evaluate_generator(generator=testing_generator) 

print(res_test)

model.save('/home/gxs372/protein/Stanford_Interface/model_1') 

# In[ ]:

