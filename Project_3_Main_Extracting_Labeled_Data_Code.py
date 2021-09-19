
# coding: utf-8

# In[1]:


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


# In[8]:


# Create a dictionary describing the features.
protein_feature_description = {
    'atom_names0': tf.io.VarLenFeature(tf.string), 
    'chains0': tf.io.VarLenFeature(tf.string),
    'elements0': tf.io.VarLenFeature(tf.string),
    'models0': tf.io.VarLenFeature(tf.string),
    'complex': tf.io.VarLenFeature(tf.string),
    'atom_names1': tf.io.VarLenFeature(tf.string), 
    'chains1': tf.io.VarLenFeature(tf.string),
    'elements1': tf.io.VarLenFeature(tf.string),
    'models1': tf.io.VarLenFeature(tf.string),
    'neg_idx': tf.io.VarLenFeature(tf.int64),
    'pos_idx': tf.io.VarLenFeature(tf.int64),
    'positions0': tf.io.VarLenFeature(tf.float32),
    'positions1': tf.io.VarLenFeature(tf.float32),
    'residues0': tf.io.VarLenFeature(tf.string), 
    'residues1': tf.io.VarLenFeature(tf.string),
    'resnames0': tf.io.VarLenFeature(tf.string),
    'resnames1': tf.io.VarLenFeature(tf.string),
    'pdb_names0': tf.io.VarLenFeature(tf.string),
    'pdb_names1': tf.io.VarLenFeature(tf.string),
}


def _parse_protein_function(example_protein):
  # Parse the input tf.train.Example proto using the dictionary above.
  return tf.io.parse_single_example(example_protein, protein_feature_description)


def pointdataFromRecord(record, starting_index, bl = 20):
    '''
    Parameters
    ----------         
    record: tfrecord 
    bl: surfacelet box length (Default = 20 Angs)
    
    
    Returns
    ----------         
    pairs = 
    targets = 
    
    '''
    #Atom-Names -
    atom_names0 = tf.sparse.to_dense(record['atom_names0'])
    atom_names1 = tf.sparse.to_dense(record['atom_names1'])
    atom_names0a = np.array(atom_names0)
    atom_names0b = [x.decode('utf-8') for x in atom_names0a]
    atom_names1a = np.array(atom_names1)
    atom_names1b = [x.decode('utf-8') for x in atom_names1a]

    #Elements - 
    elements0 = tf.sparse.to_dense(record['elements0'])
    elements1 = tf.sparse.to_dense(record['elements1'])
    elements0a = np.array(elements0)
    elements0b = [x.decode('utf-8') for x in elements0a]
    elements1a = np.array(elements1)
    elements1b = [x.decode('utf-8') for x in elements1a]

    #Positions - 
    positions0 = tf.sparse.to_dense(record['positions0'])
    positions1 = tf.sparse.to_dense(record['positions1'])
    positions0a = np.array(positions0)
    positions1a = np.array(positions1)
    x0 = positions0a[0::3]
    y0 = positions0a[1::3]
    z0 = positions0a[2::3]
    x1 = positions1a[0::3]
    y1 = positions1a[1::3]
    z1 = positions1a[2::3]
    xCA0 = []
    yCA0 = []
    zCA0 = []
    iCA0 = []

    for i in range(0,len(atom_names0b)):
        if (atom_names0b[i] == 'CA') :
            xCA0.append(x0[i])
            yCA0.append(y0[i])
            zCA0.append(z0[i])
            iCA0.append(i)


    xCA1 = []
    yCA1 = []
    zCA1 = []
    iCA1 = []

    for i in range(0,len(atom_names1b)):
        if (atom_names1b[i] == 'CA') :
            xCA1.append(x1[i])
            yCA1.append(y1[i])
            zCA1.append(z1[i])        
            iCA1.append(i)        

    #Index-Labels
    neg_idx = tf.sparse.to_dense(record['neg_idx'])
    pos_idx = tf.sparse.to_dense(record['pos_idx'])
    pos_idx = np.array(pos_idx)
    neg_idx = np.array(neg_idx)

    index = starting_index  
    for k in range(0,len(pos_idx),2):

        ### For POSITIVE ###################################################### 
        r = np.array( [x0[pos_idx[k]],y0[pos_idx[k]],z0[pos_idx[k]] ])

        rmax0 = r + bl/2
        rmin0 = r - bl/2 

        xList0 = []
        yList0 = []
        zList0 = []
        sList0 = []

        for i in range(0,len(x0)):
            xarr = np.array([x0[i],y0[i],z0[i]])
            if ( (sum(rmax0>xarr) + sum(rmin0<xarr) ) == 6 ): 
                #atom inside current surfacelet
                if ( elements0b[i] == 'C' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(1)
                if ( elements0b[i] == 'N' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(2)    
                if ( elements0b[i] == 'O' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(3)    
                if ( elements0b[i] == 'S' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(4)
                    
        r = np.array( [x1[pos_idx[k+1]],y1[pos_idx[k+1]],z1[pos_idx[k+1]] ])

        rmax1 = r + bl/2
        rmin1 = r - bl/2 

        xList1 = []
        yList1 = []
        zList1 = []
        sList1 = []

        for i in range(0,len(x1)):
            xarr = np.array([x1[i],y1[i],z1[i]])
            if ( (sum(rmax1>xarr) + sum(rmin1<xarr) ) == 6 ): 
                #atom inside current surfacelet
                if ( elements1b[i] == 'C' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(1)
                if ( elements1b[i] == 'N' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(2)    
                if ( elements1b[i] == 'O' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(3)    
                if ( elements1b[i] == 'S' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(4)
                    
                    
        points = np.zeros(np.array([len(xList0),3]) )

        points[:,0] = np.array(xList0)
        points[:,1] = np.array(yList0)
        points[:,2] = np.array(zList0)

        x_y_z = [bl, bl, bl]
        
        sList = np.zeros(np.array([len(xList0),1]))
        sList[:,0] = np.array(sList0)
        
        data = np.concatenate((points,sList), axis=1)
        with open('/home/gxs372/protein/Stanford_Interface/pointdata/p1/'+repr(index)+'.csv','a') as f:
            np.savetxt(f,data, fmt='%s')

        points = np.zeros(np.array([len(xList1),3]) )

        points[:,0] = np.array(xList1)
        points[:,1] = np.array(yList1)
        points[:,2] = np.array(zList1)

        x_y_z = [bl, bl, bl]
        
        sList = np.zeros(np.array([len(xList1),1]))
        sList[:,0] = np.array(sList1)
        
        data = np.concatenate((points,sList), axis=1)
        with open('/home/gxs372/protein/Stanford_Interface/pointdata/p2/'+repr(index)+'.csv','a') as f:
            np.savetxt(f,data, fmt='%s')

        index = index+1    
            
        ### For NEGATIVE ##################################################################### 
        r = np.array( [x0[neg_idx[k]],y0[neg_idx[k]],z0[neg_idx[k]] ])

        rmax0 = r + bl/2
        rmin0 = r - bl/2 

        xList0 = []
        yList0 = []
        zList0 = []
        sList0 = []

        for i in range(0,len(x0)):
            xarr = np.array([x0[i],y0[i],z0[i]])
            if ( (sum(rmax0>xarr) + sum(rmin0<xarr) ) == 6 ): 
                #atom inside current surfacelet
                if ( elements0b[i] == 'C' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(1)
                if ( elements0b[i] == 'N' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(2)    
                if ( elements0b[i] == 'O' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(3)    
                if ( elements0b[i] == 'S' ):
                    xList0.append(xarr[0])
                    yList0.append(xarr[1])
                    zList0.append(xarr[2])
                    sList0.append(4)

        r = np.array( [x1[neg_idx[k+1]],y1[neg_idx[k+1]],z1[neg_idx[k+1]] ])

        rmax1 = r + bl/2
        rmin1 = r - bl/2 

        xList1 = []
        yList1 = []
        zList1 = []
        sList1 = []

        for i in range(0,len(x1)):
            xarr = np.array([x1[i],y1[i],z1[i]])
            if ( (sum(rmax1>xarr) + sum(rmin1<xarr) ) == 6 ): 
                #atom inside current surfacelet
                if ( elements1b[i] == 'C' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(1)
                if ( elements1b[i] == 'N' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(2)    
                if ( elements1b[i] == 'O' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(3)    
                if ( elements1b[i] == 'S' ):
                    xList1.append(xarr[0])
                    yList1.append(xarr[1])
                    zList1.append(xarr[2])
                    sList1.append(4)
                    
        points = np.zeros(np.array([len(xList0),3]) )

        points[:,0] = np.array(xList0)
        points[:,1] = np.array(yList0)
        points[:,2] = np.array(zList0)

        x_y_z = [bl, bl, bl]

        sList = np.zeros(np.array([len(xList0),1]))
        sList[:,0] = np.array(sList0)
        
        data = np.concatenate((points,sList), axis=1)
        with open('/home/gxs372/protein/Stanford_Interface/pointdata/p1/'+repr(index)+'.csv','a') as f:
            np.savetxt(f,data, fmt='%s')

        points = np.zeros(np.array([len(xList1),3]) )

        points[:,0] = np.array(xList1)
        points[:,1] = np.array(yList1)
        points[:,2] = np.array(zList1)

        x_y_z = [bl, bl, bl]

        sList = np.zeros(np.array([len(xList1),1]))
        sList[:,0] = np.array(sList1)
        
        data = np.concatenate((points,sList), axis=1)
        with open('/home/gxs372/protein/Stanford_Interface/pointdata/p2/'+repr(index)+'.csv','a') as f:
            np.savetxt(f,data, fmt='%s')

        index =index+1
        
    return index
        


# In[52]:


filenames =np.array(["/home/gxs372/protein/Stanford_Interface/tf_files/dataset_000.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_001.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_002.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_003.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_004.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_005.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_006.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_007.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_008.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_009.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_010.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_011.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_012.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_013.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_014.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_015.tfrecord",
                     "/home/gxs372/protein/Stanford_Interface/tf_files/dataset_016.tfrecord"])
counter = 0

for filename in filenames:
    raw_protein_dataset = tf.data.TFRecordDataset(filename)

    for record in raw_protein_dataset.take(-1).map(_parse_protein_function):
        newcounter = pointdataFromRecord(record, counter, bl = 20)
        counter =  newcounter     
    print(counter)                 
        
print('final count : ', counter)
                    


# In[53]:


with open('/home/gxs372/protein/Stanford_Interface/counter_11_9.csv','a') as f:
    np.savetxt(f, np.array([counter]), fmt='%s')

