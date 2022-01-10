#!python
import tensorflow as tf
#tf.config.gpu.set_per_process_memory_growth(True)
gpus = tf.config.experimental.list_physical_devices('GPU')
for gpu in gpus:
	tf.config.experimental.set_memory_growth(gpu, True)

import sys
from functions import read_data, one_hot_label, roc_file
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2

####################Model building####################
from tensorflow.keras.layers import Dense, Dropout, Activation
from tensorflow.keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from tensorflow.keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.models import Model
from tensorflow.keras.models import load_model, Sequential, model_from_json

model1 = model_from_json(open('../3.model/model.alldata.json').read())
model1.load_weights('../3.model/model.alldata.h5')

max_len = 2000

nlspros = []
with open("../0.database/seqnls.uniprot.nls.txt") as f:
	for line in f:
		data = line.split('\t')
		nlspros.append(data[0])

##NLS.dataset
all_data, all_label, all_line = [], [], []
with open("../0.datadeal/eukaryota.whole.fixedlength.dataset") as f:
	tittle = f.readline().strip('\n')
	for line in f:
		data = line.strip('\n').split('\t')
		if data[0] in nlspros:
			all_line.append(line.strip('\n'))
			all_data.append(data[2])
			all_label.append(data[1])


train_data = all_data
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)

score = model1.predict(train_ohe)
t_label, t_score = [], []
for i in range(len(all_label)):
	if all_label[i] == "NonNucle":
		t_label.append(0)
		t_score.append(score[i][1])
	elif all_label[i] == "Nucleus":
		t_label.append(1)
		t_score.append(score[i][1])
	else:
		continue


print(roc_auc_score(t_label,t_score))
roc_file(t_score, t_label, num="pre_nls")
