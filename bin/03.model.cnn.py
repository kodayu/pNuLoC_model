#!python
from functions import read_data, one_hot_label
import numpy as np
from sklearn.model_selection import train_test_split
from keras.utils import np_utils
from keras.layers.core import Dense, Dropout, Activation
from keras.layers import Masking, Embedding, Input, LSTM, Bidirectional, Conv1D
from keras.layers import Conv1D, Concatenate, MaxPooling1D, BatchNormalization, Add, Flatten
from tensorflow.keras.layers import Attention
from sklearn.preprocessing import OneHotEncoder
from keras import optimizers
from keras.preprocessing import sequence
from keras.callbacks import EarlyStopping
from keras.models import Sequential, Model
from sklearn.metrics import roc_auc_score
from keras.models import Sequential
from keras.regularizers import l1, l2

max_len = 2000
file1 = '../0.datadeal/eukaryota.train.fixedlength.dataset'
train_data, train_label = read_data(file1)
train_data = one_hot_label(train_data)
train_data = sequence.pad_sequences(train_data, maxlen=max_len, padding='post',value=0)
train_data = np.array(train_data)
file2 = '../0.datadeal/eukaryota.valid.fixedlength.dataset'
vali_data, vali_label = read_data(file2)
vali_data = one_hot_label(vali_data)
vali_data = sequence.pad_sequences(vali_data, maxlen=max_len, padding='post',value=0)
vali_data = np.array(vali_data)

# One hot encoding of sequences
train_ohe = np_utils.to_categorical(train_data, num_classes=22)
val_ohe = np_utils.to_categorical(vali_data, num_classes=22)
print(train_ohe.shape, val_ohe.shape)
#all_data = all_data.astype(np.int)
label_train = np.array(train_label, dtype = int)
label_vali = np.array(vali_label, dtype = int)
y_train = np_utils.to_categorical(label_train, num_classes=2)
y_vali = np_utils.to_categorical(label_vali, num_classes=2)

####################Model building####################
##BiLSTM
##embed size, LSTM nodes, Dropout rate
x_input = Input(shape=(2000,))
emb = Embedding(22, 128, input_length=max_len)(x_input)
bi_rnn = Bidirectional(LSTM(64, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(emb)
x = Dropout(0.3)(bi_rnn)
# softmax classifier
x_output = Dense(2, activation='softmax')(x)
model1 = Model(inputs=x_input, outputs=x_output)
model1.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
history1 = model1.fit(train_data, y_train, batch_size=256, epochs=50, validation_data=(vali_data, y_vali),callbacks=cb)

####################MultiCNN####################
def residual_block(data, filters, sizes):
	bn1 = BatchNormalization()(data)
	act1 = Activation('relu')(bn1)
	conv1 = Conv1D(filters, sizes[0], padding='valid', kernel_regularizer=l2(0.001))(act1)
	#bottleneck convolution
	bn2 = BatchNormalization()(conv1)
	act2 = Activation('relu')(bn2)
	conv2 = Conv1D(filters, sizes[1], padding='valid', kernel_regularizer=l2(0.001))(act2)
	#third convolution
	bn3 = BatchNormalization()(conv2)
	act3 = Activation('relu')(bn3)
	conv3 = Conv1D(filters, sizes[2], padding='valid', kernel_regularizer=l2(0.001))(act3)
	#skip connection
	return conv3

x_input = Input(shape=(2000,22))
cnn1 = residual_block(x_input, 128, (4,8,12))
cnn2 = residual_block(x_input, 128, (12,4,8))
cnn3 = residual_block(x_input, 128, (8,12,4))
con123 = Concatenate()([cnn1, cnn2, cnn3])
bn11 = BatchNormalization()(con123)
act11 = Activation('relu')(bn11)
maxpool11 = MaxPooling1D(1979)(act11)
drop = Dropout(0.5)(maxpool11)
# softmax classifier
fla = Flatten()(drop)
dense = Dense(512)(fla)
x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)

model2 = Model(inputs=x_input, outputs=x_output)
model2.summary()
model2.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
history2 = model2.fit(train_ohe, y_train, batch_size=256, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)


####################MultiCNN+BiLSTM####################
def residual_block(data, filters, sizes):
	bn1 = BatchNormalization()(data)
	act1 = Activation('relu')(bn1)
	conv1 = Conv1D(filters, sizes[0], padding='valid', kernel_regularizer=l2(0.001))(act1)
	#bottleneck convolution
	bn2 = BatchNormalization()(conv1)
	act2 = Activation('relu')(bn2)
	conv2 = Conv1D(filters, sizes[1], padding='valid', kernel_regularizer=l2(0.001))(act2)
	#third convolution
	bn3 = BatchNormalization()(conv2)
	act3 = Activation('relu')(bn3)
	conv3 = Conv1D(filters, sizes[2], padding='valid', kernel_regularizer=l2(0.001))(act3)
	#skip connection
	return conv3

x_input = Input(shape=(2000,22))
cnn1 = residual_block(x_input, 128, (4,8,12))
cnn2 = residual_block(x_input, 128, (12,4,8))
cnn3 = residual_block(x_input, 128, (8,12,4))
con123 = Concatenate()([cnn1, cnn2, cnn3])
bn11 = BatchNormalization()(con123)
act11 = Activation('relu')(bn11)
bi_rnn = Bidirectional(LSTM(64, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)

drop = Dropout(0.5)(bi_rnn)
# softmax classifier
dense = Dense(512)(drop)
x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)

model3 = Model(inputs=x_input, outputs=x_output)
model3.summary()
model3.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
history3 = model3.fit(train_ohe, y_train, batch_size=256, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)


####################MultiCNN+BiLSTM+Attention####################
def residual_block(data, filters, sizes):
	bn1 = BatchNormalization()(data)
	act1 = Activation('relu')(bn1)
	conv1 = Conv1D(filters, sizes[0], padding='valid', kernel_regularizer=l2(0.001))(act1)
	#bottleneck convolution
	bn2 = BatchNormalization()(conv1)
	act2 = Activation('relu')(bn2)
	conv2 = Conv1D(filters, sizes[1], padding='valid', kernel_regularizer=l2(0.001))(act2)
	#third convolution
	bn3 = BatchNormalization()(conv2)
	act3 = Activation('relu')(bn3)
	conv3 = Conv1D(filters, sizes[2], padding='valid', kernel_regularizer=l2(0.001))(act3)
	#skip connection
	return conv3

x_input = Input(shape=(2000,22))
cnn1 = residual_block(x_input, 128, (4,8,12))
cnn2 = residual_block(x_input, 128, (12,4,8))
cnn3 = residual_block(x_input, 128, (8,12,4))
con123 = Concatenate()([cnn1, cnn2, cnn3])
bn11 = BatchNormalization()(con123)
act11 = Activation('relu')(bn11)
bi_rnn = Bidirectional(LSTM(64, kernel_regularizer=l2(0.01), recurrent_regularizer=l2(0.01), bias_regularizer=l2(0.01)))(act11)
drop = Dropout(0.5)(bi_rnn)
# softmax classifier
atten = Attention()([drop, drop])
dense = Dense(512)(atten)
x_output = Dense(2, activation='softmax', kernel_regularizer=l2(0.0001))(dense)

model4 = Model(inputs=x_input, outputs=x_output)
model4.summary()
model4.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
cb = [EarlyStopping(monitor='val_accuracy',patience=2)]
history4 = model4.fit(train_ohe, y_train, batch_size=256, epochs=50, validation_data=(val_ohe, y_vali),callbacks=cb)



score = model2.predict(train_ohe)
t_label = [int(i) for i in train_label]
t_score = [i[1] for i in score]
print(roc_auc_score(t_label,t_score))

score1 = model2.predict(val_ohe)
t_label1 = [int(i) for i in vali_label]
t_score1 = [i[1] for i in score1]
print(roc_auc_score(t_label1,t_score1))


