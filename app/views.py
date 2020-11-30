from django.shortcuts import render
from .forms import SeqForm

import numpy as np
import tensorflow as tf
from keras.models import load_model

def one_hot_encode(seq, base_map):
    seq = seq.upper()
    mapping = dict(zip(base_map, range(4))) 
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

# Create your views here.
one_bp_model = load_model('models/1bpins_manual60.h5')
delfreq_model = load_model('models/delfreq_manual60.h5')

def get_input_view(request):
    if request.method == 'POST':
        form = SeqForm(request.POST)
        context = {'form': form}
        if form.is_valid(): 
            input_seq = form['input_seq'].value()
            input_seq = one_hot_encode(input_seq, 'ACGT')
            input_seq = np.reshape(input_seq, (1, 60, 4))
            
            one_bp_pred = one_bp_model.predict(input_seq)
            one_bp_pred = one_bp_pred.flatten().tolist()[0] * 100
            one_bp_pred = str(round(one_bp_pred, 2)) + ' %'

            delfreq_pred = delfreq_model.predict(input_seq)
            delfreq_pred = delfreq_pred.flatten().tolist()[0] * 100
            delfreq_pred = str(round(delfreq_pred, 2)) + ' %'

            context = {
                'form': SeqForm(),
                'input_seq': form['input_seq'].value(),
                '1bp_pred': one_bp_pred,
                'delfreq_pred': delfreq_pred
            }
    
    else: # if a GET (or any other method), create a blank form
        context = {'form': SeqForm()}
    
    return render(request, 'seqform.html', context)