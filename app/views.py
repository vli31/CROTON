from django.shortcuts import render
from .forms import SeqForm

import numpy as np
import tensorflow as tf
from tensorflow.keras.models import load_model

def one_hot_encode(seq, base_map):
    seq = seq.upper()
    mapping = dict(zip(base_map, range(4))) 
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

# Create your views here.

def get_input_view(request):
    if request.method == 'POST':
        form = SeqForm(request.POST)
        context = {'form': form}
        if form.is_valid(): 
            input_seq = form['input_seq'].value()
            input_seq = one_hot_encode(input_seq, 'ACGT')
            input_seq = np.reshape(input_seq, (1, 60, 4))
            model = load_model('models/CROTON.h5') #load multitask model
            pred = model.predict(input_seq)
            
            delfreq = pred[:,0].flatten().tolist()[0] * 100
            prob_1bpins = pred[:,1].flatten().tolist()[0] * 100
            prob_1bpdel = pred[:,2].flatten().tolist()[0] * 100
            onemod3_freq = pred[:,3].flatten().tolist()[0] * 100
            twomod3_freq = pred[:,4].flatten().tolist()[0] * 100
            frameshift_freq = pred[:,5].flatten().tolist()[0] * 100

            delfreq = str(round(delfreq, 2)) + ' %'
            prob_1bpins = str(round(prob_1bpins, 2)) + ' %'
            prob_1bpdel = str(round(prob_1bpdel, 2)) + ' %'
            onemod3_freq = str(round(onemod3_freq, 2)) + ' %'
            twomod3_freq = str(round(twomod3_freq, 2)) + ' %'
            frameshift_freq = str(round(frameshift_freq, 2)) + ' %'
            
            context = {
                'form': SeqForm(),
                'input_seq': form['input_seq'].value(),
                'prob_1bpins': prob_1bpins, 
                'prob_1bpdel': prob_1bpdel, 
                'delfreq': delfreq, 
                'onemod3_freq': onemod3_freq, 
                'twomod3_freq': twomod3_freq, 
                'frameshift_freq': frameshift_freq,
            }

    else: # if GET (or any other method), create a blank form
        context = {'form': SeqForm()}
    
    return render(request, 'seqform.html', context)