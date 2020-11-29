from django import forms
from django.core.exceptions import ValidationError

def validate_base(input_seq):
    errors = []
    base_lst = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']
    if any(x not in base_lst for x in input_seq):
        errors.append(ValidationError(
            "Error: The input should only contain A, T, C, and G"
        ))
    
    if len(input_seq) != 79:
        errors.append(ValidationError(
            "Error: The input has %i (not 79) characters" % (len(input_seq))
        ))

    if errors:
        raise ValidationError(errors)

class SeqForm(forms.Form):
    input_seq = forms.CharField(label='',  validators=[validate_base], 
        widget=forms.Textarea(attrs={'cols': 50, 'placeholder': 'Your sequence...'}),
        error_messages={'required': 'Error: Nothing was inputted'})