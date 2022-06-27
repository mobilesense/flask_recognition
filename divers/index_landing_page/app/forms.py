from flask.ext.wtf import Form
from wtforms import BooleanField, TextField, PasswordField, SelectField, StringField, IntegerField, validators
from flask.ext.wtf.html5 import URLField
        
class SubscriptionForm(Form):
    mail = TextField('Email address', validators=[
           validators.Required('Please provide a valid email address'),
           validators.Length(min=4, message=(u'Email address too short')),
           validators.Length(max=240, message=(u'Email address too long')),
           validators.Email(message=(u'That\'s not a valid email address.'))
    ])
