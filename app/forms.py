from flask.ext.wtf import Form
from wtforms import BooleanField, TextField, PasswordField, SelectField, StringField, IntegerField, validators
from flask.ext.wtf.html5 import URLField
from app.globals import country_choices, font_choices

class LoginForm(Form):
    username = TextField('Username', [
        validators.Length(min=6, message=(u'Username too short. Min. 6 characters')),
        validators.Length(max=32, message=(u'Username too long. Max. 32 characters'))
    ])    
    password = PasswordField('Password', validators=[
        validators.Length(min=6, message=(u'Password too short. Min. 6 characters')),
        validators.Length(max=64, message=(u'Password too long. Max. 64 characters')),
        validators.Required()
    ])
    remember_me = BooleanField('remember_me', default=False)    

class ResetForm(Form):
    mail = TextField('Email address', validators=[
           validators.Required('Please provide a valid email address'),
           validators.Length(min=4, message=(u'Email address too short')),
           validators.Length(max=240, message=(u'Email address too long')),
           validators.Email(message=(u'That\'s not a valid email address.'))
           ])     

class EditUserInfoForm(Form):
    username = TextField('Username', validators=[
           validators.Required('Please provide a valid username'),
           validators.Length(min=6, message=(u'username too short. min. 6 characters')),
           validators.Length(max=32, message=(u'username too long. max. 32 characters'))
           ])
    mail = TextField('Email address', validators=[
           validators.Required('Please provide a valid email address'),
           validators.Length(min=4, message=(u'Email address too short')),
           validators.Length(max=240, message=(u'Email address too long')),
           validators.Email(message=(u'That\'s not a valid email address.'))
           ])
    country = SelectField('Country', validators=[
           validators.Length(min=2, message=(u'Please select your country'))],
           choices = country_choices)  
    
class NewProjectForm(Form):
    name = TextField('Project name', validators=[
           validators.Required('Please provide a valid project name'),
           validators.Length(min=3, message=(u'name too short. min. 3 characters')),
           validators.Length(max=32, message=(u'name too long. max. 32 characters'))
           ])
    description = TextField('Project description', validators=[                      
           validators.Length(max=512, message=(u'username too long. max. 32 characters'))
           ])

class EditProjectForm(Form):
    name = TextField('Project name', validators=[
           validators.Required('Please provide a valid project name'),
           validators.Length(min=3, message=(u'name too short. min. 3 characters')),
           validators.Length(max=32, message=(u'name too long. max. 32 characters'))
           ])
    description = TextField('Project description', validators=[                      
           validators.Length(max=512, message=(u'username too long. max. 512 characters'))
           ])

class EditVisualForm(Form):
    name = TextField('Visual name', validators=[
           validators.Required('Please provide a valid visual name'),
           validators.Length(min=3, message=(u'name too short. min. 3 characters')),
           validators.Length(max=32, message=(u'name too long. max. 32 characters'))
           ])
    visual_meta = URLField('Url meta', validators=[
           validators.Required('Please provide a valid meta url'), 
           validators.URL()
           ])   

class SearchForm(Form):
    search = StringField('search', validators=[validators.DataRequired()])
    #http://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-x-full-text-search

                                                
class ChangePassForm(Form):
    old_password = PasswordField('Old Password', validators=[
        validators.Length(min=6, message=(u'Please give a longer password. min. 6 characters')),
        validators.Length(max=240, message=(u'Please give a shorter password (max. 240 characters)')),
        validators.Required()        
    ])
    new_password = PasswordField('New Password', validators=[
        validators.Length(min=6, message=(u'Please give a longer password. min. 6 characters')),
        validators.Length(max=240, message=(u'Please give a shorter password (max. 240 characters)')),
        validators.Required(),
        validators.EqualTo('new_password_confirm', message='Passwords must match')
    ])    
    new_password_confirm = PasswordField('Repeat Password')
    
            
class RegistrationForm(Form):
    username = TextField('Username', validators=[
           validators.Required('Please provide a valid username'),
           validators.Length(min=6, message=(u'username too short. min. 6 characters')),
           validators.Length(max=32, message=(u'username too long. max. 32 characters'))
           ])
    mail = TextField('Email address', validators=[
           validators.Required('Please provide a valid email address'),
           validators.Length(min=4, message=(u'Email address too short')),
           validators.Length(max=240, message=(u'Email address too long')),
           validators.Email(message=(u'That\'s not a valid email address.'))
           ])
    country = SelectField('Country', validators=[
           validators.Length(min=2, message=(u'Please select your country'))],
           choices = country_choices)       
    plan = TextField('Plan', validators=[])
    password = PasswordField('New Password', validators=[
        validators.Length(min=6, message=(u'Please give a longer password. min. 6 characters')),
        validators.Length(max=64, message=(u'Please give a shorter password (max. 64 characters)')),
        validators.Required(),
        validators.EqualTo('password_confirm', message='Passwords must match')
    ])
    password_confirm = PasswordField('Repeat Password')
    agree = BooleanField('I agree all your Terms of Services',
           validators=[validators.Required(u'You must accept our Terms of Service')])

class EditMobileForm(Form):
    splash_color = TextField('Splash color', validators=[
           validators.Length(min=4, message=(u'Wrong string. This have to be a hex code (#xxxxxx or #xxx)')),
           validators.Length(max=7, message=(u'Wrong string. This have to be a hex code (#xxxxxx or #xxx)'))           
           ])  
    website = TextField('Website address', validators=[           
           validators.Length(max=240, message=(u'Website address too long'))           
           ])
    font = TextField('font', validators=[
            validators.Required()
            ])      
    font_size = IntegerField('Font size', [validators.NumberRange(min=5, max=50)])           
    bold = BooleanField('bold', default=False)    
    font_color = TextField('Font color', validators=[
           validators.Length(min=4, message=(u'Wrong string. This have to be a hex code (#xxxxxx or #xxx)')),
           validators.Length(max=7, message=(u'Wrong string. This have to be a hex code (#xxxxxx or #xxx)'))           
           ])
    logo_url = TextField('Logo url', validators=[])
    icon_url = TextField('Icon url', validators=[])  
                                      
