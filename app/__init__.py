"""The main routine that starts yarilab """

import flask
from flask import Flask
from flask.ext.login import LoginManager
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.mail import Mail    
#from geopy.geocoders import Nominatim
#from flask_googlelogin import GoogleLogin

#yarilab core api
#import core

# Obtain the flask app object
app = Flask(__name__)
app.config.from_object('app_config')

db = SQLAlchemy(app)

#website login manager
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.session_protection = None

#google login
#googlelogin = GoogleLogin(app, login_manager)

#mail
mail = Mail(app)

#geolocator 
#geolocator = Nominatim()

from app import views, models

# api constructor
from api.manager import Manager, UserManager
import pdb
ApiYarilabManagers = {}
manager = Manager()
for user in models.User.query.all():
    ApiYarilabManagers[user.id] = UserManager(user.id, manager.cached)  
    #except: api_yarilab_managers[user.id] = None



