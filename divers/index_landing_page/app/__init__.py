"""The main routine that starts yarilab """

import flask
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy
from flask.ext.mail import Mail    
from flask.ext.storage import get_default_storage_class
from demo_api.manager import engine

app = Flask(__name__)
app.config.from_object('config')
db = SQLAlchemy(app)
engine = engine('config_yari:DemoYarik20k(%d)'%db_id)  

from app import views, models
