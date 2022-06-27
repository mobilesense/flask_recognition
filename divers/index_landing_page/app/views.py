import flask
from flask import url_for, request, abort, jsonify, g, session, flash, redirect
from itsdangerous import URLSafeSerializer, BadSignature    
from app import app, db
from .models import Subscription
from .forms import SubscriptionForm
from pytz import timezone

#divers

import datetime
import pdb

@app.route('/', methods=['GET', 'POST'])
def index():   
    if request.method == 'POST':                
        form = SubscriptionForm(request.form)                
        if form.validate():        
            subscription = Subscription()
            form.populate_obj(subscription)            
            mail_exist = Subscription.query.filter_by(mail=form.mail.data).first()            
            if mail_exist:
                form.mail.errors.append('You are already subscribed.')           
                return flask.render_template('index.html', form = form)
            else:                     
                now = datetime.datetime.now(timezone('UTC'))   
                subscription.date = now                              
                db.session.add(subscription)                                                
                db.session.commit() 
                form.mail.errors.append('E-mail saved. You will get notified.')           
                return flask.render_template('index.html', form = form)
        else:
            form.mail.errors.append('Please enter a valid mail')
            return flask.render_template('index.html', form = form)      
    return flask.render_template('index.html', form = SubscriptionForm() )
