from app import app, db
from itsdangerous import (TimedJSONWebSignatureSerializer
                          as Serializer, BadSignature, SignatureExpired)
from passlib.apps import custom_app_context as pwd_context      
import flask.ext.whooshalchemy as whooshalchemy
import datetime
import pdb    

class User(db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(32), unique=True)
    password = db.Column(db.String(120))    
    mail = db.Column(db.String(120), unique=True) #TODO check lengths
    country = db.Column(db.String(120), index=True) 
    date_joined = db.Column(db.DateTime(timezone=False), default=datetime.datetime.utcnow, index=True)
    active = db.Column(db.Boolean(), index=True)
    last_session = db.Column(db.DateTime(timezone=False), default=datetime.datetime.utcnow, index=True)    
    #excess_queries_number = db.Column(db.Integer)
    #excess_queries_started = db.Column(db.DateTime(timezone=False), index=True)
    projects = db.relationship('Project', backref='user', lazy='dynamic')
    visuals = db.relationship('Visual', backref='user', lazy='dynamic')
    mobile = db.relationship('Mobile', backref='user', lazy='dynamic')
        
    def hash_password(self, password):
        return pwd_context.encrypt(password)

    def verify_password(self, password): #password is plan, self.password is encrypted
        return pwd_context.verify(password, self.password)        
                            
    def is_authenticated(self):
        return True

    def is_active(self):
        return True

    def is_anonymous(self):
        return False

    def get_id(self):
        try:
            return unicode(self.id)  # python 2
        except NameError:
            return str(self.id)  # python 3

    def __repr__(self):
        return '<User %r>' % (self.username)

    def generate_auth_token(self, expiration=600):
        s = Serializer(app.config['SECRET_KEY'], expires_in=expiration)
        return s.dumps({'id': self.id})

    @staticmethod
    def verify_auth_token(token):
        s = Serializer(app.config['SECRET_KEY'])
        try:
            data = s.loads(token)
        except SignatureExpired:
            return None    # valid token, but expired
        except BadSignature:
            return None    # invalid token
        user = User.query.get(data['id'])
        return user

class Project(db.Model):
    __tablename__ = 'projects'
    __searchable__ = ['name']
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(32), index=True)    #unique only per user
    description = db.Column(db.String(512), index=True)    
    creation_date = db.Column(db.DateTime(), default=datetime.datetime.utcnow, index=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))    
    visuals = db.relationship('Visual', backref='project', lazy='dynamic')

    def __repr__(self):
        return '<Project Name %r>' % (self.name) 
                
class Visual(db.Model):
    __tablename__ = 'visuals'
    __searchable__ = ['name']
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(120), index=True) 
    visual_url = db.Column(db.String(256), unique=True) #TODO check size
    visual_meta = db.Column(db.String(1000), index=True) #TODO check size    
    creation_date = db.Column(db.DateTime(), default=datetime.datetime.utcnow, index=True)
    project_id = db.Column(db.Integer, db.ForeignKey('projects.id'))
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))   #TODO is it right ?
    # One2many relationship: on project, many visuals, foreign key always in the many side 
    queries = db.relationship('Query', backref='visual', lazy='dynamic')
    
    def __repr__(self):
        return '<Visual Name %r>' % (self.name)         
            
class Query(db.Model):
    __tablename__ = 'queries'
    id = db.Column(db.Integer, primary_key=True)
    phone_id = db.Column(db.String(100), index=True)
    #phone_os = db.Column(db.String(7), index=True) #android or ios
    latitude = db.Column(db.Float)
    longitude = db.Column(db.Float)
    city = db.Column(db.String(100), index=True)    
    quering_date = db.Column(db.DateTime(), default=datetime.datetime.utcnow, index=True)
    visual_id = db.Column(db.Integer, db.ForeignKey('visuals.id'))     

    def __repr__(self):
        return '<Query id %r>' % (self.id) 

class Mobile(db.Model):
    __tablename__ = 'mobile'
    id = db.Column(db.Integer, primary_key=True)    
    website = db.Column(db.String(256), default="footer", index=True)
    font = db.Column(db.String(256), default="Arial,Arial,Helvetica,sans-serif", index=True)
    font_size = db.Column(db.String(32), default="12", index=True)                                
    bold = db.Column(db.Boolean(), default=False, index=True)
    font_color = db.Column(db.String(7), default="#ffffff", index=True)                                
    splash_color = db.Column(db.String(7), default="#6d3353", index=True)        
    logo_url = db.Column(db.String(256), unique=True)
    icon_url = db.Column(db.String(256), unique=True)
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))    

    def __repr__(self):
        return '<Mobile Id %r>' % (self.id) 

class Usage(db.Model):
    __tablename__ = 'usage'
    # a row every first of the month --> monthly usage
    id = db.Column(db.Integer, primary_key=True)    
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    date = db.Column(db.DateTime(), index=True) # current month 
    plan = db.Column(db.String(32), default="freemium", index=True)        
    plan_expiring_date = db.Column(db.DateTime(), index=True)            
    queries = db.Column(db.Integer)    

class Billing(db.Model):
    __tablename__ = 'billing'
    id = db.Column(db.Integer, primary_key=True)        
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    paiement_date = db.Column(db.DateTime(), index=True)    
    bill_type = db.Column(db.String(5), index=True)    # extra or plan
    plan = db.Column(db.String(35), index=True)
    extra_queries = db.Column(db.Integer) #nb of queries exceeded if bill_type == extra
    total = db.Column(db.Float)        
    state = db.Column(db.Boolean(), index=True)

class Notification(db.Model):
    __tablename__ = 'notification'
    id = db.Column(db.Integer, primary_key=True)        
    user_id = db.Column(db.Integer, db.ForeignKey('users.id'))
    creation_date = db.Column(db.DateTime(), index=True)  
    fresh = db.Column(db.Boolean(), default=True, index=True)
    view_date = db.Column(db.DateTime(), index=True)      
    notif_type = db.Column(db.String(256), index=True) 
    post = db.Column(db.String(256), index=True) 
 
                    
# whooshalchemy search
import sys
if sys.version_info >= (3, 0):
    enable_search = False
else:    
    import flask.ext.whooshalchemy as whooshalchemy
    #whooshalchemy.whoosh_index(app, Visual)
    whooshalchemy.whoosh_index(app, Project)    
    
    
    
    
        
    
